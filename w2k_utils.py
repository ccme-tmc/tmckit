#!/usr/bin/env python
import sys,os,shutil,subprocess,commands,string,glob
from math import *
from latt_utils import *
from chem_utils import *
from constants  import *
from io_utils   import *
import gparams
#import w2k_caseutil

from common_caseutil import f_GetCaseName
from w2k_caseutil import Wien2K_File_in0,Wien2K_Structure,Wien2K_Case

w2k_save_cmd = "w2k_save -F -s"

def W2k_refine_struct(case_name,rmt=5.0,sgroup=False,out=None):
  """
  refine the struct file starting from a raw struct
  """
  # first check whether R_MT's are properly set
  ierr = W2k_run("x_lapw -f "+case_name+" nn",inp="2",out=out)
  if ierr != 0:
    print "WARNING: R_MT's are not set properly"
    ierr = W2k_run( "setrmt_lapw "+case_name + " -r %f"%(rmt), out=out)
    if ierr != 0:
      print "ERROR: something wrong ! when running setrmt"
      return -1
    else:
      shutil.copy(case_name+".struct_setrmt",case_name+".struct")

  if sgroup:
    # use sgroup to generate a new struct file
    ierr = W2k_run( "x_lapw -f " + case_name + " sgroup", out=out )
    if ierr != 0:
      print "ERROR: something wrong ! when running sgroup"
      return -2
    else:
      shutil.copy(case_name+".struct_sgroup",case_name+".struct")
  return ierr

def W2k_read_energy(fn,nat,debug=True,para=False):

  """
  Read band energies in the wien2k case.energy format
  """
  if para is False:
    fn_all = [fn]
  else:
    fn_all = glob.glob(fn+"_[1-9]*")

  print "Read band energies from ",fn_all

  enk_all = []
  kw_all = []
  kv_all = []
  ik=0
  n_files = len(fn_all)

  for i in range( n_files ):
    if n_files == 1:
      fname = fn
    else:
      fname = fn+"_%d"%(i+1)

    ifile = open(fname,'r')

    # skip the heads related to the linearization energy
    f_Skip_Lines(ifile,nat*2)
    while 1:
      line = ifile.readline()
      if not line: break # end of file readed

      ik += 1
      # k-vectors, number of bands at each k, and the weight of this k
      kv = [ float(line[0:19]), float(line[19:38]), float(line[38:57]) ]
      nbk = int(line[73:79])
      kw =  float(line[79:84])

      if debug:
        print "Read the band energy for ik=%5d,k=(%8.4f, %8.4f, %8.4f), nbk=%5d, wk=%5.1f"%(ik,kv[0],kv[1],kv[2],nbk,kw)
        print "%5s %12s"%("n","Enk")

      enk = []
      for ie in range(nbk):
        s_en = ifile.readline().split()
        n = int(s_en[0])
        en = float(s_en[1])
        enk.append(en)
        if debug: print "%5d %12.6f"%(n,en)

      kv_all.append(kv)
      kw_all.append(kw)
      enk_all.append(enk)
    # end of reading a file
  # end of read all files
  return enk_all,kw_all,kv_all

def W2k_run(w2k_cmd,inp=None,out=None ):
  """
  Run a wien2k job
  """
  if not inp is None:
    tmp_sh = "w2k-tmp.sh"

    otmp = open(tmp_sh,'w')
    otmp.write(w2k_cmd + " <<EOF\n")
    tmp = inp.split(';')
    for i in range(len(tmp)):
      otmp.write(tmp[i]+"\n")
    otmp.write("EOF\n")
    otmp.close()
    job_cmd = "bash " + tmp_sh
  else:
    job_cmd = w2k_cmd

  if not out is None:
    job_cmd = job_cmd + " 1>> " + out +" 2>&1 "

  ierr = os.system(job_cmd)

  return ierr

def W2k_run_batch(w2k_cmd,out=None ):
  """
  Run a series w2k jobs
  """
  jobs = w2k_cmd.split(";")
  for i in range( len(jobs) ):
    job_cmd = jobs[i]
    print "  Run the job ",job_cmd
    if not out is None:
      job_cmd = job_cmd + " 1>> " + out +" 2>&1 "

    ierr = os.system(job_cmd)
    if ierr != 1: break

  return ierr


def w2k_get(case_name,flag):
  """
  Extract parameters from wien2k input/output files
  """

  if flag == "nat":
    return f_Get_Natom(case_name,0)
  elif flag == "nat_all":
    return f_Get_Natom(case_name,1)
  elif flag == "etot":
    return f_Read_Etot(case_name)
  elif flag == "vol":
    return f_Read_Vol(case_name)
  elif flag == "efer":
    return f_Read_Efer(case_name,0)
  elif flag == 'egap':
    return w2k_get_egap(case_name)
  else:
    print "WARNING in w2k_get: Unsupported flag "+flag
    return None

def w2k_get_nproc(case_name,sp):
  """
  Get the number of processed used in k-parallel calculations
  """

  if sp==0:
    sys_cmd = "ls "+case_name+".scf1_* | wc -w"
  else:
    sys_cmd = "ls "+case_name+".scf1up_* | wc -w"

  failure,output = commands.getstatusoutput(sys_cmd)
  if failure:
    print "WARNING: fail to run " + sys_cmd
    print "set nproc=0"
    nproc = 0
  else:
    nproc = int(output.split()[0])

  return nproc


def f_w2k_atom_index(case_name,atom,debug=True):
  """
  Get the indices for a particular atom
  """
  (latt_type,latt, basis_eq) = f_Read_Struct_w2k(case_name,mode=1)
  nat = len(basis_eq)

  ind_at=[]
  for iat in range(nat):
    if basis_eq[iat][0][0:2] == atom[0:2]:
      ind_at.append(iat+1)
  return ind_at

def w2k_occlsj(ne,l,ik):
  """
  set the occupation number for kappa = l and kappa = -l - 1 for (ik= 1 or 2)
  it is assumed that the kappa=l state is first occupied before kappa = -l -1
  """
  if ik==1:
    occ=max(min(ne,2.0*l),0.0001)
  else:
    occ=max(0.0,ne-2.0*l)
  return occ

def w2k_set_opencore(case_name,iat_oc,np_oc,l_oc,nel,sh_inc,tag):
  """
  set up an open-core case.inc
  """
  nat   = w2k_get(case_name,'nat')
  ifile = open(case_name+".inc",'r')
  ofile = open(case_name+".inc"+tag,'w')

  error = False
  for iat in range(nat):
    line = ifile.readline().split()
    norb = int(line[0])
    sh_old = float(line[1])

    if iat == iat_oc:  # the impurity atom
      ofile.write("%2d%5.2f    NUMBER OF ORBITALS (EXCLUDING SPIN), SHIFT\n"%(norb+2,sh_inc))
      for iorb in range(norb):
        line = ifile.readline()
        n = int(line[0:1])
        l = int(line[2:4])

        if n == np_oc and l == l_oc:
          error = True
          break
        ofile.write(line)

      if l_oc == 0:
        ofile.write("%1d,%2d,%7.5      (N,KAPPA,OCCUP)\n"%(np_oc,-1,occ) )
      else:
        ofile.write("%1d,%2d,%7.5      (N,KAPPA,OCCUP)\n"%(np_oc, l_oc,  w2k_occlsj(nel,l_oc,1) ) )
        ofile.write("%1d,%2d,%7.5      (N,KAPPA,OCCUP)\n"%(np_oc,-l_oc-1,w2k_occlsj(nel,l_oc,2) ) )
    else:   # other atoms
      ofile.write("%2d%5.2f    NUMBER OF ORBITALS (EXCLUDING SPIN), SHIFT\n"%(norb,  sh_old))
      for iorb in range(norb):
        line = ifile.readline()
        ofile.write(line)

  ifile.close()
  ofile.close()
  if error:
    print "ERROR: the target orbital is already included in the core!!!"
    sys.exit(1)

def w2k_band_analysis(enk,efer,kvec,emin=-2.0,emax=2.0,ib0=1,debug=False):
  print \
  """
  Analyze band energies
  """
  nsp = len(enk)
  nk  = len(kvec)
  nband = len(enk[0][0])
  print "\tNumber of spin=",nsp
  print "\tNumber of k-vector=",nk
  print "\tNumber of bands=",nband
  print "\n"

  nomax = []
  numin = []
  ikvm =  []
  ikcm =  []

  lmetal = False
  for isp in range(nsp):
    nomax.append(0)
    numin.append(nband)
    ikvm.append(0)
    ikcm.append(0)
    for ik in range(nk):
      nbk = len(enk[isp][ik])

      if debug:
        print "  %5d bands for ik=%5d"%(nbk,ik)

      if nbk < nband: nband = nbk
      nv = 0
      nc = nbk-1
      for ib in range(nbk):
        if enk[isp][ik][ib] <= efer + 0.00001:
          if ib > nv: nv = ib
        else:
          if ib < nc: nc = ib

      if nv > nomax[isp]:  nomax[isp] = nv
      if nc < numin[isp]:  numin[isp] = nc

      nv  = nomax[isp]
      nc  = numin[isp]
      ikv = ikvm[isp]
      ikc = ikcm[isp]
      if enk[isp][ik][nv] > enk[isp][ikv][nv]: ikvm[isp] = ik
      if enk[isp][ik][nc] < enk[isp][ikc][nc]: ikcm[isp] = ik
    # loop over ik

    if nomax[isp] >= numin[isp]: lmetal = True

  # end loop over isp

  # set evbm, which is used as the energy zero
  if lmetal :
    evbm = efer
  else:
    if nsp==1:
      evbm = enk[0][ikvm[0]][nomax[0]]
    else:
      evbm = max( enk[0][ikvm[0]][nomax[0]], enk[1][ikvm[1]][nomax[1]] )

  for isp in range(nsp):
    nv = nomax[isp]
    nc = numin[isp]
    ikv = ikvm[isp]
    ikc = ikcm[isp]
    if nsp ==2: print "\nBand analysis for spin",isp
    print "Band index for VBM and CBM=",nv+ib0,nc+ib0

    if not lmetal:
      print "Insulating system:"
      egap = (enk[isp][ikc][nc] - enk[isp][ikv][nv])*Ry2eV
      if ikvm[isp] == ikcm[isp]:
        print ":BandGap(d) =%8.3f eV "%(egap)
        print "  direct gap at k=(%8.3f,%8.3f,%8.3f)"%(kvec[ikv][0],kvec[ikv][1],kvec[ikv][2])
      else:
        print ":BandGap(i) =%8.3f eV "%(egap)
        egvbm = (enk[isp][ikv][nc] - enk[isp][ikv][nv] )*Ry2eV
        egcbm = (enk[isp][ikc][nc] - enk[isp][ikc][nv] )*Ry2eV

        egm = (enk[isp][0][nc] - enk[isp][0][nv])*Ry2eV
        ikm = 0
        for ik in range(nk):
          egk = (enk[isp][ik][nc] - enk[isp][ik][nv])*Ry2eV
          if egk < egm:
            egm = egk
            ikm = ik

        print ":Eg_d(min) =%8.3f eV, at      k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egm,  kvec[ikm][0], kvec[ikm][1], kvec[ikm][2],ikm+1)
        print ":Eg_direct =%8.3f eV  at  VBM k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egvbm,kvec[ikv][0], kvec[ikv][1], kvec[ikv][2],ikv+1)
        print ":Eg_direct =%8.3f eV  at  CBM k=(%8.3f,%8.3f,%8.3f) (ik=%4d)"%(egcbm,kvec[ikc][0], kvec[ikc][1], kvec[ikc][2],ikc+1)
    else:
      print "Metallic system"

    print "Range of each band with respect to VBM (eV):"
    print "%5s%12s%12s%12s"%('n','Bottom','Top','Width')
    for ib in range(nband):

      ebmin=enk[isp][0][ib]
      ebmax=enk[isp][0][ib]

      for ik in range(nk):
        if enk[isp][ik][ib] < ebmin:
          ebmin = enk[isp][ik][ib]
        if enk[isp][ik][ib] > ebmax:
          ebmax = enk[isp][ik][ib]

      ebmin = ebmin - evbm
      ebmax = ebmax - evbm

      if ebmin > emin  and ebmax < emax:
        print "%5d%12.3f%12.3f%12.3f"%(ib+ib0,ebmin*Ry2eV,ebmax*Ry2eV,(ebmax-ebmin)*Ry2eV)
      # end loop over ib
   # end loop over isp

  return

def w2k_get_egap(case_name,iop=0):
  """
  Get the band gap
  """
  fn_scf2 = case_name+".scf2"
  fn_scf2dos = case_name+".scf2-dos"

  read_egap = io_grep_lines( fn_scf2, ':GAP',-1,  6)
  if read_egap is None:
    print "!!! Metallic: check to ensure this is right!"
    egap = 0.0
  else:
    egap = float( read_egap )
    if egap < 0.0:
      print "WARNING: negative gap is read from "+fn_scf2
      print " -- set the gap to be 0.0"
      egap = 0.0

  return egap

def w2k_get_ecore(casename,iat_core,nl_symbl,sp_tag='',debug=False):
  """
  Get the energy of the core state, denoted as "nl_symbl" on the atom "iat"
  """
  f_name = "w2k_get_ecore"
  ifile = open(casename+".scfc"+sp_tag,'r')

  iat = 0
  found_none = False
  while (1):
    line = ifile.readline()

    if not line :
      found_none = True
      break

    if "CORE STATES" in line: # the first line for each atom
      iat += 1
      if iat > iat_core:
        found_none = True
        break

      line_str = line.split()
      atom = line_str[1]
      ncore = int(line[41:43])

      ic = 0
      while ic < ncore:
        line = ifile.readline()
        if 'CORE-FORCE' in line: # skip the lines related to the force
          for i in range(3): line = ifile.readline()
        ic += 1
        if iat == iat_core and line[1:3].strip() == nl_symbl.strip() :
          ecore = float(line[20:40])
          return ecore

  if found_none:
    print f_name + ": Fail to find the core state " + nl_symbl + " on %s-th atom"%(iat_core)
    sys.exit(1)

def f_Read_Etot(name):
  """
  Read total energy from case.scf
  """
  case_name = f_Check_Name(name)
  grep_cmd = "grep \":ENE  :\" " + case_name.strip() + ".scf | tail -n 1"
  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    sys.exit(1)
  return float(output.split()[8])


def f_Read_Efer(name,iop,sp=0,tag=''):
  """
  Read Fermi energy from case.scf/qtl(iop=0/1)
  """
  case_name = f_Check_Name(name)
  nsp =f_Check_Spin(case_name)

  if sp == 0:
    sptag=''
  else:
    sptag = 'up'

  if iop == 1:
    file=case_name.strip()+".scf2"+sptag+tag
    efer = float( io_grep_lines(file,":FER  :",1,10) )
  elif iop == 2:
    file=case_name.strip()+".qtl"+sptag+tag
    efer = float( io_grep_lines(file,"FERMI ENERGY",1,8) )
  else:
    file=case_name.strip()+".scf"+tag
    efer = float( io_grep_lines(file,":FER  :",-1,10) )

  return efer

def f_Read_Vol(name):
  """
  Read volume from case.scf
  """
  case_name = f_Check_Name(name)
  grep_cmd = "grep \":VOL  :\" " + case_name.strip() + ".scf | tail -n 1"
  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    sys.exit(1)

  return float(output.split()[6])

def f_Read_Energy(fname,nat,debug=True):
  """
  Read band energies in the wien2k case.energy format
  """
  print "Read band energies from " + fname
  ifile = open(fname,'r')

  # skip the heads related to the linearization energy
  f_Skip_Lines(ifile,nat*2)
  ik=0
  enk_all = []
  kw_all = []
  kv_all = []
  while 1:
    line = ifile.readline()
    if not line: break # end of file readed

    ik += 1
    # k-vectors, number of bands at each k, and the weight of this k
    kv = [ float(line[0:19]), float(line[19:38]), float(line[38:57]) ]
    nbk = int(line[73:79])
    kw =  float(line[79:84])

    if debug:
      print "Read the band energy for ik=%5d,k=(%8.4f, %8.4f, %8.4f), nbk=%5d, wk=%5.1f"%(ik,kv[0],kv[1],kv[2],nbk,kw)
      print "%5s %12s"%("n","Enk")

    enk = []
    for ie in range(nbk):
      s_en = ifile.readline().split()
      n = int(s_en[0])
      en = float(s_en[1])
      enk.append(en)
      if debug: print "%5d %12.6f"%(n,en)

    kv_all.append(kv)
    kw_all.append(kw)
    enk_all.append(enk)

  return enk_all,kw_all,kv_all

def f_Write_Log(info):
  """
  Write information to the default log file
  """
  ofile = open(gparams.pyw2k_logfile,'a')
  ofile.write(info+'\n')
  ofile.close()

def f_Check_Name(name=''):
  """
  Check whether <name> is a valid case name for the wien2k struct, i.e. check <name>.struct exists.
  If it is a null string, get the case name from the name of the current directory.
  the valid case name is returned, or if exit with an error message
  """
  if name == '' :
    cwd=os.getcwd()
    case_name = os.path.basename(cwd)
  else:
    case_name = name

  struct_file = case_name+".struct"

  if not os.path.isfile(struct_file) :
    print "ERROR: struct file " + struct_file + "not existing "
    sys.exit(1)
  return case_name

def f_Check_Complex(name,iop=-1):
  """
  Check whether Kohn-Sham vectores are complex, i.e. whether the system has the inversion symmetry
  """
  if iop == 0:
    status, output = commands.getstatusoutput("grep INVERSION "+name.strip()+".outputs")
    nw=len(output.split())
    if nw == 7: cmplx='c'
    else: cmplx=''
  else:
    if os.path.isfile(name+".in1c") and not os.path.isfile(name+".in1"):
      cmplx='c'
    else:
      cmplx=''

  return cmplx

def f_Check_Spin(case_name):
  """
  Check whether the current calculation is spin polarized or with spin orbit couplint (soc)
   return
      0 -- spin unpolarized
      2 -- spin polarized
      4 -- with SOC
  """
  nsp=0
  if os.path.isfile(case_name+"scfso"):
    nsp = 4
  elif not os.path.isfile(case_name+".scf1") and os.path.isfile(case_name+".scf1up") :
    nsp = 2
  else:
    nsp = 0
  return nsp


def f_Read_Latt(input_file):
  """
  Read lattice constant from the struct file
  """
  try:
    ifile = open(input_file,'r')
  except:
    print "ERROR: Fail to open " + input_file
    sys.exit(1)

  f_Skip_Lines(ifile,3)
  line = ifile.readline()
  a = float( line[0:10]  )
  b = float( line[10:20])
  c = float( line[20:30])
  alf = float( line[30:40])
  bet = float( line[40:50])
  gam = float( line[50:60])
  return [a,b,c,alf,bet,gam]

def f_Save_Min(case_name,save_dir):
  """
  This funciton saves important files after min_lapw
  """
  call ("mkdir -p " + save_dir, shell = True)
  call ("mv "+ case_name+"_[0-9]*.clm* " + save_dir, shell=True )
  call ("mv "+ case_name+"_[0-9]*.scf " +save_dir, shell=True )
  call ("mv "+ case_name+"_[0-9]*.struct " +save_dir, shell=True )

  call ("mv "+ case_name+".scf_mini " + save_dir, shell = True)
  call ("cp "+ case_name+".inM " +save_dir, shell=True )
  call ("cp "+ case_name+".scf " + save_dir, shell = True)
  call ("cp "+ case_name+".struct " + save_dir, shell = True)
  call ("cp "+ case_name+".clmsum " + save_dir, shell = True)
  call ("cp "+ case_name+".outputM " + save_dir, shell = True)
  call ("w2k_clean -p", shell = True)


def f_Set_emax(case_name,emax=2.0,complex=''):
  print \
  """
  Set emax in case.in1
  """

  in1_file = case_name+".in1"+complex
  in1_bak = in1_file.strip()+'-bak'
  print "in1 file:",in1_file

  if os.path.isfile(in1_bak): os.remove(in1_bak)
  os.rename(in1_file,in1_bak)

  ifile = open(in1_bak, 'r')
  ofile = open(in1_file,'w')
  lines = ifile.readlines()
  nl = len(lines)
  for il in range(nl-1):
    ofile.write(lines[il])

  line=lines[nl-1]

  ofile.write("%31s%10.1f%s\n"%(line[0:31],emax,line[41:]))

  ifile.close()
  ofile.close()

def f_Set_XCFunc(case_name,xcfunc_i=''):
  """
  Set the exchange correlation (xc) functional in case.in0
  """

  xcfunc = xcfunc_i.lower()
  if   xcfunc == '':        return
  elif xcfunc == 'lda':     ixc = 5
  elif xcfunc == 'pbe':     ixc = 13
  elif xcfunc == 'pbesol':  ixc = 19
  elif xcfunc == 'wc06':    ixc = 11
  elif xcfunc == 'am05':    ixc = 20
  elif xcfunc == 'pw91':    ixc = 14
  elif xcfunc == 'rpbe':    ixc = 16
  elif xcfunc == 'tpss':    ixc = 27
  else:                     ixc = int(xcfunc)

  in0_file = case_name.strip()+".in0"
  in0_bak = case_name.strip()+".in0-bak"
  if os.path.isfile(in0_bak) : os.remove(in0_bak)
  os.rename(in0_file,in0_bak)

  ifile = open(in0_bak,'r')
  ofile = open(in0_file,'w')

  il=0
  for line in ifile:
    il += 1
    if il ==1:
      ll = line.split()
      ll[1] = str(ixc)
      line =''
      for i in range(len(ll)): line = line + ll[i] +"\t"
      line=line+"\n"

    ofile.write(line)

  ifile.close()
  ofile.close()

def f_Set_RKmax(case_name,rk,complex=''):
  in1_file = case_name+".in1"+complex
  in1_bak = in1_file.strip()+'-bak'
  if os.path.isfile(in1_bak): os.remove(in1_bak)
  os.rename(in1_file,in1_bak)

  ifile = open(in1_bak, 'r')
  ofile = open(in1_file,'w')
  il = 0
  for line in ifile:
    il += 1
    if il == 2:
      ll = line.split()
      rk_old = float(ll[0])
      print " Reset RKmax froma ",rk_old, " to ",rk
      ll[0] = str(rk)
      line = ''
      for i in range(len(ll)): line = line + ll[i] + ' '
      line = line + '\n'

    ofile.write(line)
  ifile.close()
  ofile.close()

def f_Set_orbUJ(case_name,iatorb,u=None,j=None):
  """
  Reset the value of U on iatorb-th atom, if zero, reset U on all atoms
  """
  inorb_file = case_name+".inorb"
  inorb_bak = inorb_file.strip()+'-bak'

  print "Make a copy of the old inorb to "+ inorb_bak
  if os.path.isfile(inorb_bak): os.remove(inorb_bak)
  os.rename(inorb_file,inorb_bak)

  ifile = open(inorb_bak, 'r')
  ofile = open(inorb_file,'w')

  il=0
  natorb = 0
  for line in ifile:
    # 1st
    il += 1
    if il == 1:
      ll= line.split()
      natorb = int(ll[1])

    if (iatorb == 0 and il > natorb+3) or (iatorb > 0 and il == natorb+3+iatorb)  :
      ll = line.split()
      u_old = float(ll[0])
      j_old = float(ll[1])

      if u is None: u_new = u_old
      else: u_new = u
      if j is None: j_new = j_old
      else: j_new = j
      line = "  %8.3f %8.3f    U J (Ry)\n"%(u_new,j_new)

    ofile.write(line)

  ifile.close()
  ofile.close()


def f_Reset_Struct(name,latt):
  """
  Reset the struct file in terms of new lattice constant
  """
  case_name = f_Check_Name(name)
  struct_file=case_name.strip()+".struct"
  struct_bak=case_name.strip()+".struct-bak"
  if os.path.isfile(struct_bak) :
    os.remove(struct_bak)

  ## check whether struct file exists
  if not os.path.isfile(struct_file):
    print "ERROR in f_Reset_Latt: struct file " + struct_file + " does not exist"
    sys.exit(1)
  os.rename(struct_file,struct_bak)

  ifile = open(struct_bak,'r')
  ofile = open(struct_file,'w')

  line_count=0
  for line in ifile:
    line_count += 1
    if line_count != 4:
      ofile.write(line)
    else:
      ofile.write("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n" %(latt[0],latt[1],latt[2],latt[3],latt[4],latt[5]))
  ifile.close()
  ofile.close()
  return

def f_Get_Natom(case_name,iop=0):
  """
  Get the information about the number of atoms from the wien2k structure file,
  iop  -- option to control what to return
    0 return the number of non-equivalent atoms
    1 total number of atoms
    2 return the total number of atoms, and the indices of representative non-equivalent atoms
  """
  ifile = open(case_name+".struct",'r')

  f_Skip_Lines(ifile,1)
  line = ifile.readline()
  nat_neq = int( line[27:30] )

  if iop == 0:
    ifile.close()
    return nat_neq

  elif iop == 1:
    f_Skip_Lines(ifile,2)
    nat = 0
    for i_neq in range(nat_neq):
      f_Skip_Lines(ifile,1)
      line = ifile.readline()
      mult = int(line[15:17])
      nat += mult
      f_Skip_Lines(ifile,mult+3)

    ifile.close()
    return nat

  elif iop == 2:
    f_Skip_Lines(ifile,2)
    nat = 0
    index_neq = []
    for i_neq in range(nat_neq):
      index_neq.append(nat+1)
      f_Skip_Lines(ifile,1)
      line = ifile.readline()
      mult = int(line[15:17])
      nat += mult
      f_Skip_Lines(ifile,mult+3)

    ifile.close()
    return nat,index_neq

def f_Read_Struct_w2k(name='',mode=0,debug=False):
  """
  Read the most important structure information from name.struct
  read_mode controls the information to be read
  for mode == 0 :
    latt_type -- lattice type
    latt      -- lattice constants
    basis     -- internal coordinates of all atoms without keeping
                 the information on the equivalence of atoms
  for mode == 1:
    now the basis keeps the equivalence information
  for mode == 2: read all information

  """
  case_name  = f_Check_Name(name)
  struct_file = case_name.strip()+".struct"
  ifile = open(struct_file,'r')

  # title line
  line = ifile.readline()
  title = line.strip()

  # read lattice types and the number of nonequvilanet atoms
  line = ifile.readline()
  latt_type = line[0:4].strip()
  nat_neq = int( line[27:30] )

  # read line on RELA
  line = ifile.readline()
  rela_flag = line[13:17]

  # read lattice constants
  line = ifile.readline()
  a = float( line[0:10] );    b = float( line[10:20]);   c = float( line[20:30])
  alf = float( line[30:40]); bet = float( line[40:50]); gam = float( line[50:60])

  if latt_type == 'R':
    ar=sqrt(3*a*a+c*c)/3.0
    alfa = 2.0*asin(3.0/(2.0*sqrt(3.0+(c/a)**2)))*Rad2Deg
    (a,b,c) = (ar,ar,ar)
    (alf,bet,gam) = (alfa,alfa,alfa)

  latt = [a,b,c,alf,bet,gam]

  # read internal coordinates
  basis = []        # contains the most elementary informat
  misc_atoms = [ ]  # more detailed information
  basis_eq = []     # basis with the information on equivalence atoms
  for i_neq in range(nat_neq):

    xyz=[]
    # coordinates
    line = ifile.readline()
    atom_index = int(line[4:8])
    x = float(line[12:22]); y=float(line[25:35]); z = float(line[38:48])
    xyz.append([x,y,z])

    # mult,isplit
    line = ifile.readline()
    mult = int(line[15:17])
    isplit = int(line[34:37])
    if mult > 1:
      for i_eq in range(mult-1):
        line = ifile.readline()
        x = float(line[12:22]); y=float(line[25:35]); z = float(line[38:48])
        xyz.append([x,y,z])

    # atom_name, npt, r0,rmt,znucl
    line = ifile.readline()

    atom_name = line[0:10]
    npt = int(line[15:20])
    r0 = float(line[25:35])
    rmt = float(line[40:50])
    znucl = int(float(line[55:60]))
    atom = f_Element_Z_to_Symbol(znucl)

    rot = []
    for i in range(3):
      line = ifile.readline()
      x=float(line[20:30]); y=float(line[30:40]); z=float(line[40:50])
      rot.append([x,y,z])

    for i_eq in range(mult):
      basis.append([atom,xyz[i_eq][:]])

    misc_atoms.append([atom_name,atom_index,mult,isplit,npt,r0,rmt,znucl,rot])

    basis_eq.append([atom,xyz[:]])


  ## end the loop over all non-equivalent atoms

  # read symmetry operations
  if mode==2:
    line = ifile.readline()
    nsym = int(line[0:4])
    sym = []
    for i_sym in range(nsym):
      mat=[]; tau=[]
      for i in range(3):
        line = ifile.readline()
        mat.append( [int(line[0:2]), int(line[2:4]), int(line[4:6])] )
        tau.append( float( line[6:16]) )
      sym.append([mat,tau])
      line = ifile.readline()
    # end of loop over isym

    misc = (title,rela_flag,misc_atoms,sym)

  ifile.close()

  if mode == 0:
    return latt_type,latt,basis
  elif mode == 1:
    return latt_type,latt,basis_eq
  else:
    return latt_type,latt,basis_eq, misc

def f_Write_Struct_w2k(file,latt_type,latt,basis,misc=None,debug=False):
  """
  Write the structure in the wien2k struct format in the simple form
  """
  out_name = file.strip()+".struct"

  bak_name = out_name+"_bak"

  if debug: print "Write the structure into " + out_name
  if os.path.isfile(out_name) :
    if debug:
      print " WARNING: the file "+out_name+" exists!"
      print "   -- make a copy to " + bak_name

    if os.path.isfile(bak_name): os.remove(bak_name)
    os.rename(out_name,bak_name)

  ofile = open(out_name,'w')

  if not misc is None:
    (title,rela_flag,misc_atoms,sym) = misc
  else:
    rela_flag = 'RELA'
    sym = []
    misc_atoms=[]
    for iat in range(len(basis)):
      title=file
      atom = basis[iat][0]
      mult = 1
      isplit = 0
      znuc = f_Element_Symbol_to_Z(atom)
      npt = 781
      if znuc < 2:
        r0 = 5.0*1.0E-4
        rmt = 0.50
      elif znuc < 10.0:
        r0 = 1.0*1.0E-4
        rmt = 1.40
      elif znuc < 18.0:
        r0 = 1.0*1.0E-4
        rmt = 1.80
      elif znuc < 36.0:
        r0 = 0.5*1.0E-4
        rmt = 2.00
      elif znuc < 54.0:
        r0 = 0.1*1.0E-4
        rmt = 2.20
      else:
        r0 = 0.05*1.0E-4
        rmt = 2.40

      rot=[[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
      misc_atoms.append([atom,iat,mult,isplit,npt,r0,rmt,znuc,rot])

  nat_neq = len(misc_atoms)

  ofile.write(title+"\n")
  ofile.write("%-4s%-23s%3d\n"%(latt_type,'LATTICE,NONEQUIV.ATOMS:',nat_neq))
  ofile.write("%-13s%4s\n"%('MODE OF CALC=',rela_flag))
  (a,b,c,af,bt,gm) = latt
  ofile.write("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n"%(a,b,c,af,bt,gm))

  iat = 0
  for i_neq in range(nat_neq):
    (x,y,z) = basis[iat][1]; iat += 1
    (atom_name,atom_index,mult,isplit,npt,r0,rmt,znuc,rot) = misc_atoms[i_neq]

    ofile.write("ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n"%(atom_index,x,y,z))
    ofile.write("          MULT=%2d          ISPLIT=%2d\n"%(mult,isplit))
    if mult > 1:
      for i_eq in range(1,mult):
        (x,y,z) = basis[iat][1]; iat += 1
        ofile.write("ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n"%(atom_index,x,y,z))

    ofile.write("%-10s NPT=%5d  R0=%10.8f RMT=%10.5f  Z: %5.2f\n"%(atom_name,npt,r0,rmt,znuc))
    ofile.write("%-20s%10.7f%10.7f%10.7f\n"%('LOCAL ROT MATRIX:',rot[0][0],rot[0][1],rot[0][2]))
    ofile.write("%20s%10.7f%10.7f%10.7f\n"%('',rot[1][0],rot[1][1],rot[1][2]))
    ofile.write("%20s%10.7f%10.7f%10.7f\n"%('',rot[2][0],rot[2][1],rot[2][2]))

  nsym = len(sym)
  ofile.write("%4d NUMBER OF SYMMETRY OPERATIONS\n"%(nsym))
  for isym in range(nsym):
    (mat,tau)=sym[isym]
    for i in range(3):
      ofile.write("%2d%2d%2d%10.7f\n"%(mat[i][0],mat[i][1],mat[i][2],tau[i]))
    ofile.write("%8d\n"%(isym+1))

  ofile.close()

#
# the subroutine to run a wien2k calculation and return the final total energy
#
def f_Calc_Energy(w2k_cmd,name,save_dir,save_flag,init_cmd=None):
  """
  Run a wien2k calculation and return the final total energy
  """
  case_name = f_Check_Name(name)
  bExternal = False
  if not init_cmd is None:
    print "run "+init_cmd + " before run a total energy calculation"
    failure,output = commands.getstatusoutput(init_cmd)
    if failure:
      print "WARNING: error when running " + init_cmd
      f_Write_Log(output)

    clmup=case_name+".clmup"
    if os.path.isfile(clmup) and os.path.getsize(clmup) > 0:
      failure,output = commands.getstatusoutput( init_cmd+" -up")
      if failure:
        print "WARNING: error when running "+init_cmd+" -up"
        f_Write_Log(output)

      failure,output = commands.getstatusoutput(init_cmd+" -dn")
      if failure:
        print "WARNING: error when running "+init_cmd+" -dn"
        f_Write_Log(output)

  failure, output = commands.getstatusoutput( w2k_cmd )
  if failure:
    print "ERROR when running " + w2k_cmd
    print " -- write error output in to " + gparams.pyw2k_errout
    ofile = open(gparams.pyw2k_errout,'w')
    ofile.write(output)
    ofile.close()
    sys.exit(1)
  else:
    f_Write_Log(output)

  etot = f_Read_Etot(case_name)

  if save_dir != '':

    if save_flag == '':
      gparams.save_counter += 1
      save_flag = str(gparams.save_counter)

    w2k_savcmd =  w2k_save_cmd + " -d " + save_dir.strip() + '/  ' + case_name.strip() + "-" + save_flag
    failure,output = commands.getstatusoutput( w2k_savcmd )
    if failure:
      print "WARNING: something wrong  when running " + w2k_savcmd
      print " calculations not saved!!!"

  return etot

def w2k_ChangeStructure_TMDC(aInputCase,fCChange,bKeepInner=True):
    '''
    Generate a new struct of TMDC by specific C axis ( also can be used to extend C axis)

    :param aInputCase: the original case to change from
    :param fCChange: the change of C axis
    :param bChangeInner: whether change inner coordination to keep innerlayer distance ( if not used ,just the same as change c)
    :return: Wien2K_Case

    '''

    fC = 1+fCChange
    aCase=copy.deepcopy(aInputCase)
    aStructure = aCase.aStructure
    aStructure.LatticeParameter[2] = aCase.aStructure.LatticeParameter[2]*fC
    if ( bKeepInner ):
        aStructure.listAtom[1].InnerCoord[0][3] = (aStructure.listAtom[1].InnerCoord[0][3] - 0.75)/fC+0.75
        aStructure.listAtom[1].InnerCoord[3][3] = (aStructure.listAtom[1].InnerCoord[3][3] - 0.75)/fC+0.75
        aStructure.listAtom[1].InnerCoord[1][3] = (aStructure.listAtom[1].InnerCoord[1][3] - 0.25)/fC+0.25
        aStructure.listAtom[1].InnerCoord[2][3] = (aStructure.listAtom[1].InnerCoord[2][3] - 0.25)/fC+0.25
    aCase.a_in0.arFFT[2] = int(aCase.a_in0.arFFT[2]*fC)

    return aCase

def w2k_MultiStructure_TMDC(stCaseNameInput,listChangePercent,bChangeInner=True):
    '''
    Generate .struct & .in0 files for different C-axis TMDC material, keep innerlayer distance of atoms
    '''
    if ( stCaseNameInput == None):
        stCaseName = f_GetCaseName()
    else:
        stCaseName=stCaseNameInput

    aStructure = Wien2K_Structure()
    aStructure.ReadFromFile(stCaseName+".struct")
    a_in0 = Wien2K_File_in0()
    a_in0.ReadFromFile(stCaseName+".in0")

    aStructureNew = copy.deepcopy(aStructure)
    a_in0_New = copy.deepcopy(a_in0)

    for nC in listChangePercent:
        fC = 1+nC/1.0
        aStructureNew.LatticeParameter[2] = aStructure.LatticeParameter[2]*fC

        if ( bChangeInner ):
            aStructureNew.listAtom[1].InnerCoord[0][3] = (aStructure.listAtom[1].InnerCoord[0][3] - 0.75)/fC+0.75
            aStructureNew.listAtom[1].InnerCoord[3][3] = (aStructure.listAtom[1].InnerCoord[3][3] - 0.75)/fC+0.75
            aStructureNew.listAtom[1].InnerCoord[1][3] = (aStructure.listAtom[1].InnerCoord[1][3] - 0.25)/fC+0.25
            aStructureNew.listAtom[1].InnerCoord[2][3] = (aStructure.listAtom[1].InnerCoord[2][3] - 0.25)/fC+0.25

# this subroutine returns the total energy for given percentage changes of the lattice constant)
# Note:
def f_Energy_vs_Latt_Change(x,latt_ini,option,name,w2k_cmd,out_file,save_dir,save_flag,init_cmd=None):
  """
  Calclate the total energy for given lattice constants
  """
  bExternal=False
  case_name = f_Check_Name(name)
  struct_file=case_name.strip()+".struct"
  latt=latt_ini[:]
  type_ini = f_Latt_Type(latt_ini)

  if option == 'vol':    # change volume with fixed a:b:c and lattice angles
    vol = f_Cell_Volume(latt_ini)*(1.0+x[0]/100.0)
    f_Reset_Latt1D('vol',vol,latt)
    latt_info="vol= %12.6f"%(vol)
    out_data = "%12.6f " %(vol)

  elif option == "coa":   # change c/a with fixed (for hexagonal and tetragonal lattice)
    if type_ini != 'hex' and type_ini != 'tetra':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    coa = latt_ini[2]/latt_ini[0]*(1.0+x[0]/100.0); val = coa
    f_Reset_Latt1D('coa',coa,latt)
    latt_info="coa= %12.6f"%(coa)
    out_data = "%12.6f " %(coa)

  elif option == "c":   # relax c only, this is used to optimize layered materials
    c = latt_ini[2]*(1.0+x[0]/100.0)
    latt[2] = c
    latt_info="c= %12.6f"%(c)
    out_data = "%12.6f " %(c)

#  elif option == "tmdc": # used for TMDC structure like Mo/W S/Se2, use external function
#    bExternal= True
#    c = latt_ini[2]*(1.0+x[0]/100.0)
    #use dstart
    #restore lattice file to read
#    aCase = w2k_caseutil.Wien2K_Case() #read current struct
    #aCase.aStructure.LatticeParameter=latt_ini[0:3]
#    aCase.aStructure = w2k_caseutil.Wien2K_Structure(f_Check_Name()+"_init.struct")
#    NewCase = w2k_ChangeStructure_TMDC(aCase,x[0]/100.0)
#    NewCase.WriteToFile(listPart=["struct"])
#    latt_info="tmdc= %12.6f"%(c)
#    out_data = "%12.6f " %(c)

  elif option == "vol_coa":
    if type_ini != 'hex' and type_ini != 'tetra':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    vol = f_Cell_Volume(latt_ini)*(1.0+x[0]/100.0)
    coa = latt_ini[2]/latt_ini[0]*(1.0+x[1]/100.0)
    f_Reset_Latt1D('vol',vol,latt)
    f_Reset_Latt1D('coa',coa,latt)
    latt_info="vol= %12.6f coa= %12.6f"%(vol,coa)
    out_data ="%12.6f %12.6f" %(vol,coa)

  elif option == 'a_c' :
    if type_ini != 'hex' and type_ini != 'tetra':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    a = latt_ini[0]*(1.0+x[0]/100.0)
    c = latt_ini[2]*(1.0+x[1]/100.0)
    latt[0] = a; latt[1] = a; latt[2] = c
    latt_info="a= %10.4f c= %10.4f" %(a,c)
    out_data = "%12.6f %12.6f" %(a,c)

  elif option == "vol_boa_coa":
    if type_ini != 'ortho' and type_ini != 'mono' and type_ini != 'tric':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    vol = f_Cell_Volume(latt_ini)*(1.0+x[0]/100.0)
    boa = latt_ini[1]/latt_ini[0]*(1.0+x[1]/100.0)
    coa = latt_ini[2]/latt_ini[0]*(1.0+x[2]/100.0)
    f_Reset_Latt1D('vol',vol,latt)
    f_Reset_Latt1D('boa_o',boa,latt)
    f_Reset_Latt1D('coa_o',coa,latt)
    latt_info="vol= %12.6f boa= %12.6f coa= %12.6f"%(vol,boa,coa)
    out_data ="%12.6f %12.6f %12.6f" %(vol,boa,coa)

  elif option == "boa_coa":
    if type_ini != 'ortho' and type_ini != 'mono' and type_ini != 'tric':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    boa = latt_ini[1]/latt_ini[0]*(1.0+x[0]/100.0)
    coa = latt_ini[2]/latt_ini[0]*(1.0+x[1]/100.0)
    f_Reset_Latt1D('boa_o',boa,latt)
    f_Reset_Latt1D('coa_o',coa,latt)
    latt_info="boa= %12.6f coa= %12.6f"%(boa,coa)
    out_data ="%12.6f %12.6f" %(boa,coa)

  elif option == 'a_b_c':
    if type_ini != 'ortho' and type_ini != 'mono' and type_ini != 'tric':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    a = latt_ini[0]*(1.0+x[0]/100.0)
    b = latt_ini[1]*(1.0+x[1]/100.0)
    c = latt_ini[2]*(1.0+x[2]/100.0)
    latt[0] = a; latt[1] = b; latt[2] = c
    latt_info="a= %10.4f b= %10.4f c= %10.4f" %(a,b,c)
    out_data = "%12.6f %12.6f %12.6f" %(a,b,c)

  elif option == 'vol_mono':
    if type_ini != 'mono':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    vol = f_Cell_Volume(latt_ini)*(1.0+x[0]/100.0)
    boa = latt_ini[1]/latt_ini[0]*(1.0+x[1]/100.0)
    coa = latt_ini[2]/latt_ini[0]*(1.0+x[2]/100.0)
    # check which lattice angle is not 90
    if   latt_ini[3] != 90.0:
      latt[3] = latt_ini[3]*(1.0+x[3]/100.0)
      bet = latt[3]
    elif latt_ini[4] != 90.0:
      latt[4] = latt_ini[4]*(1.0+x[3]/100.0)
      bet = latt[4]
    elif latt_ini[5] != 90.0:
      latt[5] = latt_ini[5]*(1.0+x[3]/100.0)
      bet = latt[5]

    f_Reset_Latt1D('vol',vol,latt)
    f_Reset_Latt1D('boa_o',boa,latt)
    f_Reset_Latt1D('coa_o',coa,latt)

    latt_info="vol= %12.4f b/a= %10.4f c/a= %10.4f beta= %10.2f" %(vol,boa,coa,bet)
    out_data = "%12.6f %12.6f %12.6f %12.6f" %(vol,boa,coa,bet)

  elif option == 'mono':
    if type_ini != 'mono':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    a = latt_ini[0]*(1.0+x[0]/100.0)
    b = latt_ini[1]*(1.0+x[1]/100.0)
    c = latt_ini[2]*(1.0+x[2]/100.0)
    latt[0] = a; latt[1] = b; latt[2] = c

    # check which lattice angle is not 90
    if   latt_ini[3] != 90.0:
      latt[3] = latt_ini[3]*(1.0+x[3]/100.0)
      bet = latt[3]
    elif latt_ini[4] != 90.0:
      latt[4] = latt_ini[4]*(1.0+x[3]/100.0)
      bet = latt[4]
    elif latt_ini[5] != 90.0:
      latt[5] = latt_ini[5]*(1.0+x[3]/100.0)
      bet = latt[5]

    latt_info="a= %10.4f b= %10.4f c= %10.4f beta= %8.2f" %(a,b,c,bet)
    out_data = "%12.6f %12.6f %12.6f %12.6f" %(a,b,c,bet)

  elif option == 'vol_tric':
    if type_ini != 'tric':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    vol = f_Cell_Volume(latt_ini)*(1.0+x[0]/100.0)
    boa = latt_ini[1]/latt_ini[0]*(1.0+x[1]/100.0)
    coa = latt_ini[2]/latt_ini[0]*(1.0+x[2]/100.0)
    alf = latt_ini[3]*(1.0+x[3]/100.0)
    bet = latt_ini[4]*(1.0+x[4]/100.0)
    gam = latt_ini[5]*(1.0+x[5]/100.0)
    latt[3]= alf; latt[4] = bet; latt[5]= gam
    f_Reset_Latt1D('vol',vol,latt)
    f_Reset_Latt1D('boa_o',boa,latt)
    f_Reset_Latt1D('coa_o',coa,latt)

    latt_info="vol= %12.4f b/a= %10.4f c/a= %10.4f alpha=%12.2f beta= %10.2f gamma= %10.2f" %(vol,boa,coa,alf,bet,gam)
    out_data = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f" %(vol,boa,coa,alf,bet,gam)

  elif option == 'tric':
    if type_ini != 'tric':
      print "ERROR: the option %s is inconsistent with the lattice type %s !!!" %(option,type_ini)
      sys.exit(1)

    a = latt_ini[0]*(1.0+x[0]/100.0)
    b = latt_ini[1]*(1.0+x[1]/100.0)
    c = latt_ini[2]*(1.0+x[2]/100.0)
    alf = latt_ini[3]*(1.0+x[3]/100.0)
    bet = latt_ini[4]*(1.0+x[4]/100.0)
    gam = latt_ini[5]*(1.0+x[5]/100.0)
    latt[0] = a; latt[1] = b; latt[2] = c
    latt[3] = alf; latt[4] = bet; latt[5] = gam
    latt_info="a= %10.4f b= %10.4f c= %10.4f alpha= %8.2f beta=%8.2f gamma= %8.2f" %(a,b,c,alf,bet,gam)
    out_data = "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f" %(a,b,c,alf,bet,gam)
  else:
    print "ERROR: unsupported option for option = " + option
    sys.exit(1)

  f_Reset_Struct(case_name,latt)
  if ( bExternal ):
     NewCase.WriteToFile(listPart=['struct'])
  etot = f_Calc_Energy( w2k_cmd, case_name,save_dir,save_flag, init_cmd)

  info = "%6d " %(gparams.save_counter) + out_data + " %16.6f" %(etot)
  print "#w2k: "+info
  if  out_file != '' :
    ofile = open(out_file,'a')
    ofile.write(info+"\n")
    ofile.close()

  f_Write_Log(latt_info + " Etot= %16.6f" %(etot))

  return etot

def f_Prepare_Para(nproc,nmpi):
  """
  This function prepare parallelization
  """
  print "Not implemented yet!"
