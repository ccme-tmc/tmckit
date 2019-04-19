#!/usr/bin/env python

import copy
import commands,sys,string,math,re,subprocess
import glob,os.path,shutil

from io_utils import *
from bandstruct import band_gap_analysis 

def vasp_set_incar(v_name,v_val,v_type=None):
  """
  set or reset the input parameter v_name in INCAR to the value of "v_val"
  """
  # make a copy of original INCAR file 

  if v_type is None:
    if isinstance(v_val,str):
      v_type = 's'
    elif isinstance(v_val,int):
      v_type = 'i'
    elif isinstance(v_val,float):
      v_type = 'f'
    else:
      print "ERROR in vasp_set_incar: fail to recongnize the type of the input value!"
      sys.exit(0)

  incar='INCAR'
  incar_bak = incar+"_old"
  if os.path.isfile(incar_bak): os.remove(incar_bak)
  os.rename(incar,incar_bak)

  ifile = open(incar_bak,'r')  
  ofile = open(incar,'w')  
  lines = ifile.readlines()
 
  done = False
  for line in lines:
    tmp = line.strip()
    if len(tmp) == 0 or tmp[0] == '#':  
      ofile.write(line) 
      continue 

    line_out = ''
    if v_name.lower() in line.lower(): 

      # first check whether the line contains comments 
      tmp = line.split('#') 
      if len(tmp) > 1: 
        line1 = tmp[0]
        comments = " # " + tmp[1] +'\n'
      else:
        line1 = line 
        comments = '\n'

      # check whether the line contains multiple input parameters 
      line_s = line1.split(';') 
      for i in range(len(line_s)):
        if v_name.lower() in line_s[i].lower():
          name_val = line_s[i].split('=') 
          if v_type == 'f':
            new_entry = name_val[0] + " = %f "%(v_val) 
          elif v_type == 'i':
            new_entry = name_val[0] + " = %d "%(v_val)  
          elif v_type == 's':
            new_entry = name_val[0] + " = "+ v_val
          else: 
            print "ERROR: the data type %s is not supported yet!"%(v_type)
            print "--- %s is not reset!"%(v_name) 
            new_entry = line_s[i] 
          done = True 

          line_s[i] = new_entry
          
        if i == 0: 
          line_out = line_s[i]
        else:
          line_out = line_out + " ; " + line_s[i]
        
        line_out = line_out + comments 
    else:
      line_out = line   

    ofile.write(line_out) 
      
  if not done: 
    if v_type == 'f':
      ofile.write(" %s = %f \n"%(v_name,v_val))
    elif v_type == 'i':
      ofile.write(" %s = %d \n"%(v_name,v_val))
    elif v_type == 's':
      ofile.write(" %s = %s \n"%(v_name,v_val))
    else:
      print "ERROR: the data type %s is not supported yet!"%(v_type)

  ofile.close() 

def vasp_set_vol(x):
  # make a copy of the old POSCAR file 
  poscar='POSCAR'
  poscar_bak = 'POSCAR_old'
  if os.path.isfile(poscar_bak): os.remove(poscar_bak)
  os.rename(poscar,poscar_bak)

  ifile = open(poscar_bak,'r')  
  ofile = open(poscar,'w')  

  lines = ifile.readlines()

  tmp = lines[1].split()
  volx_old = float(tmp[0])
  volx = volx_old * x

  nl = len(lines)
  for il in range(nl):
    if il == 1:
      ofile.write("%12.6f\n"%(volx))
    else:
      ofile.write(lines[il])
  ifile.close()
  ofile.close()

def vasp_read_eigval(eigfile="EIGENVAL"):
  """
  Read eigen energies from EIGENVAL
  """
  ifile = open(eigfile,'r')

  # skip the first five lines  
  for i in range(5):
    ifile.readline()
  
  line_s = ifile.readline().split()
  ne = int(line_s[0])  
  nk = int(line_s[1])
  nb = int(line_s[2])
  kvecs = []
  ebands = []
  kwts = []
  for ik in range(nk):
    ifile.readline()
    line_s = ifile.readline().split()
    kv = []
    for i in range(3):
      kv.append(float(line_s[i]))
    kvecs.append(kv)
    kwts.append(float(line_s[3]))
 
    ek = []
    for ib in range(nb):
      line_s = ifile.readline().split()
      ek.append(float(line_s[1]))
    ebands.append(ek)

  ifile.close()

  return ebands,kvecs,kwts,ne
  

def vasp_incar_band(ecut,nb,ax=0.0,scr=0.0,loptic=False):
  """
  Create a INCAR file with the defaul setting for a band structure calculation 
  """

  ofile = open('INCAR','w')
  ofile.write("general:\n")
  if ecut == 0.0:
    ofile.write(" PREC = Normal \n")
  elif ecut == -1.0:
    ofile.write(" PREC = Accurate \n")
  else:
    ofile.write(" ENCUT = %8.2f\n"%(ecut))

  ofile.write(" ISMEAR = -5\n")
  if nb > 0:
    ofile.write(" NBANDS = %d\n"%(nb))

  if loptic:
    ofile.write(" LOPTICS = .TRUE.; NEDOS = 2000 \n")

  if ax > 0.0:
    ofile.write(" LHFCALC = .TRUE.; ALGO=All; NKRED= 2; PRECFOCK = Fast; AEXX = %8.3f\n"%(ax))
    if scr > 0.0:
      ofile.write(" HFSCREEN = %8.3f\n"%(scr))
    elif scr < 0.0:
      ofile.write(" HFSCREEN = %8.3f ; LTHOMAS = .TRUE. \n"%(abs(scr)))

  ofile.close()

def vasp_sah_incar(ecut,nb,ax=0.0,scr=0.0,loptic=True,hfprec="fast"):
  """
  Create a INCAR file with the defaul setting for a band structure calculation 
  """

  ofile = open('INCAR','w')
  ofile.write("general:\n")
  if ecut == 0.0:
    ofile.write(" PREC = Normal \n")
  elif ecut == -1.0:
    ofile.write(" PREC = Accurate \n")
  else:
    ofile.write(" ENCUT = %8.2f\n"%(ecut))

  ofile.write(" ISMEAR = -5\n")
  ofile.write(" NBANDS = %d\n"%(nb))
  if loptic:
    ofile.write(" LOPTICS = .TRUE.; NEDOS = 2000 \n")

  if ax > 0.0:
    ofile.write(" LHFCALC = .TRUE.; ALGO=All; AEXX = %8.3f\n"%(ax))
    if scr > 0.0:
      ofile.write(" HFSCREEN = %8.3f\n"%(scr))
    elif scr < 0.0:
      ofile.write(" HFSCREEN = %8.3f ; LTHOMAS = .TRUE. \n"%(abs(scr)))

    if hfprec == 'fast':
      ofile.write(" NKRED= 2; PRECFOCK = Fast\n")
    elif hfprec == 'accurate':
      ofile.write(" PRECFOCK = Normal\n")
    else:
      ofile.write(" PRECFOCK = Accurate\n")

  ofile.close()

def vasp_sah_aexx(ax,emi):
  """
  Obtain the value of aexx that satisfies the requirement 
        ax = 1/epsm(ax) == emi(ax)
  by assuming 
      emi(ax) = a + b*ax + c*ax**2
  """

  if len(ax)==1:
    return emi[0]

  elif len(ax) == 2:
    return ax[1]*emi[0]/(ax[1]+emi[0]-emi[1])

  elif len(ax)==3:
    (e0,e1,e2) = (emi[0],emi[1],emi[2])
    (ax1,ax2)  = (ax[1],ax[2])
    a = e0
    c = (ax2*(e1-e0) - ax1*(e2-e0))/(ax1*ax2*(ax1-ax2))
    b = (e2 - e0 - c*ax2*ax2)/ax2
  else:
    import numpy 
    coef = numpy.polyfit(ax,emi,2)
    a = coef[2]
    b = coef[1]
    c = coef[0]

  print "a,b,c=%8.3f %8.3f %8.3f"%(a,b,c)

  as1 = ( - (b-1.0) + math.sqrt( (b-1.0)**2 - 4*a*c) )/(2*c)
  as2 = ( - (b-1.0) - math.sqrt( (b-1.0)**2 - 4*a*c) )/(2*c)
  print "as1,as2= %8.3f %8.3f"%(as1,as2)
  if (as1 < 0.0 or as1 > 1.0) and (as2 < 0.0 or as2 > 1.0):
    print "ERROR: no appropriate value of alpha has been found by quadratic fitting!"
    sys.exit(1)
  elif as1 > 0.0 and as1 < 1.0:
    return as1
  else:
    return as2

def vasp_read_pot(potfile="LOCPOT"):
  """
  Read the electron density or potential from the file 
  """
  ifile = open(potfile,'r')
  ifile.readline()

  latt_scale = float(ifile.readline())  # scale for the lattice vectors

  # read the lattice vectors 
  latt_vec=[]
  for i in range(3):
    line_s = ifile.readline().split()
    x = float(line_s[0])
    y = float(line_s[1])
    z = float(line_s[2])
    latt_vec.append([x,y,z]) 
  
  # read species
  specs = ifile.readline().split()

  nat_sp =[]
  line_s = ifile.readline().split()
  nat = 0
  for i in range(len(specs)):
    nat_sp.append(int(line_s[i]))
    nat += nat_sp[i]
  
  ifile.readline()
  for i in range(nat):
    line = ifile.readline()
  
  ifile.readline()

  line_s = ifile.readline().split()
  nx = int(line_s[0])
  ny = int(line_s[1])
  nz = int(line_s[2]) 
  ndata = nx*ny*nz

  nlines = ndata/5
  if ndata%5 != 0: nlines += 1
  data = []
  for il in range(nlines):
    line_s = ifile.readline().split()
    for i in range(len(line_s)):
      data.append(float(line_s[i]))
  ifile.close()
  return data,(nx,ny,nz),latt_vec
  

def vasp_read_gwout(outfile="OUTCAR"):
  """
  Read basic information from GW ouput 
  """

  nkp = vasp_getout('nk',outfile)
  tag_begin ='QP shifts <psi_nk| G(iteration)W_0 |psi_nk>'

  tag_stop ='-------------------------------------------'
  
  data = io_read_lines_tagged(outfile,tag_begin,-1,tag_stop,-1) 

  # remove the first three lines
  for i in range(3): data.pop(0)
  data.pop()                         # remove the last empty line 

  nbgw = len(data)/nkp - 4
  print "Number of bands",nbgw 
  kvecs = []
  ebands_ks = []
  ebands_gw = []
  
  il0 = 0 
  evbm_ks = -1.0E10
  ecbm_ks =  1.0E10
  evbm_gw = -1.0E10
  ecbm_gw =  1.0E10
  for ik in range(nkp):
    line_s = data[il0].split()
    kx=float(line_s[3])
    ky=float(line_s[4])
    kz=float(line_s[5])
    kvecs.append([kx,ky,kz])

    eb_ks=[]
    eb_gw=[]
    for ib in range(nbgw):
      line_s = data[il0+3+ib].split()
      eks = float(line_s[1])
      egw = float(line_s[2]) 

      eb_ks.append(eks)
      eb_gw.append(egw) 
      
      occ = float(line_s[-1]) 
      if occ > 0.0: 
        if eks > evbm_ks: 
          evbm_ks = eks
          nvbm_ks = ib 

        if egw > evbm_gw: 
          evbm_gw = egw 
          nvbm_gw = ib 

      else:

        if eks < ecbm_ks: 
          ecbm_ks = eks 
          ncbm_ks = ib 

        if egw < ecbm_gw: 
          ecbm_gw = egw 
          ncbm_gw = ib 

    ebands_ks.append(eb_ks)
    ebands_gw.append(eb_gw) 

    il0 = il0 + nbgw + 4

  return kvecs,ebands_gw,evbm_gw,ecbm_gw, ebands_ks,evbm_ks, ecbm_ks 
  
def vasp_read_dos(dosfile="DOSCAR",iat_pdos=None,debug=True):
  """
  Read DOS and PDOS from DOSCAR and write into the file in the format 
  that can be used directly by plotting promgrams like gnuplot or xmgrace etc
  """
  ifile = open(dosfile,'r') 
  io_skip_lines(ifile,5) 

  # get some info on DOS
  line_s= ifile.readline().split()
  print line_s 
  emax  = float(line_s[0])
  emin  = float(line_s[1])
  nedos = int(line_s[2])
  efer  = float(line_s[3]) 

  if debug: 
    print "emax=%12.6f emin=%12.6f nedos=%6d efer=%12.6f"%(emax,emin,nedos,efer)
  
  # read the total density of states
  dos=[]
  for i in range(nedos):
    line_s =  ifile.readline().split()

    # check whether it is spin-polarized 
    if i== 0: 
      if len(line_s) == 3: 
        nsp = 1
      else:
        nsp = 2

    e = float(line_s[0])
    if nsp == 1:
      d = float(line_s[1])
      dos.append([e,d])
    else:
      dup = float(line_s[1])
      ddn = float(line_s[2])
      dos.append([e,dup,ddn])

  # read PDOS
  if iat_pdos is None:
    ifile.close()
    return dos,efer,emax,emin 

  pdos =[]
  iat = 0 
  while (1):
    line = ifile.readline()
    if not line:
      ifile.close()
      return dos,pdos,efer,emax,emin 

    iat += 1
    if iat_pdos == iat or iat_pdos == 0 : 
      pdos_i = []
      for i in range(nedos): 
        line_s = ifile.readline().split()
        tmp = []
        for j in range( 1+nsp*9 ):
           tmp.append(float(line_s[j]))
        pdos_i.append(tmp)
      pdos.append(pdos_i)

def vasp_readband_out(outfile='OUTCAR',debug=False,mode=1,out=0,occ_zero=0.01):
  """
  Read band energies from OUTCAR
    mode: control which data that meets the requirements are returned 
      0 return all data 
      1 return the first set
     -1 return the last ones 
  """
  # first get some dimension parameters
  nk = vasp_getout('nk',outfile)
  nb = vasp_getout('nb',outfile)
  efer = vasp_getout('efer',outfile)

  print "Read band energies from "+ outfile
  print "  Number of k points:", nk
  print "  Number of bands:   ", nb
  print "  Fermi Energy:      ", efer

  lines = io_read_lines_tagged(outfile,'E-fermi :',-1,'-------------',mode=1)
  del lines[0:2]

  kvecs=[]
  ebands=[]
  evbm = -1000.0
  ecbm = 1000.0 
  for ik in range(nk):
    i = ik * (nb+3)
    kvecs.append([float(x) for x in lines[i].split()[3:6]])

    eb_k = []
    i += 1
    for ib in range(nb):
      i += 1
      tmp = lines[i].split()
      en = float(tmp[1])
      occ = float(tmp[2])
      eb_k.append(en) 
      if occ > occ_zero and en > evbm: evbm = en 
      if occ < occ_zero and en < ecbm: ecbm = en 
     
    ebands.append(eb_k) 
  
  print "  Valence band maximum:  ",evbm
  print "  Conduc. band minimum:  ",ecbm
  egap = ecbm - evbm 
  if evbm > ecbm: 
    print "  Metallic system!"
  else:
    if efer < evbm or efer > ecbm: 
      print "  WARNING: insulating but the Fermi energy is not within [evbm, ecbm]:"
      efer = 0.5*(evbm+ecbm)

  print ":BandGap= ",egap 

  if out == 0:
    return egap,evbm,efer
  else:
    return kvecs,ebands,efer

def vasp_get_spec(posfile="POSCAR"):
  # get the information on the composition of the system
  specs  = io_get_line(posfile,6).split()
  l_nat  = io_get_line(posfile,7).split()

  nspec = len(specs)  # number of species
  nat_spec = []       # number of atoms per species
  nat_tot = 0           # total number of atoms
  for i in range(nspec):
    nat = int(l_nat[i])
    nat_spec.append(nat)
    nat_tot += nat

  return specs,nat_spec

def vasp_get_iat(specs,nat_spec,at,iat_sp):
  """
  return the index of iat_sp-th atom "at" in the whole structure 
  """
  iat = 0
  for isp in range(len(specs)):
    if at == specs[isp]:
      iat += iat_sp
      break 
    else:
      iat += nat_spec[isp]
  return iat,isp+1 

def vasp_get_nel(poscar="POSCAR",potcar="POTCAR"):
  """
  get the total number of valence electrons 
  """
  # get the line with the number of atoms per species from POSCAR
  line = io_get_line(poscar,7) 
  tmp = line.split()

  # check whether the 6-th or 7-th line contains nat_spec
  try: 
    itmp = int(tmp[0])  
  except: 
    line = io_get_line(poscar,6)
    tmp = line.split()
  
  nat_spec = []
  for i in range(len(tmp)):
    nat_spec.append(int(tmp[i]))

  # get the number of valence electrons per species 
  zvals = io_grep_lines(potcar,'ZVAL',0,6)
  nel = 0.0 
  for i in range(len(nat_spec)):
    zval = float(zvals[i])
    nel += nat_spec[i]*zval
    print "ZVAL(%d) = %4.1f"%(i,zval)

  return nel

def vasp_getout(var0,outfile='OUTCAR',debug=False):
  """
  Extract information from vasp OUTCAR file
  """
  var = var0.lower()
  if var == 'nb':
    val = int(io_grep_lines(outfile,'NBANDS=',1,-1,debug=debug))
    if debug: print var+"= ",val
  elif var == 'nat':
    val = int(io_grep_lines(outfile,'NIONS =',1,-1,debug=debug)) 
  elif var == 'nk':
    val = int(  io_grep_lines(outfile,'NKPTS =',1,4,debug=debug))
  elif var == 'nsp':
    val = int(io_grep_lines(outfile,'ISPIN  =',1,3,debug=debug))
  elif var == 'cput':
    val = float(io_grep_lines(outfile,'Total CPU time used',-1,-1,debug=debug))
  elif var == 'mag':
    val = float(io_grep_lines(outfile,'number of electron',-1,-1,debug=debug))

  elif var == 'efer':
    val = float(io_grep_lines(outfile,'E-fermi',-1,3,debug=debug))
  elif var == 'epsm':
    lines = io_read_lines_tagged(outfile,"REAL DIELECTRIC FUNCTION",3)[0]
    tmp = lines[-1].split()
    val = []
    for i in range(3):
      val.append(float(tmp[i+1]))
  elif var == 'epsm_pol': # read the epsm from the polarization approach
    lines = io_read_lines_tagged(outfile,"MACROSCOPIC STATIC DIELECTRIC TENSOR",4,mode=1)
    val = []
    for i in range(3):
      tmp = lines[i+1].split()
      val.append(float(tmp[i]))

  elif var == 'etot':
    val = float(io_grep_lines(outfile,'energy(sigma->0)',-1,-1,debug=debug))
  elif var == 'vol':
    val = float(io_grep_lines(outfile,'volume of cell',-1,-1,debug=debug))
  
  elif var == 'lat':
    lines = io_read_lines_tagged(outfile,"direct lattice vectors",iop=3,mode=-1)
    val = []
    for i in range(3):
      xyz=[]
      for j in range(3): 
        xyz.append(float(lines[i].split()[j]))
      val.append(xyz)
  elif var == 'gap':
    egap,evbm,efer = vasp_readband_out(outfile,debug=debug,out=0)
    val = egap 

  elif var == 'evbm':
    (egap,evbm,efer) = vasp_readband_out(outfile,debug=debug,out=0)
    val = evbm

  elif var == 'efer':
    (egap,evbm,efer) = vasp_readband_out(outfile,debug=debug,out=0)
    val = efer

  elif var == 'band':
    (kvecs,ebands,efer) = vasp_readband_out(outfile,debug=debug,out=1)
    val = band_gap_analysis([ebands],efer,kvecs)
  
  return val

def vasp_write_kpoints(nks,mode='G',sh=None):
  
  if mode == 'G':
    ofile = open('KPOINTS','w')
    ofile.write("K-Points\n")
    ofile.write("0 \n")
    ofile.write("Gamma \n")
    ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
    if sh is None:
      ofile.write("0 0 0\n")
    else:
      ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))

  elif mode == 'M':
    ofile = open('KPOINTS','w')
    ofile.write("K-Points\n")
    ofile.write("0 \n")
    ofile.write("Monkhorst-Pack \n")
    ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
    if sh is None:
      ofile.write("0 0 0\n")
    else:
      ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))

    ofile.close()

def vasp_run(np=1,out=None):
  """
  Run vasp
  """

  ierr = 0
  if np == 1:
    vasp_cmd = "vasp"
  else:
    vasp_cmd = "mpirun -np %d vasp"%(np)

  if not out is None:
    vasp_cmd = vasp_cmd +">& "+out

  failure,output = commands.getstatusoutput( vasp_cmd )

  if failure:
    print "ERROR when running " + vasp_cmd
    ierr = 1

  print output
  return ierr

def vasp_get_data(dir_root='.',var='etot',tag=''):
    """
    read the data from all subdirectories in dir_run and return the data as a directory 
    with the name of each subdirectory as the key
    """
    # get the current working directory
    wdir = os.getcwd()

    cases_orig = glob.glob(dir_root+ "/*"+tag+"*")

    cases_all = []
    for i in range(len(cases_orig)):
        case = cases_orig[i]
        case_name = os.path.basename(case)

        if os.path.isdir(case):  
            cases_all.append(case_name)

    n_cases = len(cases_all) 

    dict_out = {}
    for case in cases_all:
        dir_case = dir_root + "/" + case
        os.chdir(dir_case)

        try:
            val =  vasp_getout(var)
            dict_out[case] = val 
        except:
            print "WARNING: fail to read %s from %s"%(var,case) 
        os.chdir(wdir) 

    return dict_out 



def vasp_run1(np=1,out=None,err=None):
  """
  Run vasp 
  """

  ierr = 0
  if np == 1:
    vasp_cmd = "vasp"
  else:
    vasp_cmd = "mpirun -np %d vasp"%(np)

  if out is None:
    ofile = subprocess.PIPE
  else:
    ofile = open(out,'w')

  if err is None:
    efile = subprocess.PIPE
  else:
    efile = open(err,'w')

  p=subprocess.Popen(vasp_cmd,stdout=ofile,stderr=efile)
  p.wait()

  if not out is None: ofile.close()
  if not err is None: efile.close()

