#!/usr/bin/env python
import sys,os,shutil
from math import *
from chem_utils import *
from struct_utils import * 

mod_name = "siesta_utils"

def f_fdf_reset(fdf,name,val):
  """
  Reset the value of a paramter in a fdf file  
  """ 
  sname = "f_fdf_reset"

  fdf_bak=fdf+"-bak"
  os.rename(fdf,fdf_bak)
  ifile = open(fdf_bak,'r')
  ofile = open(fdf,'w') 
  while(1):
    line = ifile.readline()
    if not line: break 
  
    line_sp =line.split()
    if len(line_sp) > 1 and line_sp[0] == name: 
      line_sp[1] = val
      line_out=''
      for i in range(len(line_sp)):
        line_out=line_out+line_sp[i]+' '
      line_out=line_out+"\n"
    else:
      line_out = line 
    ofile.write(line_out) 

  ofile.close()
  ifile.close()

def f_Read_XV(sys_label,debug=True):
  """
  Read XV file and return the structural information
  """
  try: 
    ifile = open(sys_label.strip()+".XV",'r')
  except:
    print "Fail to open ", sys_label.strip(),".XV"
  
  # get lattice vectors
  latt_vec = [] 
  for i in range(3):
    line = ifile.readline()
    data = line.split()
    x = float(data[0]) 
    y = float(data[1])
    z = float(data[2])
    latt_vec.append([x,y,z])

  if debug: 
    print "Lattice Vectors:"
    for i in range(3):  print latt_vec[i]
   
  # get number of atoms 
  line = ifile.readline()
  nat = int(line.split()[0])
  if debug: 
    print "Number of atoms=",nat
    print "Coordinates:"

  itype=[]
  iznuc=[]
  basis=[]
  for i in range(nat):
    line  = ifile.readline()
    data  = line.split()
    znucl = int(data[1])
    x = float(data[2])
    y = float(data[3])
    z = float(data[4])
    atom = f_Element_Z_to_Symbol(znucl)
    basis.append( [atom,[x,y,z]] )
    if debug: 
      print " %4s %12.6f %12.6f %12.6f" %(atom,x,y,z)
  return basis,latt_vec

def f_Read_Eig(sys_label):
  """
  Read sys.EIG file 
  """
  ifile = open(sys_label+".EIG",'r')

  data = ifile.readline().split() 
  efermi = float(data[0]) 
 
  data = ifile.readline().split()
  norb = int(data[0])
  nsp  = int(data[1]) 
  nk   = int(data[2]) 
  
  eigk = []
  kpoints = []
  for ik in range(nk):
    data = ifile.readline().split()
    kp = [ float(data[0]), float(data[1]), float(data[2]) ] 
    kw = float(data[3])
    kpoints.append([kp,kw]) 

    eig = []
    for isp in range(nsp): 
      for iorb in range(norb):
        data = ifile.readline().split()
        eig.append(float(data[0]))
    eigk.append(eig) 
  return eigk,efermi,nsp,kpoints

def f_Plot_Bands(sys_label,erange=(-20.0,20.0),out='gnp',ezero='vbm'):
  """
  Read sys_label.bands corresponding to BandLines and export band energies 
  into the format that is suitable for plotting 
  """
  sname = 'f_Plot_Bands' 

  (ebmin,ebmax) = erange 

  ifile = open(sys_label + ".bands",'r') 
  data = ifile.readline().split()
  efermi = float(data[0])
  
  data = ifile.readline().split()
  (kmin,kmax) = (float(data[0]),float(data[1]))

  data = ifile.readline().split()
  (emin,emax) = (float(data[0]),float(data[1]))
  
  data = ifile.readline().split()
  (nb,nsp,nk) = (int(data[0]),int(data[1]),int(data[2]))

  data = ifile.read().split()
  
  eb = []
  for isp in range(nsp): 
    eb_s=[] 
    for ib in range(nb): 
      ebk = [] 
      for ik in range(nk):
        ebk.append(0.0) 
      eb_s.append(ebk) 
    eb.append(eb_s) 

  i = 0
  kp= []
 
  (nvm,kvm,evm) = (0,    0,-10000.0) 
  (ncm,kcm,ecm) = (1000, 0, 10000.0) 
  for ik in range(nk): 
    kp.append(float(data[i]))
    i += 1
    for isp in range(nsp):
      for ib in range(nb):
        eb[isp][ib][ik] = float(data[i]) 
        et = eb[isp][ib][ik]
        if et > efermi and et < ecm: 
          (ncm,kcm,ecm) = (ib,ik,et) 
        if et < efermi and et > evm:
          (nvm,kvm,evm) = (ib,ik,et) 
        i += 1 
      # ib
    # isp
  # ik

  print " Fermi Energy = %12.4f"%(efermi)  
  print " Valence    band maximum: n=%5d, ik==%5d, E_VBM= %12.4f "%(nvm,kvm,evm) 
  print " Conduction band maximum: n=%5d, ik==%5d, E_CBM= %12.4f "%(ncm,kcm,ecm) 
  print " Band Gap = %12.4f"%(ecm-evm)  
  
  ofile = open(sys_label+'-band-'+out+".dat",'w')
  print "Write band energies into the file",sys_label+'-band-'+out+".dat" 
  ofile.write("# Band energies data written by "+mod_name+'::'+sname+'\n')
  if ezero == 'fermi':
    eshift = efermi 
    ofile.write("# All bands shifted so that Fermi energy (middle gap) is zero\n") 
  else:
    eshift = evm 
    ofile.write("# All bands have been shifted so that VBM is zero\n") 

  for isp in range(nsp):
    for ib in range(nb): 
      for ik in range(nk):
        et = eb[isp][ib][ik]-eshift
        if et > ebmin and et < ebmax: 
          ofile.write("%f12.6 %f12.6\n"%(kp[ik],et))
      ofile.write("\n")

  ofile.close()   


def f_Read_Rho(sys_label,ftag,mode='',debug=True):
  """
  Read the potential file as written out 

  """
  ifile_name = sys_label.strip()+"."+ftag.strip()
  try:
    ifile = open(ifile_name,'r'+mode)
  except:
    print "Fail to open ", ifile_name
    sys.exit(1)

  # get lattice vectors
  latt_vec = []
  if debug: print "Lattice Vectors:"
  for i in range(3):
    data = ifile.readline().split()
    x = float(data[0])
    y = float(data[1])
    z = float(data[2])
    latt_vec.append([x,y,z])

    if debug: print "a(%i)= (%12.6f,%12.6f,%12.6f)" %(i,x,y,z)

  # get nx,ny,nz,nspin
  data = ifile.readline().split()
  nx = int(data[0])
  ny = int(data[1])
  nz = int(data[2])
  nsp= int(data[3]) 

  rho=[]
  for isp in range(nsp):
    for iz in range(nz):
      for iy in range(ny): 
        for ix in range(nx):
          line = ifile.readline()
          rho.append(float(line)) 

  ifile.close()
  return rho,(nx,ny,nz),latt_vec,nsp    

def f_Create_fdf(name):
  sname = "f_Create_fdf"

  ofile = open(name.strip()+".fdf",'w')
  print "Create new SIESTA fdf: "+name.strip()+".fdf"

  # check the number of species
  ofile.write("# SIESTA fdf file generated by %s @ %s\n" %(sname,mod_name))
  ofile.write("\n# system info\n")
  ofile.write("SystemName %s\n"  %(name))
  ofile.write("SystemLabel %s\n" %(name))

  ofile.write("""
# General parameters
UseSaveData yes

# SCF parameters
MaxSCFIterations 200
ElectronicTemperature 300.0 K      # increases this for difficult cases
OccupationFunction FD              # changing to MP may be helpful for difficult cases
OccupationMPOrder 4                # effective only for OccupationFunction=MP
DM.NumberPulay 6
DM.MixingWeight 0.25               # reduce this number for difficcult cases, but not too small

# relaxation paramters
MD.NumCGsteps 0                    # the number of steps for structural optimization
MD.TypeOfRun  cg

# Basis set parameters
PAO.BasisType    split
PAO.BasisSize    DZP

# Output options
#SaveElectrostaticPotential yes    # needed for workfunction
#SaveDeltaRho               yes    # total density minus overlapped atomic density
WriteCoorXmol               yes    # write molecular structure in the SystemLabel.xyz
  """)
  ofile.write("\n")
  ofile.close()


# write lattice structure into the siesta format
def f_Write_Struct_siesta(name,mol,latt_const=None,latt_vec=None,full_fdf=True):
  """
  Write the crystal structure information to the siesta fdf format 
  """
  sname = "f_Write_Struct_siesta"
  
  fdf = name.strip()+".fdf"
  if full_fdf is True:
    f_Create_fdf(name)
    ofile = open(name.strip()+".fdf",'a')
  else:
    ofile = open(fdf,'w')

  # check the number of species 
  nat = len(mol)
  nsp,species,sp_index = f_Check_Species(mol)

  ofile.write("# SIESTA fdf file generated by %s @ %s\n" %(sname,mod_name))
 
  ofile.write("\n# structural information \n")
  ofile.write("\nNumberOfAtoms\t %6d\n" %(nat))
  ofile.write("NumberOfSpecies\t %6d\n" %(nsp))

  ofile.write("%block ChemicalSpeciesLabel\n")
  for isp in range(nsp):
    znucl= f_Element_Symbol_to_Z(species[isp])
    ofile.write("%6d\t%6d\t%6s\n"%(isp+1,znucl,species[isp]))
  ofile.write("%endblock ChemicalSpeciesLabel\n\n") 

  if latt_const is None: 
    xyz_scale = 1.0 
    ofile.write("AtomicCoordinatesFormat\t Ang\n")
  else: 
    xyz_scale=latt_const
    a0=latt_const
    ofile.write("LatticeConstant\t %12.6f Ang\n"%(a0))
    ofile.write("kgrid_cutoff\t %12.6f Ang\n\n"%(a0*2))
    ofile.write("%block LatticeVectors\n")
    for i in range(3): 
      ofile.write("%10.6f %10.6f %10.6f\n"%(latt_vec[i][0]/a0,latt_vec[i][1]/a0,latt_vec[i][2]/a0))
    ofile.write("%endblock LatticeVectors\n\n")
    ofile.write("AtomicCoordinatesFormat\t ScaledCartesian\n")

  ofile.write("%block AtomicCoordinatesAndAtomicSpecies\n")
  for iat in range(nat):
    xyz  = mol[iat][1][:]
    for i in range(3): xyz[i] /= xyz_scale 
    atom = mol[iat][0]
    isp = sp_index[iat]
    ofile.write("%12.6f %12.6f %12.6f %6d # %6s %6d\n" %(xyz[0],xyz[1],xyz[2],isp,atom,iat+1))
  ofile.write("%endblock AtomicCoordinatesAndAtomicSpecies\n\n")

  ofile.close()

def f_SIESTA_Read_Struct_Out(name,type='m'):
  """
  Read structure information from name.STRUCT_OUT 
  """
  ifile = open(name.strip()+".STRUCT_OUT",'r')
 
  # read lattice vectors
  latt_vec=[]
  for i in range(3):
    line = ifile.readline().split()
    vec=[]
    for j in range(3): 
      vec.append(float(line[j]))
    latt_vec.append(vec) 

  # get the number of atoms 
  line = ifile.readline()   
  natom = int(line)

  mol=[]
  for ia in range(natom): 
    line = ifile.readline().split()
    znuc = int(line[1]) 
    atom = f_Element_Z_to_Symbol(znuc)
    xyz = [float(line[2]),float(line[3]),float(line[4])]

    if type == 'm': ## a molecule structure is read 
      xyz_old = xyz[:]
      xyz = [0.0, 0.0, 0.0]
      for j in range(3):
        for i in range(3):  
          xyz[j] += xyz_old[i]*latt_vec[i][j]
    mol.append([atom,xyz])

  if type == 'm':
    return mol
  else:
    return mol,latt_vec

