#!/usr/bin/env python
import sys,os,shutil,re
from chem_utils import *
from constants import *
from struct_utils import * 
from latt_utils import *
from math import *

qe_default_psp = ".pbe-van"
kgrid_cutoff = 40.0


def f_qe_read_struct_out(file,iop_out=0):
  """ 
  Read struct from QE input or output file
    iop_out = 0 -- return lattice vectors and basis (internal coordinates) 
            = 1 -- return Cartesian coordinates in units of angstrom  
  """

  ifile= open(file,'r')

  found_none = False
    
  ibrav = -1
  natom = 0 
  lat_vec = None
  vol = None 
  while (1):
    line = ifile.readline()
    if not line: 
      found_none = True
      break

    # ibrav
    if ibrav == -1 and "bravais-lattice index" in line:
      info = line.split()
      ibrav = int(info[3])    
      print "bravais-lattice index (ibrav) = %2d "%(ibrav)
  
    # alat 
    if "lattice parameter" in line:
      data = line.split()
      alat = float(data[4])
      print "lattice parameter (alat) = %12.6f Bohr"%(alat)

    # volume 
    if vol is None and "unit-cell volume" in line:
      data = line.split() 
      vol = float(data[3])
      print "unit-cell volume = %12.6f Bohr^3"%(vol) 
    
    # natom 
    if natom == 0 and "number of atoms/cell" in line:
      info = line.split()
      natom = int(info[4]) 
      print "the number of atoms/cell (natom) =%4d"%(natom) 

    # latt_vectors
    if lat_vec == None and "crystal axes:" in line:
      print "Lattice Vectors for the initial structure:"
      lat_vec = []
      for i in range(3):
        data = ifile.readline().split()
        (x,y,z) = ( float(data[3]), float(data[4]), float(data[5]) )
        print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,x,y,z)
        lat_vec.append( [x,y,z] )
     
    if "Begin final coordinates" in line: break 

  if found_none: 
    print "ERROR: no info on the optimized structure has been found!"
    sys.exit(1)      

  print ""

  line = ifile.readline()

  if "new unit-cell volume" in line:
    print "Read output from a 'vc-relax' calculation"
    irelax = 1 
    data = line.split()
    vol_new = float(data[4]) 
    print "  new volume = %12.6f Bohr^3"%(vol_new)

    # get lattice vectors 
    line = ifile.readline()  # skip a blank line 

    line = ifile.readline()  # "CELL_PARAMETERS ..." 
    print "  lattice vectors in terms of the old lattice parameter(alat=%12.6f Bohr):"%(alat) 
    lat_vec_old = lat_vec[:] 
    lat_vec = [] 
    for i in range(3):
      data = ifile.readline().split()
      (x,y,z) = (float(data[0]), float(data[1]), float(data[2]) )
      lat_vec.append([x,y,z]) 
      print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,x,y,z)

    # print lattice vectors in the new alatt 
    alat_new = alat*(vol_new/vol)**(1.0/3)
    print "  lattice vectors in terms of the new lattice parameter(alat=%12.6f Bohr):"%(alat_new)
    for i in range(3):
      xyz=lat_vec[i][:]
      for j in range(3):
        xyz[j] *= alat/alat_new
      print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,xyz[0],xyz[1],xyz[2])

    line = ifile.readline()  # skip the blank line 
  else:
    irelax = 0 
    print "Read output from a 'relax' calculation"
  
  line = ifile.readline().strip()  #"ATOMIC_POSITIONS (...)"
  if   "alat"     in line:  mol_type = 'alat'
  elif "bohr"     in line:  mol_type = 'bohr' 
  elif "angstrom" in line:  mol_type = 'angstrom'
  elif "crystal"  in line:  mol_type = 'crystal' 
  mol = []   
  print line 
  for ia in range(natom):
    line = ifile.readline()
    data = line.split()
    atom = data[0]
    (x,y,z) = (float(data[1]), float(data[2]), float(data[3])) 
    print "%-6s %12.6f %12.6f %12.6f"%(atom,x,y,z)
    mol.append([atom,[x,y,z]])

  ifile.close()
  # Up to now, all information has been read, the following will manipulate the information a little bit before return 
  # the convention is that 
  #  latt_vec are returned in units of angstrom
  #  basis return atomic positions in terms of lattice vectors 
  #  mol  returns positions of atoms in units of angstrom 

  if mol_type == "alat" and irelax ==1 :  # vc-relax 
    # get the ratio between new and old lattice parameters
    alen_old = 0.0
    alen = 0.0
    for i in range(3):
      alen_old += lat_vec_old[0][i]**2
      alen += lat_vec[0][i]**2

    alat_old = alat
    alat = alat_old * sqrt(alen/alen_old)

    print "\n  new lattice vectors in terms of new lattice parameter (alat=%12.6f):"%(alat)
    for i in range(3):
      for j in range(3): lat_vec[i][j] *= alat_old/alat

      print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,lat_vec[i][0],lat_vec[i][1],lat_vec[i][2])

    print "\natomic positions in terms of new lattice parameter:"
    for ia in range(natom):
      for i in range(3):
         mol[ia][1][i] *= alat_old/alat
      atom = mol[ia][0]
      (x,y,z) = (mol[ia][1][0],mol[ia][1][1],mol[ia][1][2])
      print "%-6s %12.6f %12.6f %12.6f"%(atom,x,y,z)
      
  for iv in range(3):
    for i in range(3):
      lat_vec[iv][i] *= alat

  # convert alat and lat_vecs in units of angstrom instead of bohr
  alat *= Bohr2Ang
  for iv in range(3):
    for i in range(3):
      lat_vec[iv][i] *= Bohr2Ang

  for ia in range(natom):
    if mol_type == 'alat': 
      for i in range(3): mol[ia][1][i] *= alat

    elif mol_type == 'bohr': 
      for i in range(3): mol[ia][1][i] *= Bohr2Ang

    elif mol_type =='crystal':
      xyz_int = mol[ia][1][:]
      mol[ia][1][:] = f_Intern_to_Cart(xyz_int,lat_vec) 

  if iop_out == 0: 
    basis = f_Latt_Cart_to_Crys(mol,lat_vec)
    return basis,lat_vec,ibrav 
  else:
    return mol
  
def f_qe_write_struct(prefix,mol,latt=None,ibrav=0,pos_type='crys'):
  """ 
  Using the input molecular (and lattice constants and type ) to generate 
  an input file for PWscf with a few default options.

  When latt is defined, pos_type indidates the type of coodinates stored in "mol"
    'angs'  -- cartesian coordinates in angstrom
    'alat'  -- cartesian coordinates scaled by a (==latt[0])
    'crys'  -- use the crystal coodinates in the QE convention  
  """
  sname = "f_QE_Write_Struct"

  print sname+": create QE input for given structure \n!!!  the default options is most likely not appropriate !!! "

  name = prefix.strip()+"-pw.in"
  ofile = open(name,'w')
  print "Create new PWscf inpup file: " + name 

  # check how many types of atoms are present in the molecule 
  nat = len(mol)
  ntyp,nat_typ,species,sp_index = f_Check_Species(mol) 

  outdir = './out'
  pspdir = os.getenv("ESPRESSO_PSEUDO",'./') 
 
  # default input for control block 
  ofile.write("""&control
  calculation = 'scf'
  restart_mode = 'from_scratch'
""")
  ofile.write("  outdir = '"+ outdir+"'\n") 
  ofile.write("  pseudo_dir = '"+ pspdir+"'\n") 
  ofile.write("  prefix = '"+ prefix+"'\n") 
  ofile.write("/\n\n")

  ofile.write("\n&system\n")
  ofile.write("  nat = %d, ntyp = %d\n"%(nat,ntyp)) 

  xyz_shift = [ 0.0, 0.0, 0.0 ]  

  if latt is None : # finite systems 

    sc = f_Struct_Size(mol)
    ofile.write("  ibrav = %d\n"%(8))
    ofile.write("  celldm(1) = %12.6f\n"%(sc[0]*Ang2Bohr))
    ofile.write("  celldm(2) = %12.6f\n"%(sc[1]/sc[0]))
    ofile.write("  celldm(3) = %12.6f\n"%(sc[2]/sc[0]))

  else:           # periodic systems 
    print "  Lattice Constants:"
    print "\t    a =%12.6f,    b =%12.6f,     c=%12.6f"%(latt[0],latt[1],latt[2])
    print "\talpha =%12.6f, beta =%12.6f, gamma=%12.6f"%(latt[3],latt[4],latt[5])
    a = latt[0]
    boa=latt[1]/latt[0]
    coa=latt[2]/latt[0]
    cos_bc= cos(latt[3]*Deg2Rad)
    cos_ac= cos(latt[4]*Deg2Rad) 
    cos_ab= cos(latt[5]*Deg2Rad) 
    
    ofile.write("  ibrav = %d\n"%(ibrav))
    ofile.write(  "  celldm(1) = %12.6f\n"%(a*Ang2Bohr))
    if ibrav == 4 or ibrav == 6 or ibrav == 7:  # hex or trigonal P
      ofile.write("  celldm(3) = %12.6f\n"%(coa))
    elif ibrav == 5:  # triangle R
      ofile.write("  celldm(4) = %12.6f\n"%(cos_bc))
    else:
      ofile.write("  celldm(2) = %12.6f\n  celldm(3) = %12.6f\n"%(boa,coa))
      if ibrav == 12 or ibrav == 13:  # Monoclinic
        ofile.write("  celldm(4) = %12.6f\n"%(cos_ab))
      if ibrav == 14: 
        ofile.write("  celldm(4) = %12.6f\n  celldm(5) = %12.6f\n  celldm(6) = %12.6f\n"%(cos_bc,cos_ac,cos_ab))

  # default input for the block "system"
  ofile.write("""
  ecutwfc =  60.0  
  ecutrho = 480.0 

!* uncomment the following for metallic systems 
! occupations = 'smearing'
! smearing = 'gauss'
! degauss = 0.005

/
""") 

  # block electrons 
  ofile.write("""&electrons
  diagonalization='cg'
  mixing_mode = 'plain'
  mixing_beta = 0.7
  conv_thr =  1.0d-6
/
""")

  # block ions  
  ofile.write("""&ions
  ion_dynamics = 'bfgs'
  ion_positions = 'default'
/
""")

  # block cell
  ofile.write("""&cell

/
""")
  
  # card ATOMIC_SPECIES 
  ofile.write("ATOMIC_SPECIES\n")
  for i in range(ntyp):
    atom= species[i]
    mass = f_Element_Symbol_to_Mass(atom) 
    psp = atom + qe_default_psp.strip() + ".UPF"
    ofile.write("  %-4s %6.2f %s\n"%(atom,mass,psp)) 
 
  # card ATOMIC_POSITIONS

  scale = [1.0,1.0,1.0] 
  if latt is None: 
    ofile.write("ATOMIC_POSITIONS angstrom\n")
  else:
    if pos_type == 'crys':
      ofile.write("ATOMIC_POSITIONS crystal\n")
    elif pos_type == 'alat':
      ofile.write("ATOMIC_POSITIONS alat\n")
    else:
      ofile.write("ATOMIC_POSITIONS angstrom\n")

  for i in range(nat):
    atom=mol[i][0]
    x = mol[i][1][0]
    y = mol[i][1][1]
    z = mol[i][1][2] 
    ofile.write("  %-4s %12.6f %12.6f %12.6f\n"%(atom,x,y,z))

  # card K_POINTS
  ofile.write("K_POINTS automatic\n")
  if latt is None:
    ofile.write(" 1 1 1 0 0 0 \n")
  else:
    # estimate the number of k-points needed 
    nk=[1,1,1]
    for i in range(3):
      nk[i] = int(kgrid_cutoff/latt[i])+1
    ofile.write("%2d %2d %2d  0 0 0  \n"%(nk[0],nk[1],nk[2]))
  
  ofile.close()  

