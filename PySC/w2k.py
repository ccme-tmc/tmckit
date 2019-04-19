#!/usr/bin/env python
from math import *
from chem_utils import *
from io_utils import * 
#
# This file contains subroutines that are used to read or write wien2k struct file 
#  currently the following subroutines are provided:
#   W2k__read/write_struct

def W2k_struct_c2p(bas,latt_type):
  """
  Convert the basis in terms of the conventional unit cell vectors to that 
  in terms of the primitive cell vectors (in the wien2k convention)
  """
  nat = len(bas)
  bas_new = []
  for i in range(nat):
    (x1,x2,x3) = (bas[i][1][0], bas[i][1][1], bas[i][1][2])

    xyz = [x1,x2,x3]
    if latt_type == 'F':
      xyz[0] = - x1 + x2 + x3
      xyz[1] =   x1 - x2 + x3
      xyz[2] =   x1 + x2 - x3
    elif latt_type == 'B':
      xyz[0] = x2 + x3
      xyz[1] = x1 + x3
      xyz[2] = x1 + x2

    for j in range(3):
      if xyz[j] < 0.0: xyz[j] += 1.0
      if xyz[j] > 1.0: xyz[j] -= 1.0

    bas_new.append([bas[i][0],xyz])

  return bas_new

def W2k_read_struct(name='',mode=0,debug=False):
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
  case_name  = W2k_check_name(name)
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
  r0rmt = [] 
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
      r0rmt.append([r0,rmt]) 

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
    return latt_type,latt,basis,r0rmt
  elif mode == 1:
    return latt_type,latt,basis_eq
  else:
    return latt_type,latt,basis, misc

def W2k_write_struct(file,latt_type,latt,basis,misc=None,debug=False):
  """
  Write the structure in the wien2k struct format in the simple form
  """
  out_name = file.strip()+".struct"

  bak_name = out_name+"_bak"

  if debug: 
    print "Write the structure into " + out_name

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
    sym = [ [ [[1,0,0],[0,1,0],[0,0,1]],[0,0,0]] ]
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
      misc_atoms.append([atom,-iat-1,mult,isplit,npt,r0,rmt,znuc,rot])

  nat_neq = len(misc_atoms)

  ofile.write(title+"\n")
  ofile.write("%-4s%-23s%3d\n"%(latt_type,'LATTICE,NONEQUIV.ATOMS:',nat_neq))
  ofile.write("%-13s%4s\n"%('MODE OF CALC=',rela_flag))
  latt_new = latt[:]

  for i in range(3,6):
    if abs( latt_new[i] - 90.0 ) < 1.e-3:
      latt_new[i] = 90.0
    if abs( latt_new[i] - 120.0) < 1.e-3: 
      latt_new[i] = 120.0
   
  (a,b,c,af,bt,gm) = latt_new

  ofile.write("%10.6f%10.6f%10.6f%10.6f%10.6f%10.6f\n"%(a,b,c,af,bt,gm))
  iat = 0
  for i_neq in range(nat_neq):
    xyz = basis[iat][1][:]; iat += 1
    (atom_name,atom_index,mult,isplit,npt,r0,rmt,znuc,rot) = misc_atoms[i_neq]

    for i in range(3):
      if   abs(xyz[i]) < 1.e-4 or abs( xyz[i] - 1.0 ) < 1.e-4:
        xyz[i] = 0.0
      elif abs(xyz[i] - 1.0/3) < 1.e-4:
        xyz[i] = 1.0/3
      elif abs(xyz[i] - 0.5 )  < 1.e-4: 
        xyz[i] = 0.5 
      elif abs(xyz[i] - 2.0/3) < 1.e-4:
        xyz[i] = 2.0/3
      elif abs(xyz[i] - 5.0/6) < 1.e-4:
        xyz[i] = 5.0/6
      elif abs(xyz[i] - 1.0/6) < 1.e-4:
        xyz[i] = 1.0/6

    (x,y,z) = xyz
 
    ofile.write("ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n"%(atom_index,x,y,z))
    ofile.write("          MULT=%2d          ISPLIT=%2d\n"%(mult,isplit))
    if mult > 1:
      for i_eq in range(1,mult):
        (x,y,z) = basis[iat][1]; iat += 1
        ofile.write("ATOM%4d: X=%10.8f Y=%10.8f Z=%10.8f\n"%(atom_index,x,y,z))

    ofile.write("%-10s NPT=%5d  R0=%10.8f RMT=%10.4f   Z:%5.1f\n"%(atom_name,npt,r0,rmt,znuc))
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

def W2k_check_name(name=''):
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

