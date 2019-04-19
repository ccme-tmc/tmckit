#!/usr/bin/env python
from math import *
from chem_utils import *
from io_utils import * 
from list_utils import *

def Struct_display(latt_type,bas,alatt=None,latt_vec=None): 

  natom = len(bas) 
  if latt_type == 'mol':
    print "\nMolecular structure (in cartesian coordinates):"
    for i in range(natom):
      print "\t%4s(%12.6f,%12.6f,%12.6f)"%(bas[i][0],bas[i][1][0],bas[i][1][1],bas[i][1][2])
  else: 

    print "\nCrystal structure:"
    print "  Lattice Type:     ",latt_type
    print "  Lattice Constants:\n\t a=%10.4f, b=%10.4f, c=%10.4f\n\t alpha=%6.2f, beta =%6.2f, gamma=%6.2f"%(latt[0],latt[1],latt[2],latt[3],latt[4],latt[5])
    print "  Lattice Vectors:"
    for i in range(3):
      print "\ta(%1d)=(%12.6f,%12.6f,%12.6f)"%(i,latt_vec[i][0],latt_vec[i][1],latt_vec[i][2])

    print "  Internal coordinates:"
    for i in range(natom):
      print "\t%4s(%12.6f,%12.6f,%12.6f)"%(bas[i][0],bas[i][1][0],bas[i][1][1],bas[i][1][2])


def Struct_substitute_atom(bas,atom_old,atom_new,scheme):
  """
  Create the new structure by replacing atom_old with the atom_new in terms of "scheme"
  "scheme" is an array of integers indicating the indices of atoms to be replaced  
   if atom_new == 'Vac', a vacancy is created 
  """
  nat = len(bas)
  nsb = len(scheme)

  ia_sub = 0
  bas_new = []
  for ia in range(nat):
    at = bas[ia][0]
    xyz = bas[ia][1][:]
    if bas[ia][0] == atom_old:
      ia_sub += 1
      if ia_sub in scheme:
        at = atom_new
    if at != 'Vac': bas_new.append([at,xyz])

  return bas_new


def Struct_parse_info(crys_info):
  """
  Extract crystal structure from a string 
  """
  info = crys_info.split(':')
  ctype = info[0]
  bas = []
  if ctype == 'diamond' or ctype == 'diam':     # diamond structure 
    atom = info[1]
    alatt = float(info[2]) 
    latt_vec =[[0.0,0.5,0.5],\
               [0.5,0.0,0.5],\
               [0.5,0.5,0.0]]
    bas=[[atom,[0.125,0.125,0.125]],\
         [atom,[0.875,0.875,0.875]]]

  elif ctype == 'fcc':       # fcc metal (single-element) 
    atom = info[1]
    alatt = float(info[2])
    latt_vec =[[0.0,0.5,0.5],\
               [0.5,0.0,0.5],\
               [0.5,0.5,0.0]]
    bas=[[atom,[0.0, 0.0, 0.0]]]

  elif ctype == 'hcp':       # hcp metal (single-element) 
    atom = info[1]
    tmp = info[2].split(',') 
    alatt = float(tmp[0])
    c = float(tmp[1])
    
    latt_vec =[[1.0,0.0,0.0],\
               [0.5,0.0,0.5],\
               [0.5,0.5,0.0]]
    bas=[[atom,[0.0, 0.0, 0.0]]]

  elif ctype == 'zb' or ctype == 'zinc-blende':
    atoms = info[1].split(',') 
    alatt = float(info[2])

    latt_vec =[[0.0,0.5,0.5],\
               [0.5,0.0,0.5],\
               [0.5,0.5,0.0]]
    bas=[[atoms[0],[0.125, 0.125, 0.125]], \
         [atoms[1],[0.875, 0.875, 0.875]]]

  elif ctype == 'NaCl':
    atoms = info[1].split(',')
    alatt = float(info[2])

    latt_vec =[[0.0,0.5,0.5],\
               [0.5,0.0,0.5],\
               [0.5,0.5,0.0]]

    bas=[[atoms[0],[0.0, 0.0, 0.0]], \
         [atoms[1],[0.5, 0.5, 0.5]]]

  elif ctype == 'perovskite' or ctype == 'ABO3': 
    atoms = info[1].split(',')
    alatt = float(info[2])

    latt_vec =[[1.0,0.0,0.0],\
               [0.0,1.0,0.0],\
               [0.0,0.0,1.0]]

    bas=[[atoms[0],[0.0, 0.0, 0.0]], \
         [atoms[1],[0.5, 0.5, 0.5]], \
         [atoms[2],[0.0, 0.5, 0.5]], \
         [atoms[2],[0.5, 0.0, 0.5]], \
         [atoms[2],[0.5, 0.5, 0.0]], \
        ]

  elif ctype == 'FeS2' or ctype == 'pyrite': 
    atoms = info[1].split(',')
    tmp = info[2].split(',') 

    alatt = float(tmp[0])
    u = float(tmp[1])
    h = 0.5 
  
    latt_vec =[[1.0,0.0,0.0],\
               [0.0,1.0,0.0],\
               [0.0,0.0,1.0]]
    bas=[\
         [ atoms[0], [0.0, 0.0, 0.0] ], \
         [ atoms[0], [0.0, 0.5, 0.5] ], \
         [ atoms[0], [0.5, 0.0, 0.5] ], \
         [ atoms[0], [0.5, 0.5, 0.0] ], \
         [ atoms[1], [  u,   u,   u] ], \
         [ atoms[1], [ -u,  -u,  -u] ], \
         [ atoms[1], [u+h, h-u,  -u] ], \
         [ atoms[1], [h-u, u-h,   u] ], \
         [ atoms[1], [ -u, u+h, h-u] ], \
         [ atoms[1], [  u, h-u, u-h] ], \
         [ atoms[1], [h-u,  -u, u+h] ], \
         [ atoms[1], [u-h,   u, h-u] ], \
        ]
  else:
    print "ERROR: unsupported structure type "+ctype 
    return None 

  return alatt,latt_vec,bas
  

def f_Struct_Center(mol,iop=0):
  """
  Return the center of the molecule 
  """
  natom = len(mol) 
  if iop == 0: ## find geometrical center of the molecule 
    xyz_min = mol[0][1][:]
    xyz_max = mol[0][1][:]
    for ia in range(1,natom):
      xyz = mol[ia][1][:]
      for i in range(3):
        if xyz[i] < xyz_min[i]: xyz_min[i] = xyz[i]
        if xyz[i] > xyz_max[i]: xyz_max[i] = xyz[i]

    mol_center  = [0.0,0.0,0.0]
    for i in range(3):
      mol_center[i] = (xyz_max[i] + xyz_min[i])/2.0

  return mol_center 

def f_Struct_Size(mol, edge=4.0, ratio=1.0):
  """
  Return the size of the molecule, this is estimated by calculating 
  the spatial scope of the molecular coordinates plus some additional extension at the outer edge
  (edge)

  """
  xyz_min = mol[0][1][:]
  xyz_max = mol[0][1][:]

  ## find the maximum and minimum of x_i,y_i and z_i
  for ia in range(1,len(mol)):
    xyz = mol[ia][1][:]
    for i in range(3):
      if xyz[i] < xyz_min[i]: xyz_min[i] = xyz[i]
      if xyz[i] > xyz_max[i]: xyz_max[i] = xyz[i]

  xyz_size    = [0.0,0.0,0.0]
  for i in range(3):
    xyz_size[i] = (xyz_max[i] - xyz_min[i] + 2.0*edge )*ratio
  return xyz_size

def f_Struct_mol2sc(mol,iop_sc=0,sc_info=None):
  """
  Create a supercell for a given molecular structure 
  """
  mol_center = f_Struct_Center(mol)

  # first set the supercell size
  if iop_sc == 1 :
    sc_size = sc_info[0:2]
  elif iop_sc == 2:
    sc_size = f_Struct_Size(mol, edge=sc_info[0], ratio=1.0)
  elif iop_sc == 3:
    sc_size = f_Struct_Size(mol, edge=0.0, ratio=sc_info[0])
  else:
    sc_size = f_Struct_Size(mol) 

  ## set lattive vectors and constants
  latt_vec=[[sc_size[0],0.0,0.0], \
            [0.0,sc_size[1],0.0], \
            [0.0,0.0,sc_size[2]]]
  latt_type = 'P'
  latt = [ sc_size[0],sc_size[1],sc_size[2], 90.0, 90.0, 90.0 ]

  ## set basis in internal coordinates
  basis=[]
  for ia in range(len(mol)):
    atom=mol[ia][0]
    xyz=[]
    for i in range(3):
      xyz.append( (mol[ia][1][i]-mol_center[i])/sc_size[i] + 0.50)

    basis.append([atom,xyz])

  return latt_type,latt,basis 

def f_Check_Species(basis):
  """
  Check the number of different atoms (species) in a molecule 
  return nat_sp(1..nsp), species(1..nsp), and sp_index(1..nat) 
  """
  nat = len(basis)
  sp_index = []
  for iat in range(nat):
    if iat == 0:
      species = [ basis[iat][0] ]
      nsp = 1
      sp_index.append(1)
      continue

    l_new_sp = True
    for isp in range( nsp ):
      if species[isp] == basis[iat][0]:
        l_new_sp = False
        sp_index.append(isp+1)
        break

    if l_new_sp:
      nsp += 1
      species.append( basis[iat][0] )
      sp_index.append(nsp)

  nat_sp = f_List_New(0,nsp) 
  for iat in range(nat):
    for isp in range(nsp):
      if sp_index[iat] == isp+1: nat_sp[isp] += 1
  
  return nsp,nat_sp,species,sp_index

def f_Check_Formula(mol):
 # count the number of atoms of each type
  for iat in range(len(mol)):
    if iat==0:
      spec = [ mol[0][0] ]
      nat_spec =[1]
      continue

    l_new_sp = True
    for isp in range( len(spec) ):
      if spec[isp] == mol[iat][0]:
        l_new_sp = False
        nat_spec[isp] += 1
        break

    if l_new_sp:
      spec.append( mol[iat][0] )
      nat_spec.append(1)
  chem_form=''
  for isp in range(len(spec)):
     chem_form += "%s_%d "%(spec[isp].strip(),nat_spec[isp])
  return chem_form

def f_Sort_Struct(mol,ic):
  """
  Sort the atoms in a molecule in terms of one of their coordinates
  """
  nat = len(mol)
  if ic == 0 or abs(ic) > 3: 
    if abs(ic) > 3 : print "f_Sort_Struct: ic = ",ic," -- molcule not sorted "  
    return  

  for i in range(nat):
    for j in range(i+1,nat): 
      atm_i = mol[i][0]
      xyz_i = mol[i][1][:]
      xyz_j = mol[j][1][:]
      atm_j = mol[j][0]
      if ic > 0: 
        if xyz_j[ic-1] < xyz_i[ic-1]: 
          mol[i][0] = atm_j  
          mol[i][1] = xyz_j[:]
          mol[j][0] = atm_i
          mol[j][1] = xyz_i[:] 
      else:
        if xyz_j[-ic-1] > xyz_i[-ic-1]:
          mol[i][0] = atm_j  
          mol[i][1] = xyz_j[:]
          mol[j][0] = atm_i
          mol[j][1] = xyz_i[:]

def f_Struct_Permute(n,mol,latt=None):
  """
  Permute coordinates in a molecule by once(n=1) or twice(n=2)
  n==1: (x,y,z) -> (z,x,y) 
  n==2: (x,y,z) -> (y,z,x) 
  """
  nat=len(mol)
  for i in range(nat):
    x=mol[i][1][0]
    y=mol[i][1][1]
    z=mol[i][1][2]
    if n == 1: 
      mol[i][1][0] = z 
      mol[i][1][1] = x
      mol[i][1][2] = y 
    elif n == 2:
      mol[i][1][0] = y
      mol[i][1][1] = z
      mol[i][1][2] = x 
  if not latt is None:
    (a,b,c) = (latt[0],latt[1],latt[2]) 
    (bc,ca,ab) = (latt[3],latt[4],latt[5]) 
    if n == 1:
      latt[0] = c 
      latt[1] = a 
      latt[2] = b 
      latt[3] = ab 
      latt[4] = bc 
      latt[5] = ca 
    
    elif n == 2:
      latt[0] = b
      latt[1] = c
      latt[2] = a
      latt[3] = ca
      latt[4] = ab
      latt[5] = bc

def f_NN_Generator(r0,atom,bond_len,mode):
  """
  This function generates nearest neighbouring coordniates around r0
  """
  myname = "f_NN_Generator"
  
  xyznew=[]
  s3d2 = sqrt(3.0)/2.0
  s3d6 = sqrt(3.0)/6.0 
  s3d3 = sqrt(3.0)/3.0
  s2d3 = sqrt(2.0/3.0) 
  if mode == "fcc" or mode =="hcp":
    xyznew.append(   [ 1.0,  0.0,  0.0] )    # 1
    xyznew.append(   [ 0.5, s3d2,  0.0] )    # 2
    xyznew.append(   [-0.5, s3d2,  0.0] )    # 3
    xyznew.append(   [-1.0,  0.0,  0.0] )    # 4
    xyznew.append(   [-0.5,-s3d2,  0.0] )    # 5
    xyznew.append(   [ 0.5,-s3d2,  0.0] )    # 6
 
    xyznew.append(   [ 0.5, s3d6, s2d3] )    # 7
    xyznew.append(   [-0.5, s3d6, s2d3] )    # 8
    xyznew.append(   [ 0.0,-s3d3, s2d3] )    # 9 
    if mode == "hcp": 
      xyznew.append( [ 0.5, s3d6,-s2d3] )    # 10
      xyznew.append( [-0.5, s3d6,-s2d3] )    # 11
      xyznew.append( [ 0.0,-s3d3,-s2d3] )    # 12 
    else:
      xyznew.append( [ 0.0, s3d3,-s2d3] )    # 10
      xyznew.append( [-0.5,-s3d6,-s2d3] )    # 11
      xyznew.append( [ 0.5,-s3d6,-s2d3] )    # 12 
  else:
    print "ERROR in " + myname + ": unsupported mode=%s" %(mode)
    sys.exit(1) 
  
  for iat in range(len(xyznew)):
    for i in range(3): xyznew[iat][i] = r0[i] + xyznew[iat][i]*bond_len

  return xyznew 

def f_Cluster_Tetra(atom,acub,nord):
  """
  Create a tetrahedral cluster (Au_20 type) 
  """
  alat=[ [0.0,0.5,0.5],\
         [0.5,0.0,0.5],\
         [0.5,0.5,0.0] ]
  for iv in range(3):
    for i in range(3): alat[iv][i] *= acub 
  
  xyz=[0.0,0.0,0.0]
  mol=[]
  for iz in range(nord+1):
    for iy in range(nord-iz+1):
      for ix in range(nord-iz-iy+1):
        for i in range(3): xyz[i] = ix*alat[0][i]+iy*alat[1][i]+iz*alat[2][i]
        mol.append([atom,xyz[:]])
  return mol

def f_Triangle_Area(a,b,c):
  p = (a+b+c)/2.0
  area = sqrt(p*(p-a)*(p-b)*(p-c))
  return area

def f_Distance(p1,p2):
  """
  Calculate the distance between points 
  """
  dist = 0.0 
  for i in range(3): dist += (p2[i]-p1[i])**2
  dist = sqrt(dist)
  return dist 

def f_NN_Distance(atoms_xyz):
  """
  Return the nearest neighbouring distance between atoms in the lattice
  """
  nat = len(atoms_xyz)
  d_min = 1.e5
  for i in range(nat):
    r_i = atoms_xyz[i][1][:]
    for j in range(nat):
      if i>=j: continue
      r_j = atoms_xyz[j][1][:]
      d_ij= f_Distance(r_i, r_j)
      if d_ij < d_min : d_min = d_ij
  return d_min

def f_Add_Atom(mol,atom,iop,args):
  """
  Add a new atom to the molecular controlled by "option"
  the content of args varies with "iop"
    iop: 
       0 --  directly with given coordinates, 
           args -> [x,y,z]  
       1 --  on top of an atom along ix-th direction, 
           args -> [i1,ix,rbond]   
       2 --  along the direction determined by two atoms, 
           args -> [i1,i2,rbond]
       3 --  on bridge of two atoms (along ix-th axis)"
           args -> [i1,i2,rbond, ix ]
       4 -- on bridge of two atoms on the is-th plane
           args -> [i1,i2,rbond, is] 

  """
  xyz=[0.0, 0.0, 0.0]
  if   iop == 0: 
    xyz=args[:]

  elif iop == 1: 
    i1 = args[0]
    ix = args[1] - 1
    rbond = args[2] 
    xyz = mol[i1][1][:]
    xyz[ix] += rbond
    
  elif  iop == 2: 
    i1  = args[0]
    i2  = args[1]
    rbond = args[2] 

    r12 = f_Distance(mol[i2][1],mol[i1][1])
    for i in range(3): xyz[i] = mol[i2][1][i]+(mol[i2][1][i]-mol[i1][1][i])*rbond/r12

  elif iop == 3: 
    i1 = args[0]; i2 = args[1]; rbond = args[2]; ix = args[3] - 1

    a = f_Distance( mol[i1][1], mol[i2][1])
    dz = sqrt(dbond**2-a**2/4.)

    for i in range(3): xyz[i] = (mol[i1][1][i] + mol[i2][1][i])*0.5
    xyz[ix] += dz

  elif iop == 5:
    i1 = args[0]; i2 = args[1]; i3 = args[2]; rbond= args[3]

    a = f_Distance( mol[i1][1], mol[i2][1] )   
    dz=(4.0/sqrt(3.)/a)*f_Triangle_Area(rbond,sqrt(3.0)/2*a,sqrt(rbond**2-a**2/4.0))

    for i in range(3): xyz[i]= (mol[i1][1][i]+mol[i2][1][i]+mol[i3][1][i])/3.0
    xyz[2] += dz 

  mol.append([atom,xyz])

