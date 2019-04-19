#!/usr/bin/env python
from math import *
from numpy import *
from constants import *
from list_utils import *
from struct_utils import *

zero_dist = 1.e-10

def f_Latt_reset_axis(mode,latt,basis):
  """
  Reset the axis 
  """
  latt0 = latt[:]
  if mode == 1: # rotate x,y,z so that the original x axis becomes z axis
    latt[0] = latt0[1]
    latt[1] = latt0[2]
    latt[2] = latt0[0]
    latt[3] = latt0[4]
    latt[4] = latt0[5]
    latt[5] = latt0[3] 
    for ia in range(len(basis)):
      xyz0 = basis[ia][1][:]
      basis[ia][1][0] = xyz0[1]
      basis[ia][1][1] = xyz0[2]
      basis[ia][1][2] = xyz0[0]

  elif mode == 2: # rotate x,y,z so that the original x axis becomes z axis
    latt[0] = latt0[2]
    latt[1] = latt0[0]
    latt[2] = latt0[1]
    latt[3] = latt0[5]
    latt[4] = latt0[3]
    latt[5] = latt0[4]
    for ia in range(len(basis)):
      xyz0 = basis[ia][1][:]
      basis[ia][1][0] = xyz0[2]
      basis[ia][1][1] = xyz0[0]
      basis[ia][1][2] = xyz0[1]

  elif mode == -1:  # switch x and z axis
    latt[0] = latt0[2]
    latt[2] = latt0[0]
    latt[3] = latt0[5]
    latt[5] = latt0[3]
    for ia in range(len(basis)):
      xyz0 = basis[ia][1][:]
      basis[ia][1][0] = xyz0[2] 
      basis[ia][1][2] = xyz0[0]

  elif mode == -2:  # switch y and z axis
    latt[1] = latt0[2]
    latt[2] = latt0[1]
    latt[4] = latt0[5]
    latt[5] = latt0[4]
    for ia in range(len(basis)):
      xyz0 = basis[ia][1][:]
      basis[ia][1][1] = xyz0[2]
      basis[ia][1][2] = xyz0[1]
  else:
    print "WARNING: mode = %d in f_Latt_reset_axis is not supported ! "%(mode) 


def f_Latt_Get_ibrav(type,latt,basis):
  """
  Get the Bravais index of a crystal structure following the QuantumEspresso convention from 
  and in some cases the lattice vectors are permuted to be consisitent with 
  QE convention 
  """ 
  ibrav = 0 
  type_new = type 
  (a,b,c,bc,ca,ab) = latt
  if type == 'F':     # face centered
    if a==b and b==c: # cubic face-centered
      ibrav = 2
    else:             # Orthorhombic face-centered
      ibrav = 10
  elif type == 'B':   # body centered
    if a==b and b==c: # cubic body-centered
      ibrav = 3
    elif a==b:        # tetragonal body-centered
      ibrav = 7
    else:             # Orthorhombic body-centered
      ibrav = 11
  elif type == 'R':   # rhombohedral
    ibrav = 5
  elif type == 'H':
    ibrav = 4
  elif type == 'CXY' or type == 'CYZ' or type == 'CXZ': # base-centered
    if type == 'CYZ':
      print " Founnd a-base centered orthohombic lattice (latt_type== CYZ):"
      print "   -- permutate x,y,z twice so that it becomes c-base centerd:"
      f_Struct_Permute(2,basis,latt)
      (a,b,c,bc,ca,ab) = latt
      type_new = 'CXY' 
    elif type == 'CXZ':
      if ab == 90.0: 
        print " Found b-base centered orthohombic lattice (latt_type=CXZ):"
        print "  -- permutate x,y,z so that it is c-base centered"
        f_Struct_Permute(1,basis,latt)
        (a,b,c,bc,ca,ab) = latt
        type_new = 'CXY'
       
    if bc == 90.0 and ca == 90.0 and ab == 90.0: # orthorhombic
      ibrav = 9
    else:                                       # monoclinic
      ibrav = 13
  else:          # primitive lattice
    if bc==90.0 and ca == 90.0 and ab == 90.0:
      if a==b and b==c:  # simple cubic
        ibrav = 1
      elif a == b and b != c: # simple tetragonal
        ibrav = 6
      else:                   # orthombic
        ibrav = 8
    elif bc != 90.0 and ca != 90.0 and ab != 90.0: # triclinic
      ibrav = 14
    else:  # monoclinic
      if bc != 90.0 :
        print "Monoclinic lattice and b-c angle != 90: permute (x,y,z) axis twice!"
        f_Struct_Permute(2,basis,latt)
      elif ca != 90.0:
        print "Monoclinic lattice and c-a angle != 90: permute (x,y,z) axis !"
        f_Struct_Permute(1,basis,latt)
      ibrav = 12
  return ibrav,type_new

def f_Latt_Type(latt):
  """
  return lattice type for given lattice constants
  """
  if latt[3] == 90.0 and latt[4] == 90.0 and latt[5] == 90.0:
    if latt[0] == latt[1] and latt[0] == latt[2]:
      type = 'cub'
    elif latt[0] == latt[1]:
      type = 'tetra'
    else:
      type = 'ortho'

  elif latt[3] == 90.0 and latt[4] == 90.0 and latt[5] == 120.0:
    type = 'hex'

  elif ( latt[3] == 90.0 and latt[4] == 90.0 ) \
    or ( latt[3] == 90.0 and latt[5] == 90.0 ) \
    or ( latt[4] == 90.0 and latt[5] == 90.0 ) :
    type = 'mono'
  else:
    type = 'tric'

  return type

def f_Ortho_Latt(latt):
  """ 
  test whether the lattice is orthogonal 
  """
  if latt[3] == 90.0 and latt[4] == 90.0 and latt[5] == 90.0: return True
  else: return False

def f_Latt_Direction(dir,alatt,latt_type):
  """
  Convert the direction vector represented by lattice vectors to Cartesian coordinates and normalize it 
  """
  # set direction vector to extract the cluster
  dir_vec = [0.0, 0.0, 0.0]
  if f_Ortho_Latt(alatt):
    for i in range(3): dir_vec[i] = float(dir[i])*alatt[i]
  else:
    latt_vec= f_Latt_Vectors(latt_type,alatt)
    for i in range(3):
      for j in range(3): dir_vec[i] += float(dir[j])*latt_vec[j][i]

  norm=0.0
  for i in range(3): norm += dir_vec[i]**2
  norm = sqrt(norm)
  for i in range(3): dir_vec[i] /= norm
  return dir_vec

#
# Calculate cell volume from given lattice constants
#
def f_Cell_Volume(latt):
  """
  Calculate volume per unit cell from given lattice constants
  """
  a=latt[0];     b=latt[1];     c=latt[2]
  alf=latt[3]*Pi/180;   bet=latt[4]*Pi/180;   gam=latt[5]*Pi/180
  vol= a*b*c*sqrt(1.0+2.0*cos(alf)*cos(bet)*cos(gam)-(cos(alf))**2-(cos(bet))**2-(cos(gam))**2)
  return vol

#
# Reset lattice parameters a certain percentage
#
def f_Reset_Latt1D(option,val,latt):
  """
  Reset lattice parameters a certain percentage
  """

  if option == "vol":
    vol=f_Cell_Volume(latt)
    ratio = (val/vol)**(1.0/3);
    latt[0] *= ratio
    latt[1] *= ratio
    latt[2] *= ratio
  elif option == "a":
    latt[0] = val
  elif option == "b":
    latt[1] = val
  elif option == "c":
    latt[2] = val
  elif option == "coa":  # change c/a with fixed volume for tetragonal and hexagonal lattice (a=b)
    coa_new = val;
    a2c=latt[0]**2*latt[2]
    latt[0]=(a2c/coa_new)**(1.0/3)
    latt[1]=latt[0]
    latt[2]=latt[0]*coa_new
  elif option == "coa_o":  # change c/a with fixed volume, b/a and lattice angles ( for lattices other than tetragonal and hexagonal ones)
    coa_new = val
    boa = latt[1]/latt[0]
    abc = latt[0]*latt[1]*latt[2]
    latt[0] = (abc/(coa_new*boa))**(1.0/3)
    latt[1] = latt[0]*boa
    latt[2] = latt[0]*coa_new
  elif option == "boa_o":  # change b/a with fixed volume, c/a and lattice angles (for lattices other than tetragonal and hexagonal ones
    boa_new = val
    coa = latt[2]/latt[0]
    abc = latt[0]*latt[1]*latt[2]
    latt[0] = (abc/(coa*boa_new))**(1.0/3)
    latt[1] = latt[0]*boa_new
    latt[2] = latt[0]*coa
  else:
    print "ERROR when calling f_Reset_Latt1D: option " + option + " not supported"
    return 1
  return 0

def f_Latt_Distance(p1,p2,latt=None):
  """
  This function calculate the distance between points, represented by internal coordinates
  """
  if latt is None:
    dist = sqrt( (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2 )
  else: 
    dx = (p2[0]-p1[0])*latt[0]
    dy = (p2[1]-p1[1])*latt[1]
    dz = (p2[2]-p1[2])*latt[2]
    af = latt[3]*Deg2Rad
    bt = latt[4]*Deg2Rad
    gm = latt[5]*Deg2Rad
    d2 = dx*dx + dy*dy + dz*dz + 2.0*dx*dy*cos(gm)+ 2.0*dy*dz*cos(af)+2.0*dx*dz*cos(bt)
    dist = sqrt(d2) 

  return  dist

def f_Latt_NN_Distance(atoms_xyz,iop=0,latt_type=None,alatt=None):
  """
  Return the nearest neighbouring distance between atoms in the lattice
  """
  nat = len(atoms_xyz)
  d_min = 1.e5
  if iop == 0:
    for i in range(nat):
      r_i = atoms_xyz[i][1][:]
      for j in range(nat):
        if i>=j: continue 
        r_j = atoms_xyz[j][1][:]
        d_ij= f_Latt_Distance(r_i, r_j)
        if d_ij < d_min : d_min = d_ij 
  return d_min 

def f_Intern_to_Cart(xyz,latt):
  """ 
  This function convert internal coordinates to cartesian coordinates
  if iop == 0, latt contains lattice constants
  else, latt constains lattice vectors 
  """
  if len(latt) == 6: 
    a=latt[0];  b=latt[1];  c=latt[2]
    af=latt[3]*Deg2Rad; bt=latt[4]*Deg2Rad; gm=latt[5]*Deg2Rad
    vol= f_Cell_Volume(latt)
    x=xyz[0];   y=xyz[1];   z=xyz[2]

    xc = x*a+ y*b*cos(gm) + z*c*cos(bt)
    yc = y*b*sin(gm) - z*c*(cos(bt)*cos(gm)-cos(af))/sin(gm)
    zc = z*(vol/(a*b*sin(gm)))
    xyz_c = [ xc, yc, zc]
  else:
    xyz_c = [0.0,0.0,0.0]
    for i in range(3):
      for j in range(3): xyz_c[i] += xyz[j]*latt[j][i]

  return xyz_c

def f_Latt_Struct_I2C(bas,latt,iop):
  """ 
  Convert crystal structures in the internal coordinates into cartesian ones
  if iop == 0, latt contains lattice constants
  else, latt contains lattice vectors 
  """
  nat = len(bas)
  bas_c = []

  for iat in range(nat):
    atom = bas[iat][0]
    xyz = bas[iat][1][:]
    xyz_c = f_Intern_to_Cart(xyz,latt)
    bas_c.append( [atom, xyz_c] ) 
  return bas_c

def f_Latt_Cart_to_Crys(mol,lat_vec):
  """
  Convert Cartesian representation into crystal representation 
  """
  basis=[]
  nat=len(mol)
  
  lat_vec_inv = f_List_inv3(lat_vec) 

  for ia in range(nat):
    atom = mol[ia][0]
    xyz_cart = mol[ia][1][:] 
    xyz_crys = f_List_mv3( lat_vec_inv, xyz_cart) 
    basis.append([atom,xyz_crys]) 

  return basis 


def f_Latt_Vec2Const(vec,ibrav=0):
  """
  Extract lattice constants from lattice vectors
  """

  if ibrav == 1:  # simple cubic 
    a = vec[0][0]
    (b,c,af,bt,gm) = (a,a,90.0,90.0,90.0)  
  elif ibrav == 2: # cubic fcc
    a = 0.0
    for i in range(3):
      a += abs(vec[0][i])
    (b,c,af,bt,gm) = (a,a,90.0,90.0,90.0)

  elif ibrav == 3: #  bcc
    a = abs(vec[0][0])*2
    (b,c,af,bt,gm) = (a,a,90.0,90.0,90.0)

  elif ibrav == 4: # hexagonal or Trigonal P 
    a = 0.0
    for i in range(3):
      a += abs(vec[0][i])

    b = a
    c = vec[2][2]
    (af,bt,gm) = (90.0,90.0,120.0)

  elif ibrav == 6: #  Tetragonal P (st)
    a = vec[0][0]
    b = a
    c = vec[2][2]
    (af,bt,gm) = (90.0,90.0,90.0)

  elif ibrav ==7: # Tetragonal I (bct) 
    a = 0.0
    for i in range(2):
      a += abs(vec[0][i])
    b = a
    c = vec[2][2]
    (af,bt,gm) = (90.0,90.0,90.0)

  elif ibrav == 8: #  Orthorhombic P
    a = vec[0][0]
    b = vec[1][1]
    c = vec[2][2]
    (af,bt,gm) = (90.0,90.0,90.0)
  
  elif ibrav == 9 or ibrav == -9 : # Orthorhombic base-centered(bco) 
    a = abs(vec[0][0])*2
    b = abs(vec[1][1])*2
    c = vec[2][2]
    (af,bt,gm) = (90.0,90.0,90.0)

  elif ibrav == 10:   #  Orthorhombic face-centered
    (a,b,c) = 0.0 
    for i in range(3):
      a += vec[i][0]
      b += vec[i][1]
      c += vec[i][2]
    (af,bt,gm) = (90.0,90.0,90.0)

  elif ibrav == 11:   #  Orthorhombic body-centered
    (a,b,c) = 0.0
    for i in range(2):
      a += abs(vec[i][0])
      b += abs(vec[i][1])
      c += abs(vec[i][2])
    (af,bt,gm) = (90.0,90.0,90.0)    
  
  elif ibrav == 12: # Monoclinic P, unique axis c
    a = vec[0][0]
    b = sqrt(vec[1][0]**2 + vec[1][1]**2)
    c = vec[2][2]
    (af,bt) = (90.0,90.0)
    gm = atan(vec[1][1]/vec[1][0])*Rad2Deg

  elif ibrav == -12: # Monoclinic P, unique axis b
    a = vec[0][0]
    b = vec[1][1]
    c = sqrt(vec[2][0]**2 + vec[2][2]**2)
    af = 90.0
    bt = atan(vec[2][0]/vec[2][2])*Rad2Deg
    gm = 90.0 
  elif ibrav == 13:  #  Monoclinic base-centered 
    a = vec[0][0]*2
    b = sqrt(vec[1][0]**2 + vec[1][1]**2)
    c = vec[2][0]*2
    (af,bt) = (90.0,90.0)
    gm = atan(vec[1][1]/vec[1][0])*Rad2Deg

  else: 
    va = vec[0]
    vb = vec[1] 
    vc = vec[2]

    va_vb = 0.0
    va_vc = 0.0
    vb_vc = 0.0 
    for i in range(3):
      va_vb += va[i]*vb[i]
      va_vc += va[i]*vc[i]
      vb_vc += vb[i]*vc[i]
  
    a  = sqrt(va[0]**2+va[1]**2+va[2]**2)
    b  = sqrt(vb[0]**2+vb[1]**2+vb[2]**2)
    c  = sqrt(vc[0]**2+vc[1]**2+vc[2]**2)
    af = acos(vb_vc/(b*c))*Rad2Deg
    bt = acos(va_vc/(a*c))*Rad2Deg
    gm = acos(va_vb/(a*b))*Rad2Deg  
  
  
  return [a,b,c,af,bt,gm]

def f_Latt_Vectors(latt_type,alatt):
  """ 
  This function generate lattices vectors in terms of lattice type and lattice constants
  following the convention used in WIEN2k 
  """
  vec=f_List_New(f_List_New(0.0,3),3)

  a=alatt[0];   b=alatt[1];  c=alatt[2]
  af= alatt[3]*Deg2Rad; bt=alatt[4]*Deg2Rad; gm=alatt[5]*Deg2Rad
  sq32= sqrt(3.0)/2.0
  vcell = f_Cell_Volume(alatt)

  cosa=cos(af)
  cosb=cos(bt)
  cosg=cos(gm)   
  sing=sin(gm) 

  if latt_type == "P": 
    vec[0] = [ a, 0.0, 0.0]
    vec[1] = [ b*cosg, b*sing, 0]
    vec[2] = [ c*cosb, c*(cosa-cosb*cosg)/sing,c*sqrt(1+2*cosa*cosb*cosg-cosa**2-cosb**2-cosg**2)/sing]
  elif latt_type == "F":
    vec[0] = [ 0.0, b/2, c/2]
    vec[1] = [ a/2, 0.0, c/2]
    vec[2] = [ a/2, b/2, 0.0]
  elif latt_type == "B":
    vec[0] = [-a/2, b/2, c/2]
    vec[1] = [ a/2,-b/2, c/2]
    vec[2] = [ a/2, b/2,-c/2]
  elif latt_type == "H":
    vec[0] = [   a,   0.0, 0.0]
    vec[1] = [-a/2,sq32*a, 0.0]
#    vec[0] = [sq32*a,-a/2, 0.0]
#    vec[1] = [0.0   ,   a, 0.0]
    vec[2] = [0.0   , 0.0,   c]
  elif latt_type == 'R':
    tx=sqrt((1.0-cosa)/2.0)
    ty=sqrt((1.0-cosa)/6.0)
    tz=sqrt((1.0+2.0*cosa)/3.0)
    vec[0] = [a*tx,   -a*ty, tz*a]
    vec[1] = [0.0,   2*ty*a, tz*a]
    vec[2] = [-tx*a, -ty*a,  tz*a]

#    vec[0] = [a/sq32, -a/2, c/3]
#    vec[1] = [a/sq32, a/2,  c/3]
#    vec[2] = [-a/sqrt(3.0), 0.0, c/333] 
  elif latt_type == 'CXY': 
    vec[0] = [a/2,-b/2, 0]
    vec[1] = [a/2, b/2, 0]
    vec[2] = [  0,   0, c]
  elif latt_type == 'CYZ':
    vec[0] = [  a,   0,  0]
    vec[1] = [  0,-b/2,c/2]
    vec[2] = [  0, b/2,c/2] 
  elif latt_type == 'CXZ': 
#    vec[0] = [a*sin(gm)/2, a*cos(gm)/2, -c/2]
#    vec[1] = [0, b, 0]
#    vec[2] = [a*sin(gm)/2, a*cos(gm)/2,  c/2] 
    vec[0] = [a/2.0,  0.0, -c/2.0]
    vec[1] = [b*cosg, b*sing, 0.0]
    vec[2] = [a/2.0, 0.0,   c/2.0]
  else: 
    print "ERROR in f_Latt_Vectors: unsupported lattice type: ",latt_type
    sys.exit(1)
  
  return vec

def f_Charge_Center(charges):
  """
  This function return the charge center of a set of point charges (or atoms with apparent charges)
  """  
  nat = len(charges)
  center = [0.0,0.0,0.0]
  chtot = 0.0
  for iat in range(nat):
    xyz=charges[iat][2][:]
    ch = abs(charges[iat][1])
    center[0] += xyz[0]*ch
    center[1] += xyz[1]*ch
    center[2] += xyz[2]*ch
    chtot += ch

  center[0] /= chtot
  center[1] /= chtot
  center[2] /= chtot

  return center 

def f_Set_Range(rmax,a):
  nmax = int(rmax/a+0.49) 
  irange=[0]
  for i in range(1,nmax,1):
    irange.append(i)
    irange.append(-i)
  return irange

