#!/usr/bin/env python
import sys,os,shutil
from math import *
from chem_utils import *
from struct_utils import * 
mod_name = 'vasp.py'

def vasp_read_poscar(fname='POSCAR'):
  """
  Read structural information from POSCAR
  """
  ifile = open(fname,'r') 
  titile = ifile.readline()
  alatt = float(ifile.readline())
  latt_vec=[]
  for i in range(3):
    line_s = ifile.readline().split()
    x = float(line_s[0])
    y = float(line_s[1])
    z = float(line_s[2]) 
    latt_vec.append([x,y,z])
  
  specs = ifile.readline().split() # the list of species 
  if specs[0].isdigit():
    print "ERROR in vasp_read_poscar: The file "+fname+" does not contain species information!" 
    ifile.close()
    sys.exit(1)

  nspec = len(specs)   # the number of species 

  line_s = ifile.readline().split()

  nat = 0 
  nat_spec = []  # the number of atoms for each species 
  for i in range(nspec):
    nat_spec.append(int(line_s[i]))
    nat += nat_spec[i]
  
  line = ifile.readline()
  print line 
  if line[0] == 'S' or line[0] == 's':
    line = ifile.readline()
    print line 

  mode = line[0]
  basis = []
  for isp in range(nspec):
    atom = specs[isp]
    for ia in range(nat_spec[isp]):
      
      line_s = ifile.readline().split()
      xyz = []
      for i in range(3):
        xyz.append(float(line_s[i]))
      basis.append([atom,xyz])

  ifile.close()
  return alatt,latt_vec,specs,nat_spec,mode,basis 

def vasp_write_poscar(mode,basis,alatt,latt_vec,fname=''):
  print """
  Write structural information in the the format consistent with POSCAR 
    mode =0/1 -- direct/cartesian
    fname -- the name of the file  
  """ 
  
  nat = len(basis) 
  nsp,nat_sp,species,sp_index = f_Check_Species(basis)

  mol_name=''
  for i in range(nsp):
    mol_name += species[i]
    if nat_sp[i] > 1:
      mol_name += str(nat_sp[i])

  print "The compostion of the system: "+ mol_name 
  
  if fname=='':
    fn = mol_name+".POSCAR"
  else:
    fn = fname+".POSCAR"

  ofile = open(fn,'w') 
  if fname == '':
    info = mol_name 
  else:
    info = fname 

  ofile.write(info+"\n")
  
  # lattice vectors
  ofile.write("%12.6f\n"%(alatt))
  for ix in range(3):
    xyz=latt_vec[ix]
    ofile.write("%12.6f %12.6f %12.6f\n"%(xyz[0],xyz[1],xyz[2]))

  # number of atoms per species
  info1 = ''
  info2 = ''
  for isp in range(nsp):
    info1 +="%6s "%(species[isp])
    info2 +="%6d "%(nat_sp[isp])
  ofile.write(info1+"\n")
  ofile.write(info2+"\n")

  # the mode of coodinates 
  if mode == 0:
    ofile.write("direct\n")
  else:
    ofile.write("cart\n")

  for isp in range(nsp):
    for ia in range(nat):
      if sp_index[ia] == isp + 1: 
        xyz=basis[ia][1]
        ofile.write("%12.6f %12.6f %12.6f \t # %s \n"%(xyz[0],xyz[1],xyz[2],species[isp]))  

  ofile.close()
