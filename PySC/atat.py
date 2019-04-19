#!/usr/bin/env python
import sys,os,shutil
from math import *
from chem_utils import *
from struct_utils import * 
mod_name = 'atat.py'

def atat_read_latin(fname='lat.in'):
  """
  Read structural information from the file in the atat format
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

def atat_write_latin(mode,latt_type,latt,basis,fname='lat.in'):
  """
  Write structural information in the the format consistent with POSCAR 
    mode =0/1 --- the coordinate system is denoted by lattice constants/vectors 
    fname -- the name of the file  
  """ 
  nat = len(basis) 

  ofile = open(fname,'w') 

  if mode == 0:
    ofile.write("%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f \n"%(latt[0],latt[1],latt[2],latt[3],latt[4],latt[5]))
  else:
    for i in range(3):
      ofile.write("%12.6f %12.6f %12.6f \n"%(latt[i][0],latt[i][1],latt[i][2]))

  if latt_type == 'P' or latt_type == 'H':
    ofile.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n") 
  elif latt_type == 'F':
    ofile.write("0.0 0.5 0.5\n0.5 0.0 0.5\n0.5 0.5 0.0\n")
  elif latt_type == 'B':
    ofile.write("-0.5 0.5 0.5\n0.5 -0.5 0.5\n0.5 0.5 -0.5\n")
  else:
    print "atat_write_latin: unsupported value for latt_type=",latt_type
    sys.exit(1) 
 
  for ia in range(nat):
    atom = basis[ia][0]
    xyz=basis[ia][1]
    ofile.write("%12.6f %12.6f %12.6f \t  %s \n"%(xyz[0],xyz[1],xyz[2],atom))  

  ofile.close()
