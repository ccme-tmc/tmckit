#!/usr/bin/env python
from math import *
from chem_utils import *
from io_utils import * 
#
# This file contains subroutines that are used to read or write molecular structure in the MDL MOL format
#   f_Read_struct_mol

def f_mol_read(name,debug=False):
  """ 
  This subroutine reads the MDL MOL structure file name.mol  
  """
  mol_file = name.strip()+".mol"

  ## check whether struct file exists
  if not os.path.isfile(mol_file):
    print "ERROR in f_mol_read: struct file " + mol_file + " does not exist"
    sys.exit(1)

  ifile = open(mol_file,'r')

  # read lattice types and the number of nonequvilanet atoms
  f_Skip_Lines(ifile,3) 
  line = ifile.readline()
  nat = int(line[0:3])
  if debug: print "The Number of atoms: %d" %(nat)
  mol=[]
  for iat in range(nat):
    line = ifile.readline().split()
    atom = line[3]
    xyz = [ float(line[0]),float(line[1]),float(line[2])]
    mol.append( [atom,xyz] )

    if debug: print "%6s %12.6f %12.6f %12.6f # atom %6d"%(atom,xyz[0],xyz[1],xyz[2],iat+1)

  ifile.close()
  return mol


