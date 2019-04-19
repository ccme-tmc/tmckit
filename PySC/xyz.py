#!/usr/bin/env python
from math import *
from chem_utils import *
from io_utils import *

def f_xyz_read(name,debug=False):
  """
  This subroutine reads name.xyz
  """
  xyz_file = name.strip()+".xyz"

  ## check whether struct file exists
  if not os.path.isfile(xyz_file):
    print "ERROR in f_Read_XYZ: struct file " + xyz_file + " does not exist"
    sys.exit(1)

  ifile = open(xyz_file,'r')

  # read lattice types and the number of nonequvilanet atoms
  line = ifile.readline()
  nat = int(line.split()[0])
  if debug: print nat
  line = ifile.readline()
  if debug: print line.strip()
  mol=[]
  for iat in range(nat):
    line = ifile.readline().split()
    atom = line[0]
    xyz = [ float(line[1]),float(line[2]),float(line[3])]
    if debug: print "%-6s %12.6f %12.6f %12.6f"%(atom,xyz[0],xyz[1],xyz[2])
    mol.append( [atom,xyz] )

  ifile.close()
  return mol

def f_xyz_write(name,mol,title=None):
  """
  This subroutine write molecular structure into name.xyz in the xyz format
  """
  xyz_file = name.strip()+".xyz"

  print "Write molecular structure to xyz format: "+xyz_file

  ofile = open(xyz_file,'w')
  if title is None: title = "molecular structure written by f_xyz_write"

  # read lattice types and the number of nonequvilanet atoms
  nat = len(mol)
  ofile.write("%6d\n" %(nat))
  ofile.write("%s\n"%(title))

  for iat in range(nat):
    atom=mol[iat][0]
    xyz=mol[iat][1][:]
    ofile.write("%6s %12.6f %12.6f %12.6f\n"%(atom,xyz[0],xyz[1],xyz[2]))

  ofile.close()

