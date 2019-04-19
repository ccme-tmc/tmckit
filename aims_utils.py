#!/usr/bin/env python

import commands,sys,string,math

# read the total energy from the output file of a scf file
def f_Read_Energy(aims_out):

  # check whehter it is converged SCF calculation
  cmd = "grep \"Self-consistency cycle converged.\" " + aims_out

  failure,output = commands.getstatusoutput( cmd )
  if failure:
    print "ERROR: Convergence fails! "
    sys.exit(1)

  # get the SCF total energy
  cmd = "grep \"| Total energy corrected        :\" " + aims_out + " | tail -n 1 "
  failure,output = commands.getstatusoutput( cmd )
  if failure:
    print "ERROR: Fail to find the total energy in SCF output"
    sys.exit(1)

  str = output.split()
  etot = float(str[5])
  return etot

# this function does a scf calculation and return the SCF total energy
def f_SCF_Energy(aims_x,aims_out):
  cmd = aims_x + ">>" + aims_out
  failure,output = commands.getstatusoutput( cmd )
  if failure:
    print "ERROR when running " + aims_x
    sys.exit(1)

  etot = f_Read_Energy(aims_out)
  return etot

def f_Write_Struct_aims(name,mol,latt_vec=None):
  """
  This subroutine write molecular structure into name-geometry.in in the FHI-aims geometry.in 
  """
  nat = len(mol)
  ofile = open(name+"-geometry.in",'w')
  print "Write molecular structure to FHI-aims format: " + name + "-geometry.in"
  for iat in range(nat):
    xyz  = mol[iat][1]
    atom = mol[iat][0]
    ofile.write("%-12s %12.6f %12.6f %12.6f %6s\n"%("atom",xyz[0],xyz[1],xyz[2],atom))
  if not latt_vec is None:
    for i in range(len(latt_vec)):
      vec=latt_vec[i]
      ofile.write("%s %12.6f %12.6f %12.6f\n"%("lattice_vector",vec[0],vec[1],vec[2]))
  ofile.close()

