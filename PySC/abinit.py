#!/usr/bin/env python
import sys,os,shutil
from math import *
from chem_utils import *
from struct_utils import * 
mod_name = 'abinit.py'

def abi_write_struct(name,):
  sname = "abi_write_struct"

  ofile = open(name.strip()+".abi",'w')
  print "write the struct in the abinit-consistent format to "+name.strip()+".abi"

  # check the number of species
  ofile.write("# structure information generated by %s @ %s\n" %(sname,mod_name))
  ofile.write("\n# system info\n")
  ofile.write("SystemName %s\n"  %(name))
  ofile.write("SystemLabel %s\n" %(name))

  ofile.close()

def abi_read_struct(name,type='m'):
  """
  Read structure information from name. 
  """
  ifile = open(name.strip()+".STRUCT_OUT",'r')
 

