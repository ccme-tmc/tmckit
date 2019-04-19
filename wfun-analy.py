#!/usr/bin/env python

from data_utils import *
from constants import *

myname = sys.argv[0]

if len(sys.argv) < 2:
  print "ERROR: the file name for wfun data missing!\n!!!"
  print "  Usage: "+myname+" <wfun_data_file_name>\n!!!"
  sys.exit(1)
  
wfun = sys.argv[1]

fmt=['s','f','f','f','f','f']
data = f_Read_Data(wfun,format=fmt)
nr = len(data) 
nc = len(data[0])

# calculate work function from V_vac and evbm(slab) 
print "# wfun(1) == V_vac  - evbm_slab"
print "# wfun(2) == (V_vac - ecore_slab) -( evbm_bulk - ecore_bulk)" 
print "# default units: Ry"
print "#%-30s %12s %12s %12s %12s %12s %12s %12s"%("system","V_vac","evbm_slab","ecore_slab","evbm_bulk","ecore_bulk","wfun(1)/eV","wfun(2)/eV")
for i in range(nr):
  sys_name= data[i][0]
  V_vac   = data[i][1]
  ef_slab = data[i][2]
  ec_slab = data[i][3]
  ef_bulk = data[i][4]
  ec_bulk = data[i][5] 

  wf1 = (V_vac-ef_slab)*Ry2eV
  wf2 = ((V_vac-ec_slab)-(ef_bulk-ec_bulk))*Ry2eV
  
  print "%-30s %12.5f %12.5f %12.5f %12.5f %12.5f %12.3f %12.3f"%(sys_name,V_vac,ef_slab,ec_slab,ef_bulk,ec_bulk,wf1,wf2)

