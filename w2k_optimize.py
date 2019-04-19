#!/usr/bin/env python
import sys,os,shutil
import commands,string
from math import *
from numpy import *
from scipy.optimize import *

from eos_fitting import *
from list_utils  import *
from w2k_utils   import *
from optimize    import *
import gparams 

def f_w2k_Optimize(type = '',\
               method = '',\
               case='',\
               w2k_cmd='',\
               constraint = None, \
               xtol=1.e-2,  \
               ftol=1.e-6,  \
               fit_tol = 0.1, \
               fit_types = None, \
               step = 2.0,   \
               debug = True, \
               out_name='',\
               save_all = False,\
               use_init = False,\
               restart=False):
  """
  Optimize lattice parameters with percentage changes     
   <type> detertimes which lattice paramete to be optimized:
     "vol"         -- 1D, volume optimization with fixed a:b:c and lattice angles
     "coa"         -- 1D, optimize c/a with fixed volume (only for tetragonal and hexagonal lattices)
     "vol_coa"     -- 2D, volume and c/a (only for hexgonal and tetragonal lattices) 
     "hex"         -- 2D, a and c in hexagonal lattice 
     "tetra"       -- 2D, a and c in tetragonal lattice 
     "boa_coa"     -- 2D, optimize b/a and c/a with fixed volume (for orthombic, monoclinc and triclinic lattices) 
     "vol_boa_coa" -- 3D, volume, b/a and c/a (for orthombic, monoclinc and triclinic lattices) 
     "a_b_c"       -- 3D, optimize a,b,c (for orthombic, monoclinc and triclinic lattices)
     "vol_mono"    -- 4D, optimize volume, b/a, c/a and beta 
     "mono"        -- 4D, optimize a,b,c and beta
     "vol_tric"    -- 6D, optimize over vol, b/a, c/a, alpha, beta and gamma
     "tric"        -- 6D, optimize a,b,c, alpha, beta and gamma
  """
  case_name = f_Check_Name(case)
  struct_file = case_name+".struct"

  # if use_init is set true and case_init.struct exists, then it will be used as 
  # referece, otherwise the current case.struct will be used as the reference   

  struct_initfile = case_name+"_init.struct"

  if not use_init and not restart: 
    if os.path.isfile(struct_initfile): os.remove(struct_initfile)
    shutil.copy(struct_file,struct_initfile)

  latt_ini = f_Read_Latt(struct_initfile)
  latt_type = f_Latt_Type(latt_ini)
  print "Lattice type: " + latt_type
  print "Initial Lattice Constants:"
  print "     a=%12.6f,    b=%12.6f,     c=%12.6f" %(latt_ini[0],latt_ini[1],latt_ini[2])
  print " alpha=%12.6f, beta=%12.6f, gamma=%12.6f" %(latt_ini[3],latt_ini[4],latt_ini[5])

  # set and print w2k_cmd 
  if w2k_cmd == '': 
    w2k_cmd = 'run_lapw -ec 0.0000001'
  print "w2k_cmd=" + w2k_cmd 

  # if type is not set, then determine it in terms of lattice type
  if type == '':
    if   latt_type == 'cub'   : type = 'vol'         
    elif latt_type == 'tetra' : type = 'vol_coa'     
    elif latt_type == 'hex'   : type = 'vol_coa'     
    elif latt_type == 'ortho' : type = 'vol_boa_coa' 
    elif latt_type == 'mono'  : type = 'vol_mono'    
    elif latt_type == 'tric'  : type = 'vol_tric'    
    else : 
      print "ERROR: no appropriate optimization type available for lattice "+latt_type
      sys.exit(1)

  if   type == 'vol' or type == 'coa' or type =='c' or type == 'tmdc':        
    ndim = 1
  elif type == 'vol_coa' or type == 'a_c' or type == 'boa_coa':  
    ndim = 2
  elif type == 'vol_boa_coa' or type == 'a_b_c':                     
    ndim = 3
  elif type == 'vol_mono'    or type == 'mono':
    ndim = 4
  elif type == 'vol_tric'    or type == 'tric':   
    ndim = 6
  else:
    print "ERROR: unsupported type= "+ type
    sys.exit(1)

  p0 = f_List_New(0.0,ndim) 

# setup some output options 
#    out_file : store SCF total energy data at different lattice parameters
#    save_dir : used to save all intermediate data via w2k_save 
#    gparams.pyw2k_logfile: standard output from wien2k calculations 

  if out_name == '': 
    out_name = './opt-'+type.strip()

  out_file = out_name.strip()+".dat"

  save_flag = ''
  if save_all: 
    save_dir = out_name
  else: 
    save_dir = ''

  if restart:
    line=f_Get_Last_Line(out_file)
    print "Last line of out_file", line
    gparams.save_counter = int(line.split()[0]) 
  else:
    gparams.save_counter = 0 
    if os.path.isfile(gparams.pyw2k_logfile): os.remove(gparams.pyw2k_logfile)
    if os.path.isfile(out_file): os.remove(out_file)
    
  print "save_counter start with",gparams.save_counter
    
  args = (latt_ini,type,case_name,w2k_cmd,out_file,save_dir,save_flag)
  func = f_Energy_vs_Latt_Change

  if fit_types is None: fit_types = f_List_New("H",ndim)

  if method == '': method = "BracketFit"

  print "Optimize " + type + " by " + method + " method"

  if method == "BracketFit":
    delta_x = f_List_New(step,ndim)
    out = f_Minimize_BracketFit(func,p0,\
                 constraint=constraint, dxs=delta_x, args=args, xtol=xtol, ftol=ftol, fit_tol=fit_tol,\
                 fit_types=fit_types,itmax=10,debug=debug,restart=restart)

  elif method == "TaxiCab":
    out = f_Minimize_TaxiCab(func,p0,args=args,constraint=constraint,xtol=xtol,ftol=ftol,\
                             itmax=200,debug=debug,restart=restart)

  elif method == "Powell":
    out = f_Minimize_Powell(func,p0,args,xtol=xtol,ftol=ftol,itmax=200,debug=debug,restart=restart)

  elif method == "ScanFit":
    scan_ranges = f_List_New( (-step*2,step*2,5), ndim )
    out = f_Minimize_ScanFit(func,p0,scan_ranges,args=args,constraint=constraint,xtol=xtol,ftol=ftol,\
                             fit_tol=fit_tol,fit_typs=fit_types,itmax=10,debug=debug,restart=restart)

  elif method == "fmin_powell":
    out = fmin_powell(func,p0,args=args,full_output=True,xtol=xtol,ftol=ftol)
    nfcall = out[5]
    print "Total number of function calls:%d" %(nfcall) 

  elif method == "fmin":
    out = fmin(func,p0,args=args,full_output=True,xtol=xtol,ftol=ftol)
    nfcall = out[5]
    print "Total number of function calls:%d" %(nfcall)
  
  else:
    print "ERROR in f_Opt_nD: unsupported method " + method
    sys.exit(1)

  pmin = out[0]
  emin = out[1]
  print "Optimized paramters (in percentage changes with respect to the initial structure):"
  fmt = f_Str_New('%12.6f ',ndim)
  print fmt %(tuple(pmin)) + "%16.6f"%(emin)
  return emin

