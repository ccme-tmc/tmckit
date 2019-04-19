#!/usr/bin/env python
import sys,os,shutil
import commands,string
#from math import *
#from subprocess import *
#from scipy.optimize import *
#from eos_fitting import *
#from numpy import *
from list_utils import *
from constants  import *

def f_Read_Number(fn,dtype,irow,icol):
  """ 
  Read a number from a data file 
  """ 
  myname = "f_Read_Number"
  ifile = open(fn,'r')
  lines = ifile.readlines()
  # remove commenting lines
  for line in lines:
    if line.strip()[0]=='#': lines.remove(line)
 
  nrows = len(lines)
 
  if irow < 0: i = nrows + irow
  else: i = irow - 1
 
  if i >= nrows or i < 0:
    print myname+"==> WARNING: row index overflow "
    ret = None 
   
  line = lines[i].split()
  ncols = len(line)
 
  if icol< 0: j = ncols + icol
  else: j = icol - 1
 
  if j >= ncols or j < 0:
    print myname+"==> WARNING: column index overflow "
    ret = None 
   
  if dtype == 'f': ret = float(line[j])
  elif dtype == 'i': ret = int(line[j])
  else: ret = line[j]
  ifile.close() 
  return ret  

def f_Read_Data(data_file,format=None): 
  """
  Read multi-column data from a file 
  """ 

  ifile = open(data_file,'r') 
  data = []

  for line in ifile:
    if len(line.strip()) < 1 or line.strip()[0] == '#': continue # neglect comments   
    # extract data line by line 
    #  l_data -- the data extracted from the line 
    l_str = line.split()      # the line as a string
    l_n = len(l_str)         # the number of fields in each line  

    l_data=[]                # the data extracted from each line 

    if format is None:       # if format is None, it is assumed that all data are float 
      for i in range(l_n):
        s = l_str[i]
        if s[0] == '#': break # neglect comments   
        try: 
          d = float(s) 
        except:
          try:
            d = int(s) 
          except:
            d = s
        l_data.append(d) 
    else:
      ncol = len(format) 
      if ncol > l_n: # the number of entries in the "format" is bigger than that in the line 
        print "WARNING: the number of entries in format is bigger than that really present in the line!!"
        ncol = l_n 
      for i in range(ncol): 
        s = l_str[i] 
        fmt = format[i]
        if fmt == 'i':
          l_data.append(int(s))
        elif fmt == 'f':
          l_data.append(float(s))
        else:
          l_data.append(s)
      
    data.append(l_data) 
  ifile.close()
  return data 

def f_Col(data,ic):
  """ 
  Return ic-th colume of data as an array
  """
  col=[]
  nrow = len(data)
  for i in range(nrow):
    col.append(data[i][ic-1])
  return col

def f_Row(data,ir):
  """
  Return ic-th colume of data
  """
  return data[ir-1][:]

def f_Add_Col(data,col):

  n = len(col) 
  data_new = []
  if len(data) == 0:
    for i in range(n):
      data_new.append([ col[i] ])
  else: 
    for i in range(len(data)):
      r = data[i][:]
      r.append(col[i])
      data_new.append(r)

  return data_new

def f_Gauss(x,s):
  return exp(-x*x/(2*s*s))/sqrt(2*Pi*s*s)

def f_GaussConvl(x,y,w):

  h = x[1]-x[0]
  n = len(x)

  nw = 2*int(w/h)

  y_new = []
  for i in range(n): 
    if i < nw or i >= n - nw:
      yn = y[i]
    else: 
      yn = 0.0 
      for j in range(i-nw,i+nw+1):
        yn += y[j]*f_Gauss(x[i]-x[j],w)*h 
    y_new.append(yn) 
  return y_new 
    
  

  

     
