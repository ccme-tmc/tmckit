#!/usr/bin/env python
import sys,os,shutil
import commands,string

def f_Getopt(opflag,nval,def_val,debug=False):
  try:
    i_op = sys.argv.index(opflag)
    del sys.argv[i_op]

    if nval == 0: 
      val = True

    elif nval ==1: 
      if isinstance(def_val,int):
        val = int(sys.argv[i_op]); del sys.argv[i_op]
      elif isinstance(def_val,float):
        val = float(sys.argv[i_op]); del sys.argv[i_op]
      else:
        val = sys.argv[i_op].strip(); del sys.argv[i_op]
    else:

      val=[]
      for i in range(nval): 
        if isinstance(def_val[i],int):
          t = int(sys.argv[i_op]); del sys.argv[i_op]
        elif isinstance(def_val[i],float):
          t = float(sys.argv[i_op]); del sys.argv[i_op]
        else:
          t = sys.argv[i_op].strip(); del sys.argv[i_op]
        val.append(t)

  except: 
    val = def_val
    
  if opflag != '-h' and debug : print "option " + opflag + " = ",val 
  return val
    
def f_Get_Last_Line(file_name):
  try: 
    ifile=open(file_name,'r')
    lines=ifile.readlines()
    ifile.close()
    return lines[len(lines)-1]
  except:
    return ''

def f_Skip_Lines(ifile,n):
  i=0
  while i < n :
    ifile.readline()
    i += 1
  return

def f_Copy_Lines(ifile,ofile,n=None):

  if n is None:
    while 1:
      line=ifile.readline()
      if line:
        ofile.write(line)
      else:
        break 
  else:     
    for i in range(n):
      line=ifile.readline()
      ofile.write(line)
  return 

def f_Grep_Lines(file,tag,nl,nf):
  """
  Get the "nf"-th field (separated by space) on the "nl"-th line 
  that contains "tag"   
  """ 
 
  if nl > 0:
    grep_cmd = "grep --no-filename '" + tag + "' " + file + "| head -n " + "%4d"%(nl)
  else:
    grep_cmd = "grep --no-filename '" + tag + "' " + file + "| tail -n " + "%4d"%(-nl)

  print "grep_cmd='",grep_cmd,"'" 
 
  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    sys.exit(1)
  return output.split()[nf-1]

