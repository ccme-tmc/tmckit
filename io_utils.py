#!/usr/bin/env python
import sys,os,shutil
import commands,string,subprocess

def f_Getopt(opflag,nval,def_val,debug=False):
  try:
    i_op = sys.argv.index(opflag)
    del sys.argv[i_op]
    print i_op

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
    
  if opflag != '-h' and debug :
    print "option " + opflag + " = ",val 
  return val

def io_get_val(file_in,tag,dtype,i=0,mode=1):
  """ 
  get the value of a paramter from a free-style text file assuming it is 
  represented in the following format
   "tag = val"
  or 
   "tag val"
  or 
   "tag : val" 
  or the value is in the i-th field ( seperated by a blank space) in a line containing <tag>

  dtype -- the type of the value, f/i/s (float/int/string) 
  mode  -- control the mode of reading 
      n!=0 return the n-th value found in the file 
      0    read all values in match  
  Any lines starting with "#" or "!" will be neglected 
  """
  ifile= open(file_in,'r')
  nw = len(tag)

  val=[]

  while (1):
    line = ifile.readline()
    if not line: break    # end of the file reached 
    
    line = line.strip()
    if len(line)==0 or line[0]=='#' or line[0] == '!':  continue 

    if tag in line:

      if i != 0: 
        ltmp = line.split() 
        if i>0:
          val_s = ltmp[i-1]
        else:
          val_s = ltmp[i]
      else:
        ltmp = line[nw:].strip()
        if ltmp[0]=='=' or ltmp[0]==':': 
            ltmp = ltmp[1:]

        ltmp_s = ltmp.split()
        val_s = ltmp_s[0]  

      if dtype =='f':
        val.append(float(val_s))
      elif dtype =='i':
        val.append(int(val_s))
      else:
        val.append(val_s)
  ifile.close()

  if len(val)==0: 
    return None
  else:
    if mode == 0:
      return val
    elif mode > 0: 
      return val[mode-1]
    else:
      return val[mode]

def io_read_stdout(cmd):
  '''
  Open a Shell with specific command and return its stdout 
  '''
  ProcSub = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell=True)
  return ProcSub.stdout.read()


def io_set_val(fl,v_tag,v_val,pos=None, eq_symbl='='):
  """
  set or reset the input parameter v_tag in fl to the value of "v_val"
    v_type == 's'  a string 
              'i'  an integer 
              'f'  a float number
    if pos="ir ic", is set, the parameter in the block (in the liboctparser format) is assumed  
  """
  # make a copy of original INCAR file 
  fl_bak = fl + "_bak"
  if os.path.isfile(fl_bak): os.remove(fl_bak)
  os.rename(fl,fl_bak)

  ifile = open(fl_bak,'r')
  ofile = open(fl,'w')

  done_reset = False

  while 1: 
    line = ifile.readline() 
    if not line: break 

    tmp = line.strip()
    # copy the comment lines
    if len(tmp) == 0 or tmp[0] == '#' or tmp[0] == '!':
      ofile.write(line)
      continue 

    line_out = ''
    if v_tag in line and (not pos is None):  # reset a block parameter 

      # extract the location in the block 
      tmp = pos.split()
      ir = int(tmp[0])
      ic = int(tmp[1])

      ofile.write(line)   # write the block name 

      # read the block 
      il = 0 
      while 1: 
        line = ifile.readline() 

        if line.strip()[0] == '%': break # reach the end of the block

        il += 1 
        if ir != il: 
          ofile.write(line)
        else:
          tmp    = line.split('#') 
          params = tmp[0].split('|') 
          if len(tmp)>1: 
            info = tmp[1].strip() 
          else:
            info = ''
  
          info = info + ' Note: the %d-th paramerer is reset'%(ic) 

          params[ic-1] = v_val 
 
          nc = len(params)
          for i in range(nc):
            ofile.write(params[i])
            if i < nc-1 : ofile.write('  | ') 

          ofile.write("\t\t # " + info + "\n") 

      ofile.write('%\n') 
      done_reset = True

    elif v_tag in line:        
      # first check whether the line contains comments 
      tmp = line.split('#')
      if len(tmp) > 1:
        line1 = tmp[0]
        comments = tmp[1]
      else:
        line1 = line
        comments = ''

      # check whether the line contains multiple input parameters 
      line_s = line1.split(';')
      line_out = ''
      nv = len(line_s)
      for i in range(nv):
        name_val = line_s[i].split(eq_symbl)

        if name_val[0].strip().lower() == v_tag.lower() : 
          line_out = line_out + name_val[0] + ' ' + eq_symbl + ' ' + v_val 
          done_reset = True
        else:
          if i>0: 
            line_out = line_out + " ; " + line_s[i] 
          else:
            line_out = line_s[i]

      line_out = line_out + ' # ' + comments + '\n'
      ofile.write(line_out)

    else:
      ofile.write(line)

  if not done_reset: 
    if pos == None: 
      ofile.write(v_tag + ' ' + eq_symbl + ' ' +  v_val + '\n')
    else:
      print "ERROR: the block %s is not found!!!"%(v_tag) 

  ofile.close()
  ifile.close() 

def io_count(fl,flag):
  """
  Count the number of appearances of "flag" in the fl
  """
  grep_cmd = grep_cmd = "grep -c '" + flag + "' " + fl
  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    return None
  else:
    return int(output)
    
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

def f_Read_Lines(file,il_start,il_end=0):
  ifile = open(file,'r')

  il = 0
  lines_all = ifile.readlines()
  ifile.close()
  nlines = len(lines_all) 

  if il_start > 0: 
    il0 = il_start-1
  else:
    il0 = il_start 
 
  if il_end == 0:
    il1 = il0+1
  elif il_end < il_start: 
    print "WARNING: il_end < il_start!" 
    print " -- read all lines from il_start"
    if il_start> 0: 
      il1 = nlines 
    else:
      il1 = -1
  else:
    il1 = il_end 

  return lines_all[il0:il1] 

def io_skip_lines(ifile,n):
  i=0
  while i < n :
    ifile.readline()
    i += 1
  return


def io_read_lines_tagged(file,tag_begin,iop=0,tag_stop=None,mode=0):
  """
  Read lines from <file> between <tag_begin> and <tag_stop>  
   if iop  = 0, then it stops when a blank line is reached 
   if iop  < 0, it stops when <tag_stop> appears 
   if iop > 0, it stops after reading <iop> lines 
  mode: control which data that meets the requirements are returned 
    0 return all data 
    1 return the first set
    2 return the last ones 
  """
  ifile = open(file,'r') 

  lreach_begin = False
  lread_end = False 
  
  lines_all=[]
  while (1):
    line = ifile.readline()
    if not line: break       # reach the end of file

    if tag_begin in line: lreach_begin = True

    if lreach_begin:
      lines=[]
      if iop <= 0: 
        while(1): 
          line = ifile.readline()
          if iop == 0: 
            if line.strip() == '': break 
          else: 
            if tag_stop in line: break 

          lines.append(line)
      else:
        for i in range(iop):
          line = ifile.readline()
          lines.append(line)
           
      lreach_begin = False  
      lines_all.append(lines) 
  
  ifile.close() 
  if mode == 0: 
    return lines_all      
  elif mode == 1:
    return lines_all[0]
  elif mode == -1:
    return lines_all[-1]

def io_copy_lines(ifile,ofile,n=None):

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
  return line 

def io_get_line(file,nl):

  ifile=open(file,'r')
  io_skip_lines(ifile,nl-1)
  line = ifile.readline()
  ifile.close()

  return line 

def io_grep_lines(file,tag,nl,nf,debug=False):
  """
  Get the "nf"-th field (separated by space) on the "nl"-th line 
  that contains "tag". If nl=0, then all data on  "nf"-th field (separated by space) in the lines 
  that contains "tag"
  """ 

  if nl == 0:
    grep_cmd = "grep -a --no-filename '" + tag + "' " + file 
  elif nl > 0:
    grep_cmd = "grep -a --no-filename '" + tag + "' " + file + "| head -n " + "%4d"%(nl)
  else:
    grep_cmd = "grep -a --no-filename '" + tag + "' " + file + "| tail -n " + "%4d"%(-nl)
 
  failure,output = commands.getstatusoutput( grep_cmd )
  if failure:
    print "ERROR when running " + grep_cmd
    sys.exit(1)

  if debug: print output

  if output.strip() == '':
    val = None 
  else: 
    if nl == 0:
      val=[]
      tmp = output.split("\n") 
      for i in range(len(tmp)):
        ltmp = tmp[i].split()
        if nf > 0:
          val.append(ltmp[nf-1])
        else:
          val.append(ltmp[nf]) 
 
    else:
      if nf > 0: 
        val = output.split()[nf-1]
      else:
        val = output.split()[nf]

  return val 

def io_dat2dict(file_name): 
    """
    this function reads the data in a file to a dictionary 
    """


