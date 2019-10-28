#!/usr/bin/env python

from __future__ import print_function
import sys
import os
import copy
import subprocess as sp
import math
from math import sqrt, cos, sin
import re
import itertools
try:
    from commands import getstatusoutput
except ImportError:
    def getstatusoutput(cmds, shell=True):
        """Emulator to the old commands.getstatusoutput function
        """
        try:
            output = sp.check_output(cmds, shell=shell)
            return 0, output
        except sp.CalledProcessError as err:
            return err.returncode, err.output

#import numpy,numpy.linalg
import constants
import list_utils as lu
import warnings
import functools
import kpt_spec
from py_xmlio import XMLSerilizer

logLevel = 0

def deprecated(func):
    """This is a decorator which can be used to mark functions
    as deprecated. It will result in a warning being emitted
    when the function is used."""
    
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        warnings.warn_explicit(
            "Call to deprecated function %(funcname)s." % {
                'funcname': func.__name__,
            },
            category=DeprecationWarning,
            filename=func.func_code.co_filename,
            lineno=func.func_code.co_firstlineno + 1
        )
        return func(*args, **kwargs)
    return new_func

  
def debug_print(level, *para):
    '''
    print for debug, logLevel can be used to control whether output
    '''
    if ( level <= logLevel ):
        for i in para:
            print(i,end=" ")
        print('')

class FormatError(Exception):
    '''
    Represent file format incorrect
    '''
    def __init__(self,value):
        self.value = value
    
    def __str__(self):
        return self.value

def f_ReadFileFloatValue(stFileName):
    '''
    Read a single floating point value of one file
    If the file does not exists, then return None
    '''
    d = None
    if ( os.path.exists(stFileName)):
        fIn = open(stFileName,'r')
        d =  float(fIn.readline().strip())
        fIn.close()
    return d  

def f_ReadFileIntValue(stFileName):
    '''
    Read a single integer point value of one file
    If the file does not exists, then return None
    '''
    d = None
    if ( os.path.exists(stFileName)):
        fIn = open(stFileName,'r')
        d =  int(fIn.readline().strip())
        fIn.close()
    return d  
    
def f_WriteFileFloatValue(val,stFileName):
    '''
    Write a single floating point value to one file
    '''
    f = open(stFileName,'w')
    f.write("%f" % val)
    f.close()

def f_WriteFileIntValue(val,stFileName):
    '''
    Write a single integer value to one file
    '''
    f = open(stFileName,'w')
    f.write("%i" % val)
    f.close()

def f_write_file_string(val,filename):
    '''
    Write a single string to one file
    '''
    f = open(filename,'w')
    f.write(val)
    f.close()

def f_FindGCD(n1,n2):
    '''
    Find the greatest common divisor of two integers
    '''
    if ( n2 == 0):
        return n1
    if ( n1 > n2 ):
        return f_FindGCD(n2,n1 % n2)
    else:
        return f_FindGCD(n1,n2 % n1)

def f_FindLCM(n1,n2):
    '''
    Find the least common multiple of two integers
    '''
    return n1*n2/f_FindGCD(n1,n2)


def f_FindFraction(dIn,dErrMax=1.0e-6,dMaxDiv=10):
    '''
    Find the most appraoching fraction for given float number

    :param dIn: the floating number
    :param dErrMax: the relative maximum error allowed between fraction and dIn
    :param dMaxDiv: the maximum divider allowed
    :return: two integer to represent the fraction. If not found, return two None
    '''
    for i in xrange(1,dMaxDiv+1):
        d2=  dIn*i
        if ( abs(d2-round(d2)) < dErrMax):
            return int(round(dIn*i)),i

    return None,None

def f_CorrectFloat(dIn,dErrMax=1.0e-6):
    '''
    make a float number reach the most exact value ( like 60.00000001 to 60 )
    also make it reach 1/6, 1/5, 1/4 , 1/3 1/2 and so on

    :param dIn: The number to be detect and modify, or an array ( every number in it will be modified)
    :param dErrMax: The error can be treated as floating point error
    :return: 0 if not found, otherwise x in 1/x 
    '''
    def Correct(dOld,dErrMax):
        d2 = dOld
        bOK = False
        nMax10 = 2 #  integer from 1-10^nMax10 as divide
        for j2 in range(0,nMax10):
            for j1 in range (1,10):
                j = j1 * 10**j2
                if ( abs(dOld*j - round(dOld*j,0) ) < dErrMax*j):
                    #print("Correct with %7i : %20.14f to %20.14f"%(j,result,round(dOld*j)/j))
                    d2 = round(dOld*j)/j
                    bOK = True
                    break
            if ( bOK ):
                break
        return d2


    result = copy.deepcopy(dIn)
    if ( isinstance(result,list) ):
        for k in range(0,len(result)):
            result[k] = Correct(dIn[k],dErrMax)
           #d1 = dIn[j]
           #for i1 in range (1,10):
           #    for i2 in range (0,7):
           #        i = i1 * 10**i2
           #        if ( abs(d1*i - round(d1*i,0) ) < dErrMax*i):
           #            result[j] = round(d1*i)/i
           #            break
    elif ( isinstance(result,float)):
        result = Correct(dIn,dErrMax)
       #for i1 in range (1,10):
       #    for i2 in range(0,7):
       #        i = i1 * 10**i2
       #        if ( abs(dIn*i - round(dIn*i,0) ) < dErrMax*i):
       #            #print("Correct %20.14f to %20.14f"%(result,round(dIn*i)/i))
       #            result = round(dIn*i)/i
       #            break
    else:
        raise TypeError('Unsupported type to correct')
        
    return result
            
def f_File_IsNullOrEmpty(st):
    '''
    Determine whether a file does not exist, or is empty
    '''
    if ( os.path.exists(st) ):
        if ( os.stat(st).st_size != 0 ):
            return False 
    return True

def f_GetLibDataPath():
    '''
    Get PYTHONPATH for read some data file in the path
    '''
    return sys.path[0]


def f_GetExecFullPath(stExe):
    '''
    Get a program full path
    '''
    status,stOut= commands.getstatusoutput("which "+stExe)
    stHome = os.getenv("HOME")
    #replace home folder
    if ( stOut[0] == "~"):
        stOut = stHome+stOut[1:]
    return stOut

def f_File_GetLibraryPath():
    '''
    Get the path of python library file itself
    '''
    return os.path.dirname(__file__)

def f_file_ensuredir(dirname):
    '''
    Test whether a directory exists. If not then create it.

    :param dirname: the directory to create. If it is empty or None, then set to current directory
    :return: the diretory path guaranteed to work
    '''
    if (dirname is None or dirname.strip() == ""):
        dirname = os.getcwd()
    elif (not os.path.exists(dirname)):
        os.mkdir(dirname)
    return dirname

def f_file_ensure_no_file(filename):
    '''
    Test whether a file  exists. If yes then delete it.
    '''
    if (os.path.exists(filename)):
        os.remove(filename)

def f_ReadStdout(cmd):
    '''
    Open a Shell with specific command and return its stdout 
    '''
    ProcSub = sp.Popen(cmd,stdout=sp.PIPE,shell=True)
    return ProcSub.stdout.read()

def f_env_GetJobSystem():
    '''
    Get Job Management System Name
    '''
    result = None
    if ( os.environ.has_key("LSF_VERSION") ):
        result = "LSF"
    elif ( "pbs" in commands.getoutput("which qsub")  ):
        result = "PBS"
        
    return result

def f_env_MpirunCommand(stCommand,nProcessIn = 0,stJobSystem = None):
    '''
    create a string to run command with mpirun
    if enviroment varialbe tmckit_serial = 1, then mpirun will not be used ( depreceted )
    if nProcess = 0 and enviroment varialbe tmckit_process is set, then use $tmckit_process as nProcess will be used

    :param stCommand: The main command used after mpirun
    :param nProcess: The count of process. If not used or set as 0, default to hostfile
    :param stJobSystem: Specify Job Management System. If not used, default as auto detect. If no job system detected, the command will be used directly.
    '''
    nProcess = nProcessIn
    if (nProcess == 0 and os.environ.has_key("tmckit_process")):
        nProcess = int(os.environ["tmckit_process"])
         
    if ( nProcess != 0 ):
        if ( nProcess != 1):
            stRun = "mpirun -np %d %s" % (nProcess,stCommand)
        else:
            stRun = stCommand
    else:        
        dicJob = {"PBS":"PBS_NODEFILE","LSF":"LSF_MCPU_HOSTS"}
        if ( stJobSystem == None):
            stJobSystem = f_env_GetJobSystem()
        if ( stJobSystem != None):
            stRun = "mpirun -hostfile %s %s"%(os.environ[dicJob[stJobSystem]],stCommand)
        else:
            stRun = stCommand
            
    print("\33[34mCommand:\33[m%s" %stRun)
    return stRun   

def f_env_RunMpirunCommand(stCommand,nProcess = 0,stJobSystem = None):
    '''
    run command with mpirun
    '''
    status,stStdOut = commands.getstatusoutput(f_env_MpirunCommand(stCommand,nProcess,stJobSystem))
    return status,stStdOut    

def f_GetCaseName():
    '''
    Get current folder name, used in Wien2K as casename
    '''
    return os.path.basename(os.getcwd())




def GetDefaultCharge(nIndex):
    '''
    Give a default Charge for specific atom index
    '''
    result = -1
        
    nT = 0
    listAdd = [2,8,8,18,18,32,32]
    for nRow,nAdd in enumerate(listAdd):
        nT = nT+nAdd
        if ( nT >= nIndex):
            nLeft = nIndex - nT +nAdd
            nRight = nT-nIndex
            if ( nT == nIndex):
                result = 0 # rare gas
            elif ( nLeft <= 3 ):
                result = nLeft # IA-IIIA,IIIB
            elif ( nRight <= 7):
                if ( nRight+nRow <=5 ): #Diagnoal split
                    result = -nRight # - IIIA-VIIA
                else:
                    result = 8 -nRight # + ,IB-IIB IIIA-VIIA 
            elif ( nRight <= 10 ): # IIIV
                result = 2
            elif ( nRight <= 14 ): #IVB-VIIB
                result = 18-nRight
            else: # La/Ac series
                result = 3            
            break
    return result
        

def ReadSplitLine(fIn,stSep=""):
    '''
    Read a line from file and split it by default ( no EOL and ignore empty), or split with specific string
    '''
    if (stSep==""):
        return fIn.readline()[:-1].split()
    else:
        return fIn.readline()[:-1].split("\t")


def f_IsNumber(stInput):
    '''
    Detect if this string can be converted to a float number, return true or false
    '''
    result = True
    try:
        f = float(stInput)
    except:
        result = False
    return result
    
def f_Data_ReadTwoCol(stFileName,stComment="#",func=None):
    '''
    Read data from empty line seperated two-columns data, which can be read directly by xmgrace 
    This function require all sets in the file must have same number of rows

    :param stComment: If a line start with this character, then treat it as comment 
    :return: a list of rows
    '''
    f = open(stFileName)
    data = []
    bStarted = False
    nRow = 0
    nCol = 0
    for stLine in f.readlines():
        stLine = stLine.strip()
        if ( len(stLine) ==0):#Empty line, which means a set end
            if ( bStarted):
                if ( nRow != len(data) ):
                    raise ValueError("The %i set has different number of rows %i from previous %i" % (nCol,nRow,len(data)))
                else:
                    nRow = 0
                    nCol += 1
                    bStarted = False
                    continue
            else:#Just skip 
                continue
            
        if ( stLine[0] == "#"): #Skip comments
            continue

        #ar = [ float(x) for x in stLine.split()]
        ar = stLine.split()
        if ( len(ar) != 2):
            raise ValueError("There should be only 2 value in one line")
#Find data, mark as started        
        bStarted = True
        if (func is not None):
          if ( nCol == 0 ):
              data.append([func(x) for x in ar])
          else:
              data[nRow].append(func(ar[1]))
        else:
          if ( nCol == 0 ):
              data.append(ar)
          else:
              data[nRow].append(ar[1])
        nRow += 1
#First set

    f.close()
    return data

def f_Data_WriteTwoCol(data,stFileName):
    '''
    Write data to empty line seperated two-columns data, which can be read directly by xmgrace 
    '''
    f  = open(stFileName,'w')
    nLen = len(data[0])
    for i,row in enumerate(data):
        if ( len(row) != nLen):
            raise ValueError("Row %i has different columns from set" % i)
    for i in xrange(1,len(data[0])):
        for j in xrange(len(data)):
            f.write("%14s    %14s\n" % ( str(data[j][0]) , str(data[j][i])))
        f.write("\n")

    f.close()

    return 0

def f_Data_ReadMultiCol(stFileName,stComment="#",func=None):
    '''
    Read data from multi-columns data

    :param stComment: If a line start with this character, then treat it as comment 
    :param func: the function to convert each element string to specific object
    '''
    f = open(stFileName)
    nCol = -1
    result = []
    i = -1
    for stLine in f.readlines():
        i += 1
        stLine = stLine.strip()
        if ( len(stLine) == 0):
            continue
        if ( stLine[0] == "#"):
            continue

        ar_stLine = stLine.split()
        if ( nCol == -1):
            nCol = len(ar_stLine)
        elif ( nCol != len(ar_stLine)):
            raise ValueError("Unmatched column number in row #%i" % i)
        if (func is not None):
            ar_stLine = [func(x) for x in ar_stLine]
        result.append(ar_stLine)
    f.close()

    return result

def f_Data_WriteMultiCol(data,stFileName,delim="  ",func=None):
    '''
    Write data into multi-columns with space seperated 

    :param data: 2-D list contains the data
    :param stFileName: the output filename
    :param func: the function to convert element to string. If not specified, set to str()
    '''
    f = open(stFileName,"w")

    if (func is None):
        func = str

    for line in data:
        f.write(delim.join([func(x) for x in line]))
        f.write("\n")

    f.close()
    return 0

def f_Data_Write(list_data,filename,b_col=False,b_twocol=False,b_complex=False,b_split=True):
    '''
    This function can write data into files in the two-columns or the multi-columns format.
    The data is a nested lists, can be row first or column first.

    :param list_data: list of [x,y1,y2,y3...]
    :param b_col: the data is in the format of [[x,x,x...],[y1,y1,y1]...
    :param b_twocol: write the data in two-columns format instead of [x,y1,y2..] format
    :param b_split: write all real parts of y first and then imaginary parts: x,y1.real,y2.real,...y1.imag,y2.imag...
    :param b_complex: can be used when the data is complex, the real and imaginary parts will be in different columns
    '''
    if (b_col):
        list_data2 = list(zip(*list_data))
        f_Data_Write(list_data2,filename,b_col=False,b_twocol=b_twocol,b_complex=b_complex,b_split=b_split)
        return

    if (b_complex):
        if (b_split):
            list_data2 = [[x[0]]+[y.real for y in x[1:]]+[y.imag for y in x[1:]] for x in list_data]
        else:
            list_data2 = [[x[0]]+sum([[y.real,y.imag] for y in list_data],[]) for x in list_data]
    else:
        list_data2 = list_data

    if (b_twocol):
        f_Data_WriteTwoCol(list_data2,filename)
    else:
        f_Data_WriteMultiCol(list_data2,filename)

    return  

class NamelistIO(object):
    '''
    A Class used for Fortran namelist like file I/O
    Warning: This invokes exec on lines directly, please use with caution!
    '''
    def __init__(self,pre_line=None):
        '''
        :param pre_line: Preprocess all input line, accept a string (stripped) and returns a string
        '''
        self.dicExistName = {} #! the dictionary include &-/ part
        self.listExistName = [] #! the list contain the order of keys of dicExistName
        self.list_stLine = [] #! the file content
        self.stCurrentPart = "" #! temp
        self.pre_line = pre_line 
        self.b_warn_write = True #! Display warning when overwriting some files 

    @staticmethod
    def default_pre_line(line):
        '''
        Treat all "!" as the starting point of a comment
        Remove all redundant blanks and comments
        This is used when pre_line is None
        @todo we don't check if it is inside a string
        '''
        ix = line.find("!")
        if (ix != -1):
            line = line[:ix]
        return line.strip()

    
    def __DetectFileExist__(self,stFileName):
        '''
        Detect whether a file is exists, display warning if it is
        '''
        if ( os.path.exists(stFileName)):
            if (self.b_warn_write):
                print("\33[36mWarning: %s already exists! It is overwritten.\33[m" % stFileName)
            return True
        else:
            return False

    def __FormatValue__(self,stValue):
        '''
        Detect whether a value is a float value, process 'd' or 'D' to 'e'
        '''
        if ( '"' in stValue or "'" in stValue): #string
            return stValue
        stV2 = stValue.lower()
        i = stV2.find(".true.")
        if ( i != -1):
            stValue = stV2[:i] + "True" + stV2[i+6:]
        i = stV2.find(".false.")
        if ( i != -1):
            stValue = stV2[:i] + "False" + stV2[i+7:]
        if ( (stValue.find("d") != -1 and stValue[-1] != "d") or ( stValue.find("D") != -1 and stValue[-1] != "D" ) ): # float
            #ensure float: d must between in 2 digit or '-' and no other digit is allowed
            if ( not stValue.lower().replace("d","").replace("-","").replace(".","").replace(" ","").isdigit() ):
                return stValue
            else:
                stNext = stValue[stValue.lower().find("d")+1]
                if ( stNext != "-" and (not stNext.isdigit())):
                    return stValue
                else:
                    #print ( "Detect float : %s" % stValue)
                    return  stValue.lower().replace("d","e")
        return stValue
 
    def __GetNameValue__(self,stLine):
        '''
        Read a line contain "a=b" 
        Return a
        '''
        arLine = stLine.strip().split("=")
        stVarName = arLine[0].strip()
        if ( stLine.find("=") == -1):
#treat as empty line
#           print("Namelist Input Format Error!")
            return None,None

        stVar = self.__FormatValue__(arLine[-1])

        if ( stLine.find("(")!= -1):
            stVarName = stLine.strip().split("(")[0]
            if ( not hasattr(self,stVarName)):# new 100 length array
                setattr(self,stVarName,[-1 for x in range(0,100)])
            #set value
            #split
            reIndex = re.compile(pattern="(.+)\((.+)\)(.*)")
            aMatch = reIndex.match(arLine[0])
            if ( aMatch != None):
                stVarNameTotal = aMatch.group(1)+"["+str(int(aMatch.group(2)))+"]"
            else:
                print("Quamtum-Espresso Input Format Error!")
                raise FormatError("Quamtum-Espresso Input Format Error!")
        else:
            stVarNameTotal = stVarName

        try:
            exec("self."+stVarNameTotal.strip()+"="+stVar.strip())
        except NameError:#Treated as a string
            setattr(self,stVarNameTotal.strip(),stVar.strip())

        return stVarName, stVar
    
    def __ToStringFromName__(self,stName,nMaxIndex = 100) :
        '''
        Create a Name = Value string
        '''
        if ( not hasattr(self,stName) ):
            print("Error: Invalid property %s" % stName)
            return ""

        result = "  "+ stName + " = "
        value = getattr(self,stName)
        return self.__ToStringFromNameValue__(stName,value,nMaxIndex)


    def __ToStringFromNameValue__(self,stName,value,nMaxIndex = 100):
        '''
        Create a Name = Value string, value comes from self.name
        '''
        result = "  "+ stName + " = "
        if ( isinstance(value,bool)):
            if ( value):
                result += ".true."
            else:
                result += ".false."        
        elif ( isinstance(value,int)):
            result += str(value)
        elif ( isinstance(value,float)):
            if(value <= 0.001 and value >= -0.001): # common converge standard like 1d-10
                result += ("%2.2e" % value).replace("e","d")
            else:
                result += str(value) 
        elif ( isinstance(value,str)):
            result += "'"+value+"'"
        elif ( isinstance(value,list)): # array list, if < 0 then treat as not output
            result = ""
            for i in range(1,len(value)): #ignore value[0] 
                if( value[i] > -1): #ignore -1 float type ( in pw.x celldm )
                    result += self.__ToStringFromNameValue__("%s(%d)" % ( stName,i),value[i],nMaxIndex)
        else:
            print("Unsupported type : %s = %s" % (stName,str(value)))
        #print("line:" + result)
        if ( result[-1] != "\n"):
            result += "\n"
        return result
    
    def __WriteNameValue__(self,fOut,stPart,list_stName):
        '''
        Write a part 
        including:
        part name,  (if part name is not empty)
        all Name = value string in list_stName,
        and all name readed
        '''
        if ( stPart != ""):
            fOut.write("&%s\n" % stPart)
        #write readed
        if ( self.dicExistName.has_key(stPart)):
            listNameRead = self.dicExistName[stPart]
            #print("find %s" % stPart)
            #print(listNameRead)
        else:
            listNameRead = []
            #print("Not find %s" % stPart)
        
        for stName in listNameRead:
            #print("In Readed: " + stName)
            fOut.write(self.__ToStringFromName__(stName))
        
        for stName in list_stName:
            if ( not stName in listNameRead):
                #print("In Specific: " + stName)
                fOut.write(self.__ToStringFromName__(stName))
            else:
                pass
                #print("Duplicate: " + stName)            

        if (stPart != ""):
            fOut.write("/\n\n")

    def __AddPart__(self,name_part):
        '''
        Add a part in the index
        '''
        #print("Add Part %s" % stLine[1:].strip())
        self.stCurrentPart = name_part 
        self.dicExistName[self.stCurrentPart]= []
        self.listExistName.append(self.stCurrentPart)
        return

    def __ProcessLine__(self,i0,stLineOrg):
        '''
        Process a line in standard way. The two parameters is for convienience to edit.

        :param i0: the line number
        :param stLineOrg: the line content
        '''
        i = i0
        stLine = NamelistIO.default_pre_line(stLineOrg)
        # normal part
        while(True):
            if ( len(stLine) == 0): # skip empty line
                i += 1
                break
            if ( stLine[0] == '&' or stLine[0] == "/" ): #mark line, skip
                if ( stLine[0] == '&'):
                    self.__AddPart__(stLine[1:].strip().lower())
                i += 1
                break
            # comma split
            if ( stLine.find(",") != -1): # split ','
                del self.list_stLine[i]
                for stSplit in stLine.split(","):
                    self.list_stLine.insert(i,stSplit)
                break
            stCurrentName,var_tmp = self.__GetNameValue__(stLine)
            if (stCurrentName is not None):
                if (self.stCurrentPart == "" and not self.dicExistName.has_key("")):
                    self.__AddPart__("")
            #print("Current Name %s in part %s" % (stCurrentName,stCurrentPart))
                if ( not stCurrentName in  self.dicExistName[self.stCurrentPart]):
                    self.dicExistName[self.stCurrentPart].append(stCurrentName)
            i += 1
            break
        return i
    
    def AddValue(self,part,name,value):
        '''
        Add a name=value contained in specific &part/
        '''
        if ( not part in self.dicExistName.keys()):
            self.dicExistName[part] = []
            self.listExistName.append(part)
        if ( not name in self.dicExistName[part]):
            self.dicExistName[part].append(name)
        setattr(self,name,value)
    
    def WriteToFile(self,stFileName):
        '''
        Write content to one file in Quantum-Espresso readable format
        '''        
        self.__DetectFileExist__(stFileName)
        fOut = open(stFileName,'w')
        for stPart in self.listExistName:
            self.__WriteNameValue__(fOut,stPart,[])
    
    def ReadFromFile(self,stFileName):
        '''
        Read content from one file in Fortran namelist format
        '''              
        fIn = open(stFileName,'r')
        self.list_stLine = fIn.readlines()
        i = 0
        self.stCurrentPart = ""
        while ( i < len(self.list_stLine)):
            stLine = self.list_stLine[i].strip()
            if (hasattr(self,"pre_line") and self.pre_line is not None):
                stLine = self.pre_line(stLine)
            i = self.__ProcessLine__(i, stLine)
        #clean
        self.list_stLine = []
        pass    

class RawFileIO():
    '''
    A class represents a files which can be used for IO
    It keeps exactly every line in the file, and won't lost anything between I/O, but still support value edit and addition ("\\n" is stored)
    It is suitable to edit files that have large parts not concerned by us, but should not be modified, and it is not worth to encode a full I/O for it. 
    Its efficiency is a hell but still always faster than QM calculation so do not bother it.
    '''
    def __init__(self,filename=None):
        self.lines = []
        if (filename is not None):
            self.read(filename)

    def read(self,filename):
        '''
        Read a file
        '''
        f = open(filename,'r')
        self.lines = f.readlines()
        f.close()

    def write(self,filename):
        '''
        Write to a file
        '''
        f = open(filename,'w')
        for line in self.lines:
            f.write(line)
        f.close()

    def find_marked_line(self,mark,num_shift=0,ix_start=None):
        '''
        Get the index of line, with specific number of lines before/after its line

        :param mark: the string which the first appear should be found
        :param num_shift: the i-th line after the marked line(can be 0 or negative)
        :param ix_start: which line the search should start, start at 0
        :return: the line index, if not found return 0
        '''
        for i,line in enumerate(self.lines):
            if (ix_start is not None):
                if (ix_start > i):
                    continue
            if (mark in line):
                return i+num_shift
        return -1
    
    def __is_float__(self,st):
        '''
        Determine whether a string is float number
        '''
        st = st.strip()
        if (len(st) == 0):
            return False
        if (not st[0].isdigit() and st[0] != "-" and st[0] != "+" and st[0] != "."):
            return False

        for i in xrange(1,len(st)):
            if (not st[i].isdigit()): #Possibly A[d,e,E,D][+,-]B
                if (st[i] == "."):
                    continue
                if (st[i] in ("d","D","e","E")):
                    if not (st[i+1] == "+" or st[i+1] == "-" or st[i+1].isdigit()):
                        return False
                elif (st[i] == "+" or st[i+1] == "-"):
                    if (st[i-1] in ("d","D","e","E")):
                        return False
                return False

        return True


    def edit_value(self,ix_line,value,pos_end=None,name=None):
        '''
        Edit the value in the line, assuming only one value is presented in the line
        there are several way to guess it:
        1. after first "="
        2. first number 
        if pos_end is used then look for the position 
        if name is used then look for the name

        :param ix_line: the line index which to be edited
        :param value: the value which should be set
        :param pos_end: the end position of the value which is to be set
        :param name: the value just after name will be set ( possibly skip a "=" )
        '''
        line = self.lines[ix_line]
#Find the position of the number
        pos_start = -1 #Start of the number
        pos_space = 0 #Start of empty before the number
        if (pos_end is not None):
            pos_start = 0
            for i in xrange(pos_end,-1,-1):
                if (line[i]==" " or line[i]=="="): #Space or "=" split
                    pos_start = i+1
                    break
        else:
#Find name seperated
#pos_space start from just after "=" or the name
            pos_end = -1 #End of the number 
            if (name is not None):
                pos_space = line.find(name)
                if (pos_space == -1):
                    raise FormatError("Cannot find specific name %s" % name)
                pos_space += len(name)

#Find "="
            pos_space = line.find("=",pos_space)
            if (pos_space != -1):
                pos_space += 1
#Forward search value
                for i in xrange(pos_space,len(line)):
                    if (line[i] != " " and line[i] != "\t"):
                        pos_start = i
                        break
                if (pos_start == -1):
                    raise FormatError("Cannot find any number in line")
            else:
#No sign "="
#Find the first number
#pos_start start from the digit or "-xxx"
                pos_space = 0
                for i in xrange(pos_space,len(line)):
                    if (line[i].isdigit()):
                        pos_start = i
                        break
                if (pos_start != -1): #Found
#Check negative character
                    if (pos_start>1 and line[pos_start-1] == '-'):
                        pos_start -= 1
                else:
                    raise FormatError("Cannot find any number in line")
#Backward search empty space
                for i in xrange(pos_start-1,-1,-1):
                    if (line[i] != " "):
                        pos_space = i + 1
                        break
            #Find end
            b_inword = False
            pos_end = len(line)-1
            for i in xrange(pos_start,len(line)):
                if (b_inword and (line[i] == " " or line[i] == "\t")):
                    pos_end = i-1
                    break
                b_inword = (line[i] != " " and line[i] != "\t")

        #Use available spaces for editing
        #Try find floating format
        #Is it a floating number?
        st_org = line[pos_start:pos_end+1]
#Replace d,D,E to e
        st_org = st_org.replace("d","e").replace("D","e").replace("E","e")
        b_float = self.__is_float__(st_org)
        if (b_float):
#Try to find the format
            if ("e" in st_org):
#%e
                if ("." in st_org):
                    num_decimal = st_org.find("e")-st_org.find(".")
                else:
                    num_decimal = 1
                st_val = "%%%i.%ie" % (pos_end-pos_space+1,num_decimal)
                st_val = st_val % value
            else:
#%f
                num_decimal = pos_end-st_org.find(".")+1
                st_val = "%%%i.%if" % (pos_end-pos_space+1,num_decimal)
                st_val = st_val % value 
        else:
#Fill with new value
            st_val = str(value)

        if (len(st_val)>pos_end-pos_space):
            if (pos_end != len(line.rstrip())):
                raise FormatError("Not enough space in line")
            else:
                line = line[:pos_space] + st_val + "\n"
        else:
            line = line[:pos_end-len(st_val)+1] + st_val + line[pos_end+1:]

        self.lines[ix_line] = line

        return

    def edit_name_value(self,name,value):
        '''
        Edit name = value type string
        if name does not exist, then add one
        '''
        try:
            self.edit_value(self.find_marked_line(name),value=value,name=name)
        except FormatError:
#Not found, so we just add one
#Try to find where is "="
#           traceback.print_exc()
            pos_eq = [x.find("=") for x in self.lines]
            pos_max = max(itertools.groupby(sorted(pos_eq)),
                    key=lambda x, v:(len(list(v)),-pos_eq.index(x)))[0]
            pos_max = pos_max - 1

            if (len(name)>pos_max-1):
                name2 = name +" = "
            else:
                name2 = ("%%%is = " % pos_max) % name

            self.lines.append("%s%s\n" % (name2,str(value)))

        return

def f_Split_EnumerateRangeString(stRange,stRangeType="mms"):
    '''
    Generate a list of all possible combination from a string
    string format : PropertyA:20,40,50~100~10;PropertyB_x:4~12~2
    "_x" after name means relative value ( used for ecutrho, which is relative to ecutwfc)
    relative value will be set in the order of string, 
    if ecutrho_x is before ecutwfc, it will be used as ecutwfc in original input file * parameter 
    '''
    
    listProperty = []
    arList1 = stRange.split(";")
    for stProperty in arList1:
        (stName,stValue) = stProperty.split(":")
        listProperty.append([stName,f_Split_RangeString(stValue,stRangeType)])    
    
    result = []
    for (stName,aProperty) in listProperty:
        if ( result == []):
            for aValue in aProperty:
                result.append([[stName,aValue]])
        else:
            result2 = []
            for aValue in aProperty:
                for aPart in result:
                    aPart2 = copy.deepcopy(aPart)
                    aPart2.append([stName,aValue])
                    result2.append(aPart2)
            result = result2
#   print(result)         
    return result 


def f_Split_RangeString(stRange,stRangeType="mms"):
    '''
    Generate a list of number from string.
    String format: a,b,c,d~e~f,h,... ( every character represent a float value)
    ',' seperate different value
    'a~b~c' define a list of a,a+c,a+2*c,...a+n*c, until a+n*c > b ( if a <b) or a+n*c < b ( if a > b)
    which means 1~3~1 -> 1,2,3
    if a contains non-digit character, the value will be stored as string
    if "~" is used, the value will be stored as float; otherwise string
    
    :param stRange: the string you wanna to split
    :param stRangeType: define what a~b~c represent , 'mms' means min-max-step, 'msm' means min-step-max
    '''
    arList = stRange.split(',')
    result = []
    for aPart in arList:          
        arList2 = aPart.split("~")
        if ( len(arList2) == 3):
            if ( stRangeType == "mms"):
                (fStart,fEnd,fStep) =[float(x) for x in arList2]
                #(fStart,fEnd,fStep) =arList2
            else:
                (fStart,fStep,fEnd) =[float(x) for x in arList2]
                #(fStart,fEnd,fStep) =arList2
            if ( arList2[2] == "0" ):
                raise ZeroDivisionError("Cycle Step cannot be 0!")
            #Always use float for this 
            #
            i = float(fStart)
            d1 = fEnd-fStart
            if (abs(d1) <= 1e-8):
                raise ValueError("Start/End value must be different : %s" % aPart)
            while ( (fEnd-i)*d1 >= 0 ):
                result.append(i)
                i = i + fStep
        else:
            #if ( aPart.replace("-","").replace(".","").isdigit()):
            #    #result.append(float(aPart))
            #    result.append(aPart)
            #else:
            #    result.append(aPart)
            result.append(aPart)
    return result

def f_SolveCross(s1,s2,s3):
    '''
    Find minimum integer solution for UxV = nS
    '''
    (u1,u2,u3,v1,v2,v3,n) = (0,0,0,0,0,0,1)
    #test 0-9 for u-n
    result = []
    nMax = 3
    for n in range(1,2):
        for u1 in range(-nMax,nMax):
            for u2 in range(-nMax,nMax):
                for u3 in range(-nMax,nMax):
                    for v1 in range(-nMax,nMax):
                        for v2 in range(-nMax,nMax):
                            for v3 in range(-nMax,nMax):
                                if ( u2*v3-u3*v2 == n*s1 and u3*v1-u1*v3 == n*s2 and u1*v2-u2*v1 == n*s3 ):
                                    result.append((abs(u1)+abs(u2)+abs(u3)+abs(v1)+abs(v2)+abs(v3),[u1,u2,u3,v1,v2,v3,n] ))
    result.sort()
    return result



def f_MaxCommonSubmultiple(a, b):
    '''
    return maximum common submultiple of two integer
    '''
    if a == 0:
        return b
    else:
        return f_MaxCommonSubmultiple(b % a, a)

def f_MinCommonMultiple(a,b):
    '''
    return minimum common multiple of two integer
    '''
    return (a*b/f_MaxCommonSubmultiple(a,b))

def f_GetSurfaceFromMiller(list_index):
    '''
    Get three vectors ( base cell vectors ) perpendicular to specific miller index
    :return: first two make the surface , the third make it a cell 
    '''
    if ( len(list_index) != 3):
        raise ValueError("Millier index must has 3 integer")
    #detect if 0 exist
    #b0 = 0 in list_index
    index = [0,1,2]
    #revIndex = []
    result = []
    
    #select non-0 value in the miller index
    #select a vector from one point of miller surface cross point with axis plane to another
    #(this line is useless now)note due to symmetry, (1,0,0) -> (0,1,0) = -1,1,0 = 1,1,0, so only positive is used 
    #if one of its index is 0, use 1 as vector, otherwise use min common multiple
    for i in range(0,3):
        if ( list_index[i] != 0):
            break
    if ( i == 2):
        index = [2,0,1]
    elif ( i == 1):
        index = [1,0,2]
    for i in range(1,3):
        base = [0,0,0] 
        if ( list_index[index[i]] == 0):
            base[index[i]] = 1
        else:
            CM = f_MinCommonMultiple(list_index[index[0]],list_index[index[i]])
            base[index[0]] = CM / list_index[index[0]]
            base[index[i]] = -CM / list_index[index[i]]
        result.append(base)
    # find c axis has minimum number
    #try linear combiniation to find lowest 
    #dMin = lu.f_List_dot3(result[0], result[1])
    #listLM = copy.deepcopy(result)
    #for i in range(-1,1):
    #    for j in range(-1,1):
    #        vA = lu.f_List_Op_List(result[0],'+',lu.f_List_Op_Scalar(result[1],"*",i))
    #        vB = lu.f_List_Op_List(result[1],'+',lu.f_List_Op_Scalar(result[0],"*",j))
    #        if ( lu.f_List_IsParallel(vA, vB)):
    #            continue
    #        dV = lu.f_List_norm3(lu.f_List_cross3(vA, vB))
    #        if ( dV < dMin - 0.00001 ):
    #            print(i,j,dV)
    #            dMin = dV
    #            listLM = [vA,vB]
    #result = listLM
    
    result.append ( [int(x) for x in lu.f_List_cross3(result[0], result[1])]) 
    return result
    

def f_RotateMatrix(dAngle,vec):
    '''
    Create a rotation matrix from angle and vector ( vector will be convert to unit vector )
    '''
    vRot = lu.f_List_normalize(vec)
    mRot = lu.f_Matrix_zeros(3,3)
    cosA = math.cos(dAngle)
    sinA = math.sin(dAngle)
    mRot[0][0] = cosA + ( 1-cosA)*vRot[0]*vRot[0]
    mRot[0][1] = (1-cosA)*vRot[0]*vRot[1]-sinA*vRot[2]
    mRot[0][2] = (1-cosA)*vRot[0]*vRot[2]+sinA*vRot[1]
    mRot[1][0] = (1-cosA)*vRot[1]*vRot[2]+sinA*vRot[0]
    mRot[1][1] = cosA + (1-cosA)*vRot[1]*vRot[1]
    mRot[1][2] = (1-cosA)*vRot[1]*vRot[2]-sinA*vRot[0]
    mRot[2][0] = (1-cosA)*vRot[2]*vRot[0]-sinA*vRot[1]
    mRot[2][1] = (1-cosA)*vRot[0]*vRot[1]+sinA*vRot[2]
    mRot[2][2] =  cosA + (1-cosA)*vRot[2]*vRot[2]    
    return mRot

def ConvertLatticeToCartesian(fLatticeParameter,fLatticeAngle):
    '''
    Convert lattice parameter+angle to Bravis matrix( Numpy array )
    
    :param fLatticeAngle: should in unit of degree
    '''
    fAngleRad = [x/180.0*math.pi for x in fLatticeAngle]
    #print(fLatticeAngle)
    #print(fAngleRad)
    
    #fCosGamma0 = (math.cos(fAngleRad[2])-math.cos(fAngleRad[0])*math.cos(fAngleRad[1]))/math.sin(fAngleRad[0])/math.sin(fAngleRad[1])
    fCosGamma0 = (math.cos(fAngleRad[0])-math.cos(fAngleRad[2])*math.cos(fAngleRad[1]))/math.sin(fAngleRad[2])/math.sin(fAngleRad[1])
    
    #print(fCosGamma0)
    
    fGamma0 = math.acos(fCosGamma0)
    
    #print(fGamma0*180.0/math.pi)
    
    result = lu.f_Matrix_zeros(3,3)
    
    #set a-axis as [a,0,0], a-b in a plane of axis
#    result[2,2] = math.sin(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[2]
#    result[2,1] = math.cos(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[2]
#    result[2,0] = math.cos(fAngleRad[1])*fLatticeParameter[2]
#    result[1,2] = 0
#    result[1,1] = math.sin(fAngleRad[0])*fLatticeParameter[1]
#    result[1,0] = math.cos(fAngleRad[0])*fLatticeParameter[1]
#    result[0,2] = 0
#    result[0,1] = 0
#    result[0,0] = 1*fLatticeParameter[0]
    result[2][2] = math.sin(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[2]
    result[2][1] = math.cos(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[2]
    result[2][0] = math.cos(fAngleRad[1])*fLatticeParameter[2]
    result[1][2] = 0
    result[1][1] = math.sin(fAngleRad[2])*fLatticeParameter[1]
    result[1][0] = math.cos(fAngleRad[2])*fLatticeParameter[1]
    #print(math.cos(fAngleRad[2]),fAngleRad[2]*180.0/math.pi)
    result[0][2] = 0
    result[0][1] = 0
    result[0][0] = 1*fLatticeParameter[0]

    
    #set c-axis as [0,0,c] b-c in a plane of axis
#    result[0,0] = math.sin(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[0]
#    result[0,1] = math.cos(fGamma0)*math.sin(fAngleRad[1])*fLatticeParameter[0]
#    result[0,2] = math.cos(fAngleRad[1])*fLatticeParameter[0]
#    result[1,0] = 0
#    result[1,1] = math.sin(fAngleRad[0])*fLatticeParameter[1]
#    result[1,2] = math.cos(fAngleRad[0])*fLatticeParameter[1]
#    result[2,0] = 0
#    result[2,1] = 0
#    result[2,2] = 1*fLatticeParameter[2]
    
    
    
    #print(result)
    
    return result

def ConvertCartesianToLattice(listBravis,nLatticeType=0):
    '''
    Convert Bravis matrix to 3 parameter + 3 angle

    :param listBravis: Bravis matrix
    :param nBravis: represent bravis type, will convert primitive lattice to convetional lattice according to this.
    '''
    fBravis = listBravis
    result = fBravis
    resultAngle = [0.0,0.0,0.0]
  
    #print(numpy.dot(result[1,:],result[2,:]))
    #print(lu.f_List_norm3(result[1,:]))
    #print(fBravis[:,0])
    
    resultAngle[0] = lu.f_List_angle3(result[1],result[2])
    resultAngle[1] = lu.f_List_angle3(result[2],result[0])
    resultAngle[2] = lu.f_List_angle3(result[0],result[1])
    resultAngleDeg = [x / math.pi * 180 for x in resultAngle]
    
    resultLength = [0,0,0]
    for i in range(0,3):
        resultLength[i] = lu.f_List_norm3(fBravis[i])        
    # centered lattice; this function has been removed into Lattice class
#    if ( nLatticeType == 2 or nLatticeType == 10 ): # fc
#        resultAngleDeg = [90.0,90.0,90.0]
#        resultLength[0] = result[0,0]*2
#        resultLength[1] = result[1,1]*2
#        resultLength[2] = result[2,2]*2
#    elif ( nLatticeType == 3 or nLatticeType == 7 or nLatticeType == 11 ): # bc
#        resultAngleDeg = [90.0,90.0,90.0]
#        resultLength[0] = abs(result[0,0]*2)
#        resultLength[1] = abs(result[0,1]*2)
#        resultLength[2] = abs(result[0,2]*2)
#    elif ( nLatticeType == 9): # boc ortho
#        resultAngleDeg = [90.0,90.0,90.0]
#        resultLength[0] = abs(result[0,0]*2)
#        resultLength[1] = abs(result[0,1]*2)
#        resultLength[2] = abs(result[2,2]*2)        
#    elif ( nLatticeType == 13 ): #boc monoclinic
#        #b = lu.f_List_norm3(numpy.matrix(result[1,:]))
#        b = lu.f_List_norm3(result[1])
#        resultAngleDeg = [90.0,90.0,math.acos(result[1,0]/b)/math.pi * 180]
#        resultLength[0] = abs(result[0,0]*2)
#        #resultLength[1] = lu.f_List_norm3(numpy.matrix(result[1,:]))
#        resultLength[1] = lu.f_List_norm3(result[1])
#        resultLength[2] = abs(result[0,2]*2)
#    elif ( nLatticeType == 5):#R
#        resultAngleDeg = [90.0,90.0,120]
#        resultLength[0] = abs(result[0,0]*2)
#        resultLength[1] = abs(result[0,1]*math.sqrt(3)*2)
#        resultLength[2] = abs(result[0,2]*3)        
        
               
    return resultLength,resultAngleDeg


def ConvertInnerCoordToCartesian(fCell,listInner):
    '''
    Convert inner coordinate to cartesian coordinate in specific lattice
    '''
    fBravis = fCell
    fInner = listInner
    
    result = lu.f_Matrix_dot(fBravis,fInner)
    
    return result

def ConvertCartesianToInnerCoord(fCell,listCart):
    '''
    Convert cartesian to inner coordinate in specific lattice

    :param fCell: a 3x3  ndarray, each row is a vector of lattice axis
    :param listCart: a list of Cartesian coordinate
    '''
    fBravis = fCell
    fCart = listCart
    
    result = lu.f_Matrix_dot(lu.f_List_inv3(fBravis),lu.f_Matrix_transpose(fCart)) # convert fCart to column vector
    return lu.f_Matrix_transpose(result) # return normal list from row vector, remove dual list
    

def f_Latt_IsClinic(listAngle):
    '''
    Determine whether a cell is monoclinic or tricilinic
    '''
    if ( abs(listAngle[0]-listAngle[1]) < 0.01 and abs(listAngle[0]-listAngle[2])<0.01):
#Rho
        return False
    for dAngle in listAngle:
#Ortho or Hex
        if ( abs(dAngle-90.0) > 0.01 and abs(dAngle-120.0) > 0.01):
            return True
    return False


class Lattice(object):
    '''
    A basic representation of lattice, all basic number unit in bohr and deg
    Cartesian vectors follow this: c-axis as [0,0,c] b-c in a plane of axis
    Members start with "Create" will return a new Lattice object, otherwise itself will change
    '''
    
    coord_Cartesian = 1
    coord_Primitive = 2
    coord_Conventional = 3
    coord_Raw = 4
    
    coord_Cart = 1
    coord_Prim = 2
    coord_Conv = 3
    
    dErrMax = 0.000001 # using to <,=,> of float type

    filename_xsd = "lattice.xsd"
    
    def __init__(self):
        self.nLatticeType = 0
        self.fLatticeLength = [0,0,0]
        self.fLatticeAngle = [0,0,0]
        self.listAtom = [] # in internal coordinate of primitive cell
        self.__bRaw__ = False # in raw mode, raw cell vector will be used to determine atom position. And primitive-conventional conversion will be banned. Don't change it in other function !
        self.__RawPrimitiveVector__ = [[1,0,0],[0,1,0],[0,0,1]] # raw cell vector, maybe deviate from autocreate vector like [1,0,0]/[0,1,0] but [/2,/2,0] and [-/2,/2,0]  
        self.__RawConventionalVector__ = [[1,0,0],[0,1,0],[0,0,1]]
        #self.PrimitiveCellVector = []
        self.listFix = [] #: list of where atoms are fixed on given direction

    @classmethod
    def load_xml(cls,filename):
        '''
        Load parameters from a XML file
        '''
        return XMLSerilizer.load_xml(globals(),cls,filename)
 
    def save(self,filename):
        '''
        Save parameters to a XML file
        '''
        XMLSerilizer.save_xml(Lattice,self,filename)

    @property
    def IsRaw(self):
        '''
        Return whether Lattice object is operated in raw mode
        '''
        return self.__bRaw__
    
    @property
    def PrimitiveCellVector(self):
        '''
        Primitive cell vectors
        '''
        if ( self.__bRaw__):
            return self.__RawPrimitiveVector__
        else:
            return self.__GetPrimitiveCell__()
    
    @property    
    def ConventionalCellVector(self):
        '''
        Conventional cell vectors
        '''
        if ( self.__bRaw__):
            return self.__RawConventionalVector__
        else:
            return ConvertLatticeToCartesian(self.fLatticeLength,self.fLatticeAngle)    

    @property
    def Volume(self):
        '''
        Volume of convetional cell vectors
        '''
        if ( self.__bRaw__):
            vt = zip(*self.__RawConventionalVector__)
            #print(vt)
            return lu.f_List_dot3(lu.f_List_cross3(vt[0],vt[1]),vt[2])
        else:
            #print(self.fLatticeAngle)
            dCos = [math.cos(x/180.0*math.pi) for x in self.fLatticeAngle]
            #print(dCos)
            return self.fLatticeLength[0]*self.fLatticeLength[1]*self.fLatticeLength[2]*math.sqrt(1+2*dCos[0]*dCos[1]*dCos[2]-dCos[0]**2-dCos[1]**2-dCos[2]**2)

    def init_fix(self,b_fix=False):
        '''
        Create informations about fix

        :param b_fix: initialze with True or False
        '''
        self.listFix = [ [b_fix,b_fix,b_fix] for x in self.listAtom]
       
    def ReadFromABC(self,nLatticeType,fLatticeLength,fLatticeAngle,unit_length="bohr",unit_angle="deg"):
        if ( unit_length == "ang" or unit_length == "angstrom"):
            fLatticeLength = [x * constants.Ang2Bohr for x in fLatticeLength]
        if ( unit_angle == "rad"):
            fLatticeAngle  = [x * constants.Rad2Deg for x in fLatticeAngle]
        self.nLatticeType = nLatticeType
        self.fLatticeLength = fLatticeLength
        self.fLatticeAngle = fLatticeAngle
        #self.PrimitiveCellVector = self.__GetPrimitiveCell__()
       
        #print("Before Primitive:")
        #self.WriteToFile("")
        
        #self.__GetLatticeFromPrimitiveCell__()
        self.__GetLatticeFromPrimitiveCellVector__(self.PrimitiveCellVector)
        
        #print("After Primitive")
        #self.WriteToFile("")
    
    def ReadFromRawCellVector(self,mBravisConv,mBravisPrim=None):
        '''
        Read Lattice object from specified cell vectors
        Note one may provide BOTH conventional and primitive cell vectors

        :param mBravisConv: Bravis matrix of the conventional cell
        :param mBravisPrim: Bravis matrix of the primitive cell, optional; If negelected Conv/Prim are the same.
        '''
        if (mBravisPrim is None):
            mBravisPrim = copy.deepcopy(mBravisConv)
        self.__RawPrimitiveVector__ = copy.deepcopy(mBravisPrim)
        self.__RawConventionalVector__ = copy.deepcopy(mBravisConv)
        self.__bRaw__ = True
#Create lattice length/angle information
        (self.fLatticeLength,self.fLatticeAngle) = ConvertCartesianToLattice(self.__RawConventionalVector__)

    def ReadFromPrimitiveCellVector(self,mBravis,nLatticeType=None):
        '''
        Read Cell from cell vectors and Lattice Type ( If no lattice type is used, original one will be used)
        Atoms' crystal coordination will not moved ! 

        :param mBravis: Cell vectors
        :param nLatticeType: The lattice type, if NONE then the nLatticeType stored in the current object is used 
        :param bRaw: If set to true, will not utilize primitive-conventional convertion in the process, and put the mode of Lattice object to raw
        '''
        if ( nLatticeType != None):
            self.nLatticeType = nLatticeType
        self.__GetLatticeFromPrimitiveCellVector__(mBravis)
    

    def __GetPrimitiveCell__(self):
        
        (a,b,c)=self.fLatticeLength
        (aa,ab,ac)= [x /180.0 * math.pi for x in self.fLatticeAngle]
        #print(a,b,c,aa,ab,ac)
        
        if ( aa == 0.0):
            aa = math.pi / 2.0
        if ( ab == 0.0):
            ab = math.pi / 2.0
        if ( ac == 0.0):
            ac = math.pi / 2.0
        
        listPrim = []
        listPrim.append([[a,0,0],[0,a,0],[0,0,a]]) # Nothing
        listPrim.append([[a,0,0],[0,a,0],[0,0,a]]) # P cubic
        listPrim.append([[-a/2,0,a/2],[0,a/2,a/2],[-a/2,a/2,0]]) # F
        listPrim.append([[a/2,a/2,a/2],[-a/2,a/2,a/2],[-a/2,-a/2,a/2]]) # I
        listPrim.append(([a,0,0],[-a/2,math.sqrt(3)/2*a,0],[0,0,c])) # H simple
        tx = sqrt((1-cos(aa))/2)*a
        ty = sqrt((1-cos(aa))/6)*a
        tz = sqrt((1+2*cos(aa))/3)*a
        listPrim.append([[tx, -ty, tz],[0, 2*ty, tz], [-tx, -ty, tz]]) #R (+5)
        listPrim.append([[a,0,0],[0,a,0],[0,0,c]]) # P Tetra
        listPrim.append([[a/2,-a/2,c/2],[a/2,a/2,c/2],[-a/2,-a/2,c/2]]) # I Tetra
        listPrim.append([[a,0,0],[0,b,0],[0,0,c]]) # P Ortho
        listPrim.append([[a/2,b/2,0],[-a/2,b/2,0],[0,0,c]]) # B-center Ortho
        listPrim.append([[a/2,0,c/2],[a/2,b/2,0],[0,b/2,c/2]]) # F Ortho
        listPrim.append([[a/2,b/2,c/2],[-a/2,b/2,c/2],[-a/2,-b/2,c/2]]) # I Ortho
        listPrim.append([[a,0,0],[b*math.cos(ac),b*math.sin(ac),0],[0,0,c]]) # monoclinic
        listPrim.append([[a/2,0,-c/2],[b*math.cos(ac),b*math.sin(ac),0],[a/2,0,c/2]]) # B-c mono
        listPrim.append([[a,0,0],[b*math.cos(ac),b*math.sin(ac),0],[c*math.cos(ab),  c*(math.cos(aa)-math.cos(ab)*math.cos(ac))/math.sin(ac), c*math.sqrt( 1 + 2*math.cos(aa)*math.cos(ab)*math.cos(ac)  - math.cos(aa)**2-math.cos(ab)**2-math.cos(ac)**2 )/math.sin(ac)]]) # Tri
        
        #print("Primitive Cell: %d" % self.nLatticeType)
        #print(listPrim[self.nLatticeType])
        
        return listPrim[self.nLatticeType]
    
    def __GetLatticeFromPrimitiveCellVector__(self,mBravis):
        '''
        Get QE-type described ibrav from bravis matrix
        Note some matrix input here is not from QE and may be in other format!
        @todo trigonal,monoclinic and triclinic are different in w2k and QE
        ibrav -5 and -12 is not supported
        '''
        #PV = self.PrimitiveCellVector
        PV = mBravis
        #print(PV)
        a0 = lu.f_List_norm3(PV[0]) # for a axis = [a,0,0] cell to minimize non-symmetry influence
        b0 = lu.f_List_norm3(PV[1])
        c0 = lu.f_List_norm3(PV[2])
        a = PV[0][0]
        b = PV[1][1]
        c = PV[2][2]
        if ( a < 0 ):
            a = -a
        if ( c < 0 ):
            c = -c
        #Calc primitive cell information    
        (self.fLatticeLength,self.fLatticeAngle) = ConvertCartesianToLattice(PV)    
        #print(self.fLatticeLength,self.fLatticeAngle)
        #The ibrav has more than 1 unit, need listLatt
        arMulti = [2,3,5,7,9,10,11,13]
        arNonSpecialAngle = [12,13,14] # non-90 or 60 

#If it is monoclinic or triclinic then directly return
        if ( self.nLatticeType in arNonSpecialAngle ):
            return
        
        listLatt= []
        listLatt.append([a0,a0,a0,90.0,90.0,90.0])
        listLatt.append([a0,a0,a0,90.0,90.0,90.0])
        listLatt.append([a*2,a*2,a*2,90.0,90.0,90.0])
        listLatt.append([a*2,a*2,a*2,90.0,90.0,90.0])
        listLatt.append([a0,a0,c,90.0,90.0,120.0])
        listLatt.append([a*2,a*2,c*3,90.0,90.0,120.0])
        listLatt.append([a0,a0,c,90.0,90.0,90.0])
        listLatt.append([a*2,a*2,c*2,90.0,90.0,90.0])
        listLatt.append([lu.f_List_norm3(PV[0]),lu.f_List_norm3(PV[1]),lu.f_List_norm3(PV[2]),90.0,90.0,90.0])
        listLatt.append([a*2,b*2,c,90.0,90.0,90.0])
        listLatt.append([a*2,b*2,c*2,90.0,90.0,90.0])
        listLatt.append([a*2,b*2,c*2,90.0,90.0,90.0])
        #if ( PV[0][1] != 0 or PV[0][2] != 0 or  PV[1][0] != 0 or PV[1][2] != 0 or PV[2][0] != 0 or PV[2][1] != 0): 
        if ( False):
            # only deal with QE-type monoclinic and triclinic when the bravis matrix is not belong to orthogonal to avoid 0 error
            listLatt.append([a,lu.f_List_norm3(PV[1]),c,90.0,90.0,180.0/math.pi*math.atan(PV[1][1]/PV[1][0])])
            listLatt.append([a*2,lu.f_List_norm3(PV[1]),c*2,90.0,90.0,180.0/math.pi*math.atan(PV[1][1]/PV[1][0])])
        else:
            listLatt.append([-1,-1,-1,-1,-1,-1])
            listLatt.append([-1,-1,-1,-1,-1,-1])
        listLatt.append(self.fLatticeLength + self.fLatticeAngle)
        
        #if ( self.nLatticeType in arMulti):
        if ( listLatt[0] == -1):
#If not found, than use direct calculated result, which is already set
            pass
        else:#Use detected result
           # Indeed, some results can be obtained from direct calculation, but we prefer manually classified ones to avoid inconsistence in some non-sysmmetry-keep optimization ( like BFGS in Espresso)
            self.fLatticeLength = listLatt[self.nLatticeType][0:3]
            self.fLatticeAngle = listLatt[self.nLatticeType][3:6]
 
        

    def __ConvertCoordinate__(self,nIn,nOut,arPos):
        '''
        Convert a kind of coordinate to another
        '''
        mConventionalCellT = lu.f_Matrix_transpose(self.ConventionalCellVector)
        mPrimitiveCellT = lu.f_Matrix_transpose(self.PrimitiveCellVector)
        #mRawCellT = lu.f_Matrix_transpose(self.__RawCellVector__)
        if ( nIn == Lattice.coord_Conventional):
            mPosC =  lu.f_Matrix_dot(mConventionalCellT , lu.f_Matrix_transpose([arPos]))
        elif ( nIn == Lattice.coord_Primitive):
            mPosC = lu.f_Matrix_dot(mPrimitiveCellT , lu.f_Matrix_transpose([arPos]))
        elif ( nIn == Lattice.coord_Cartesian):
            mPosC = lu.f_Matrix_transpose([arPos])
#        elif ( nIn == Lattice.coord_Raw):
#            mPosC = lu.f_Matrix_dot(mRawCellT,lu.f_Matrix_transpose([arPos]))
        else:
            print("Unsupported input coordinate type: " + str(nIn))
        #print("arPos: ",arPos)
        #print("mPosC: ",mPosC)    
        if ( nOut == Lattice.coord_Conventional ):
            mOut = lu.f_Matrix_dot(lu.f_List_inv3(mConventionalCellT) , mPosC)
        elif ( nOut == Lattice.coord_Primitive):
            mOut = lu.f_Matrix_dot( lu.f_List_inv3(mPrimitiveCellT) , mPosC)
        elif ( nOut == Lattice.coord_Cartesian):
            mOut = mPosC
#        elif ( nOut == Lattice.coord_Raw):
#            mOut = lu.f_Matrix_dot( lu.f_List_inv3(mRawCellT) , mPosC)
        else:
            print("Unsupported output coordinate type: " + str(nIn))
        
        #print("mOut",lu.f_Matrix_transpose(mOut)[0])
        return lu.f_Matrix_transpose(mOut)[0]        
    
    def AddAtom(self,listPos,unit="bohr",latt="primitive",a=None):
        '''
        Add a atom list by Cartesian Coordinate or Crystal Coordinate, indicated by unit_length
        Note: please be cautious with latt="alat", because "a" is not fixed during cell shape change.
        
        :param listCart: in form of  [Element name,x,y,z]
        :param unit: indicate type of xyz : bohr, ang/angstrom, cry/crys/crystal, alat
        :param latt: indicate which cell xyz is related to : primtive/conventional ( only affect when unit is crystal ), raw ( stored cell vectors )
        :param a: the length of a in unit of Bohr. Default is the fLatticeLength (which is the first lattice vector if bRaw, otherwise conventional lattice vector).
        '''
        if ( unit == "bohr"):
            self.AddAtomFromCartesian(listPos)
        elif ( unit == "ang" or unit == "angstrom"):
            listPos2 = [ [ aPos[0] ] + [x*constants.Ang2Bohr for x in aPos[1:4]] for aPos in listPos ]
            self.AddAtomFromCartesian(listPos2)
        elif ( unit == "alat"):
            if (a is None):
                a = self.fLatticeLength[0]
            listPos2 = [ [ aPos[0] ] + [x*a for x in aPos[1:4]] for aPos in listPos ]
            self.AddAtomFromCartesian(listPos2)
        elif ( unit == "cry" or unit == "crys" or unit =="crystal" or unit == "internal"):
            if ( latt == "primitive" or latt == "prim" ):#internal coordinate of Primitive Cell
                for fCoord in listPos:
                    self.listAtom.append(fCoord)
            elif ( (latt == "normal") or (latt == "conv") or (latt=="convetional")): #internal coordinate of Normal cell
                #ConventionalCell = numpy.matrix(ConvertLatticeToCartesian(self.fLatticeLength,self.fLatticeLength))
                for fCoord in listPos:
                    #self.listAtom.append([fCoord[0]]+self.__GetPrimitiveCoordFromNormalCoord__(fCoord[1:4]))
                    self.listAtom.append([fCoord[0]]+self.__ConvertCoordinate__(Lattice.coord_Conventional,Lattice.coord_Primitive,fCoord[1:4]))
                #self.AddAtomFromCartesian([fCoord[0]]+ConventionalCell.T*numpy.matrix(fCoord[1:4]).tolist())
            else:
                print("AddAtom: Unrecognized mode '%s'" % unit)                    
        else:
            print("AddAtom: Unrecognized mode '%s'" % unit)
                             

    
    def AddAtomFromCartesian(self,listCart):
        '''
        Add a atom list by Cartesian Coordinate

        :param listCard: a list of atom positions, each atom positions must be in the format ['Name',x,y,z]
        '''
        for fCart in listCart:
            fInner = self.__ConvertCoordinate__(Lattice.coord_Cartesian, Lattice.coord_Prim, fCart[1:4])
            #print(fCart,fInner)
            #fInner = (numpy.matrix(fPrim).I*numpy.matrix(fCart[1:4]).T).T.tolist()
            #fInner = ConvertCartesianToInnerCoord(fPrim,fCart[1:4])
            self.listAtom.append( [fCart[0]]+fInner )
         
    def CorrectError (self):
        '''
        Set some standard fractional coordinate (like 1/3 ) to their exact value to beautify the output
        Example : 0.6666663562472133 to 0.6666666666667, 0.249999998232 to 0.25
        '''
        for i in range(0,3):
            self.fLatticeAngle[i] = f_CorrectFloat(self.fLatticeAngle[i])
        for aAtom in self.listAtom:
            for i in range(1,4):
                aAtom[i] = f_CorrectFloat(aAtom[i])            
                

    def ShowSummary(self,stFileName=None):
        '''
        Write cell information to specific file or screen

        :param stFileName: if set to None, set print to screen
        '''
        if ( stFileName != None):
            fOut = open(stFileName,"w")
        else:
            fOut = sys.stdout
            
        print("Summary of cell:",file=fOut)
        print(self.fLatticeLength,file=fOut)
        print(self.fLatticeAngle,file=fOut)
        print(self.PrimitiveCellVector,file=fOut)
        for aAtom in self.listAtom:
            print(aAtom,file=fOut)
        print("\n",file=fOut)
        if ( stFileName != None):
            fOut.close()
            
    @property        
    def listAtomInConvCell(self):
        '''
        List all atom in the conventional cell ( by conventional cell crystal coordinate )
        '''
#        listAtomNormal = []
#        for aAtom in self.listAtom:
#            for a1 in range(-1,2): #Extend to 3x3 primitive cell and list all atom in the normal cell
#                for a2 in range(-1,2):
#                    for a3 in range(-1,2):
#                        #arPos = self.__GetNormalCoordFromPrimitiveCoord__([aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
#                        arPos = self.__ConvertCoordinate__(Lattice.coord_Primitive,Lattice.coord_Conventional,[aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
#                        if ( arPos[0] < 1 and arPos[0] >= 0 and arPos[1] < 1 and arPos[1] >= 0 and arPos[2] < 1 and arPos[2] >=0):
#                            listAtomNormal.append([aAtom[0]] + f_CorrectFloat(arPos))
#        return listAtomNormal
        return self.GetAtomList("conv", "conv")
    
    @property
    def listAtomCartesian(self):
        '''
        list all atom in Cartesian Coordinate
        '''
#        listAtom = []
#        for aAtom in self.listAtom:
#            arPos = self.__ConvertCoordinate__(Lattice.coord_Primitive, Lattice.coord_Cartesian, aAtom[1:4])
#            listAtom.append([aAtom[0]]+arPos)
#        return listAtom
        return self.GetAtomList()
    
    def GetAtomList(self,unit="bohr",latt="primitive",extend=[0,0,0,0,0,0]):
        '''
        Get a list of atom in some representaion

        :param unit: the unit of coordinate, including bohr,ang,alat,primitive,conventional
        :param latt: the atom in whole convetional cell or one primitive cell
        :param extend: the lattice will be examed 
        '''
        listAtom1 = []
        #Get atoms
        
        bPrim = True
        #fMin = [ -x for x in extend] # lower boundary
        #fMax = [ x+1 for x in extend]    #upper boundary
        fMin = 0
        fMax = 1
           
        if ( latt == "normal" or latt == "conv" or latt == "convetional"):
            #extend = [ x * 2 + 1 for x in extend]
            bPrim = False
        elif ( latt == "primitive" or latt == "prim"):
            pass 
        else:
            raise ValueError("Unknown cell type: %s" % latt)
        
        arConvTrans = [[],[],[]]
        listAtom2 = []
        if ( not bPrim ): # get convtional cell transitional move vector ( [a,b,c], in unit of primitive cell )
            arConvTrans[0] = self.__ConvertCoordinate__(Lattice.coord_Conventional,Lattice.coord_Primitive,[1.0,0.0,0.0])
            arConvTrans[1] = self.__ConvertCoordinate__(Lattice.coord_Conventional,Lattice.coord_Primitive,[0.0,1.0,0.0])
            arConvTrans[2] = self.__ConvertCoordinate__(Lattice.coord_Conventional,Lattice.coord_Primitive,[0.0,0.0,1.0])
            #print(arConvTrans)
            #Get atoms in 1x1x1 conv cell
            for aAtom in self.listAtom:
                for a1 in range(-1,2): #Extend to 3x3 primitive cell and list all atom in the conventional cell
                    for a2 in range(-1,2):
                        for a3 in range(-1,2):
                            arPos = self.__ConvertCoordinate__(Lattice.coord_Primitive,Lattice.coord_Conventional,[aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
                            arPos = f_CorrectFloat(arPos)
                            #print(arPos)
                            bInCell = True
                            for i in range (0,3):
                                if ( arPos[i] >= fMax or arPos[i] < fMin ):
                                    #print(arPos)
                                    bInCell = False
                                    break
                            if ( bInCell):
                                #print(a1,a2,a3)
                                listAtom2.append([aAtom[0],aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
        else:
            listAtom2 = self.listAtom    
        #print(listAtom2)
        for aAtom in listAtom2:
            for a1 in range(extend[0],extend[1]+1): #Extend to 3x3 primitive cell and list all atom in the conventional cell
                for a2 in range(extend[2],extend[3]+1):
                    for a3 in range(extend[4],extend[5]+1):
                        #arPos = self.__GetNormalCoordFromPrimitiveCoord__([aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
                        if ( bPrim ):
                            listAtom1.append([aAtom[0],aAtom[1]+a1,aAtom[2]+a2,aAtom[3]+a3])
                        else:
                            arPos = aAtom[1:4]
                            #arPos = lu.f_List_Op_List(arPos, "+", lu.f_List_Op_Scalar(arConvTrans[0], "*", a1))
                            #arPos = lu.f_List_Op_List(arPos, "+", lu.f_List_Op_Scalar(arConvTrans[1], "*", a3))
                            #arPos = lu.f_List_Op_List(arPos, "+", lu.f_List_Op_Scalar(arConvTrans[2], "*", a3))
                            arPos = lu.f_List_Op_List(arPos,"+",lu.f_Matrix_transpose(lu.f_Matrix_dot(lu.f_Matrix_transpose(arConvTrans),[[a1],[a2],[a3]]))[0])
                            listAtom1.append([aAtom[0]] + arPos)
                            continue
        
        #Convert to coordinate
        nIn = Lattice.coord_Primitive
        listAtom2 = []
        if ( unit == "bohr"):
            for aAtom in listAtom1:
                listAtom2.append([aAtom[0]]+self.__ConvertCoordinate__(nIn, Lattice.coord_Cartesian, aAtom[1:4]))
        elif ( unit == "ang" or unit == "angstrom"):
            for aAtom in listAtom1:
                listAtom2.append([aAtom[0]]+[x*constants.Bohr2Ang for x in self.__ConvertCoordinate__(nIn, Lattice.coord_Cartesian, aAtom[1:4])])
        elif ( unit == "alat"):
            for aAtom in listAtom1:
                listAtom2.append([aAtom[0]]+[x / self.fLatticeLength[0] for x in self.__ConvertCoordinate__(nIn, Lattice.coord_Cartesian, aAtom[1:4])])
        elif ( unit == "primitive" or unit == "prim"):
            for aAtom in listAtom1:
                listAtom2.append([aAtom[0]]+ f_CorrectFloat(aAtom[1:4]))
        elif ( unit == "conventional" or unit == "conv"):
            for aAtom in listAtom1:
                listAtom2.append([aAtom[0]]+f_CorrectFloat(self.__ConvertCoordinate__(nIn, Lattice.coord_Conventional, aAtom[1:4])))
#These unit should not be supported here
#       elif (unit == "crystal" or unit == "crys"):
#           listAtom2 = listAtom1
        
        return listAtom2            
            
                          
    def MoveCell(self,vDistance=[0,0,0] ):
        '''
        Move Cell by a specific vector

        :param vDistance: The vector to move, unit of crystal coordinate of conventional vector
        :param nAtomIndex: which atom to move, if not set vector will be used, vector will be neglected otherwise
        :param nNewPos: Which position the atom will be moved to
        '''
        # Convert conventional crystal coord to primitive crystal coord
        #vMove = self.__GetPrimitiveCoordFromNormalCoord__(vDistance)
        vMove = self.__ConvertCoordinate__(Lattice.coord_Conventional, Lattice.coord_Primitive, vDistance)
        
        #reversely move  atom = move cell
        for aAtom in self.listAtom:
            aAtom[1] -= vMove[1]
            aAtom[2] -= vMove[2]
            aAtom[3] -= vMove[3]
            
            
    def ExpandCell(self,stRefAxis="001",nCell=1,nVacuum=1):
        '''
        Then extended Cell x times plus vacuum y times of original cell.
        Some high symmertry may broke for this reason. 

        :param stRefAxis: the axis to extended ( always postive ), only 100,010 and 001 is supported
        :param nCell: the number of cell layer
        :param nVacuum: the length of vacuum layer ( relative to cell length )
        '''
        nAxis = stRefAxis.find("1")
        #Get All atom and extended them
        #Duplicate all atom and store in cartesian coordinate
        listAtom2 = []
        for i in range(0,nCell):
            for aAtom in self.listAtomInConvCell:
                aAtom2 = copy.deepcopy(aAtom)
                aAtom2[nAxis+1] += i
                listAtom2.append([aAtom2[0]]+self.__ConvertCoordinate__(Lattice.coord_Conventional, Lattice.coord_Cartesian,aAtom2[1:4]))   
        # Deal with symmetry
        self.fLatticeLength[nAxis] *= (nCell+nVacuum)
        if ( self.nLatticeType == 1 or self.nLatticeType == 2 or self.nLatticeType == 3): #Cubic
            self.nLatticeType = 6
        elif ( self.nLatticeType == 6 or self.nLatticeType == 7): #
            if ( nAxis != 2):
                self.nLatticeType = 8
        #Clear all atoms and store new atoms
        self.listAtom = []
        self.AddAtomFromCartesian(listAtom2)
        
        self.CorrectError()
        
        
        
            
    def RotateCell(self,arRefPoint=[0,0,0],stRefAxis="001",stNewAxis="001"):
        '''
        Transform a cell to another orientation ( for surface SLAB )
        Fix ref point, rotate ref axis (a or b or c) to new orientation, xyz axis will be rotated as well
        Equally = all atom rotate reverse
        
        Note: the rotation axis is vertical to the plane of old ref axis and new axis
        All action is do in conventional cell. 
        '''
        arCell =self.ConventionalCellVector
        mCellT = lu.f_Matrix_transpose(arCell)
        #get new axis from miller index
        arOldAxis = [int(x) for x in stRefAxis]
        vOldAxis = lu.f_Matrix_transpose(lu.f_Matrix_dot(mCellT, lu.f_Matrix_transpose([arOldAxis])))[0]   
        arNewAxis = [int(x) for x in stNewAxis]
        vNewAxis = lu.f_Matrix_transpose(lu.f_Matrix_dot(mCellT, lu.f_Matrix_transpose([arNewAxis])))[0]
        #Get rotation matrix ( reversed, as it will affect atom )
        vRot =  lu.f_List_cross3(vOldAxis, vNewAxis)
        vRot = vRot / lu.f_List_norm3(vRot) # shoulde be unit vector
        nAngle = -lu.f_List_angle3(vOldAxis,vNewAxis)
        mRot = lu.f_Matrix_zeros(3,3)
        cosA = math.cos(nAngle)
        sinA = math.sin(nAngle)
        mRot[0][0] = cosA + ( 1-cosA)*vRot[0]*vRot[0]
        mRot[0][1] = (1-cosA)*vRot[0]*vRot[1]-sinA*vRot[2]
        mRot[0][2] = (1-cosA)*vRot[0]*vRot[2]+sinA*vRot[1]
        mRot[1][0] = (1-cosA)*vRot[1]*vRot[2]+sinA*vRot[0]
        mRot[1][1] = cosA + (1-cosA)*vRot[1]*vRot[1]
        mRot[1][2] = (1-cosA)*vRot[1]*vRot[2]-sinA*vRot[0]
        mRot[2][0] = (1-cosA)*vRot[2]*vRot[0]-sinA*vRot[1]
        mRot[2][1] = (1-cosA)*vRot[0]*vRot[1]+sinA*vRot[2]
        mRot[2][2] =  cosA + (1-cosA)*vRot[2]*vRot[2]

        #verify if rotation is right ( for old axis )
        vOldCalc = lu.f_Matrix_transpose(lu.f_Matrix_dot(mRot,lu.f_Matrix_transpose([vNewAxis])))
        print("Result verify:")
        print(vOldCalc)
        print(vOldAxis)
        #Do rotate
        self.nLatticeType = 14
        listAtom = self.listAtomInConvCell
        self.listAtom = [] # clear atom list
        for aAtom in listAtom:
            aAtom2 = copy.deepcopy(aAtom)
            aAtom2[1:4] = lu.f_Matrix_transpose(lu.f_Matrix_dot(mRot, lu.f_Matrix_transpose([aAtom2[1:4]])))[0]
            self.AddAtom([aAtom2], "crystal", "conv")
        
        self.CorrectError()
    
    def ChangeCell(self,stMode="c:1"):
        '''
        Change a specific Lattice object with specific method
        Note: any change on crystal axis will not affect the internal coordinate of atom
        In clayer mode, if there exist a distance much larger than others, it will be treated as Slab vacuum and not extended

        :params stMode: the way to change lattice, "c" and "tmdc" is supported. in format "mode:para;mode:para"
        '''
        listMode = f_Split_EnumerateRangeString(stMode)
        aMode = listMode[0]# only first mode is useds
        for aNamePara in aMode:
            stName = aNamePara[0]
            aPara = aNamePara[1]
            if ( stName == "c" ): #modify c axis
                aPara = float(aPara)
                self.fLatticeLength[2] *= aPara
            elif ( stName == "tmdc"): # modify c axis while keep M-F distance in c axis
                aPara = float(aPara)
                listAtom = self.listAtomInConvCell
                self.listAtom = []
                self.fLatticeLength[2] *= aPara
                for aAtom in listAtom: # test whether in 0.25 part or 0.75 part, and keep distance
                    if ( abs(aAtom[3]- 0.75) < abs(aAtom[3]-0.25)  ):
                        nRef = 0.75
                    else:
                        nRef = 0.25
                        
                    aAtom[3] = (aAtom[3]-nRef)/aPara+nRef
                    self.AddAtom([aAtom], "internal", "conv")
            elif (stName == "clayer"): #along c axis, modify distance between layers while keep all non-vdw bond ( detected by bond length ) not change.
                aPara = float(aPara)
                listAtom = self.GetAtomList("bohr", "conv")

                #self.listAtom = []
                #self.fLatticeLength[2] *= aPara
                # detect all layer ref point 
                # note the layer near c=0 is not moved at all
                #generate atom layer    
                listLayer  = self.GetLayerAtom()
                listDistance = []
                for i in range(0,len(listLayer)-1):
                    listDistance.append(listLayer[i+1][0][0][3]-listLayer[i][-1][-1][3])
                #last -> first
                listDistance.append(listLayer[0][0][0][3]+self.fLatticeLength[2]-listLayer[-1][-1][-1][3])
                
                listCAdd = [0]
                for i in range(0,len(listDistance)):
                    listCAdd.append(listCAdd[-1]+listDistance[i])
                    
                #Create new structure
                self.listAtom = []
                self.fLatticeLength[2] += listCAdd[-1]*(aPara-1)
                #ignore slab vacuum 
                if ( len(listDistance) > 1 ):
                    if ( listDistance[-1] / listDistance[-2] > 1.5):
                        print("Slab detected, the vacuum layer distance won't change.")
                        self.fLatticeLength[2] += (listCAdd[-2]-listCAdd[-1])*(aPara-1)
                
                #print(listLayer)
                #print(listCAdd)
                for i,aLayer in enumerate(listLayer):
                    for listAtomLayer in aLayer:
                        for aAtom2 in listAtomLayer:
                            aAtom = copy.deepcopy(aAtom2)
                            aAtom[3] += listCAdd[i] * (aPara-1)
                            self.AddAtom([aAtom], "bohr", "conv")
            else:
                raise ValueError("Unsupported method to change cell!")

    def GetLayerAtom(self,dThreshold = -1.0):
        '''
        Get atom layer along c axis

        :param dThreshold: The minimum distance between 2 layer. Auto detect if < 0. For example, all atom in range of 4 bohr is treated as same atom layer ( van der Waals ). For a ionic crystal, it should less than 0.1. 
        '''
        #detect whether the slab is suggested by arrange all atom along c axis
        #sort atom by c axis
        listTmp = [(x[3],i,x) for i,x in enumerate(self.GetAtomList("bohr", "conv"))]
        listTmp.sort()
        listAtom = [x[2] for x in listTmp]
        listAtomLayer = []
        #print(listAtom)
        
        #print("there are %d / %d atom in a cell" % (len(listAtom),len(result.listAtom)))
        #find can distance be arraged
        listDistance = [] # the distance between different atom layer 
        for i in range(0,len(listAtom)-1):
            listDistance.append(listAtom[i+1][3]-listAtom[i][3])
        listDistance.sort()
        if ( len(listDistance) < 1):
            raise ValueError("Only 1 atom in cell, seems not a suitable cell.")
        dThresholdGuess = listDistance[0] * 0.9
        #print(listDistance)
        for i in range(0,len(listDistance)-1):
            if ( listDistance[i]*1.5 < listDistance[i+1] and listDistance[i+1] > Lattice.dErrMax): # too large difference between atom distance, means there may be physical split
                dThresholdGuess = listDistance[i+1] * 0.99
                #break
        
        #print("Layer distance: ", listDistance)
        print("Detected layer split distance: %f" % dThresholdGuess)
        if ( dThreshold < 0):# use detected
            dThreshold = dThresholdGuess
            
        #generate atom layer    
        for aAtom in listAtom:
            if ( len(listAtomLayer) == 0):
                listAtomLayer.append([aAtom])
            else:
                if ( aAtom[3]-listAtomLayer[-1][-1][3] > 0.1):
                    listAtomLayer.append([aAtom])
                else:
                    listAtomLayer[-1].append(aAtom)

        #generate layer by threshold
        listLayer = []    
        for aAtomLayer in listAtomLayer:
            if ( len(listLayer) == 0):
                listLayer.append([aAtomLayer])
            else:
                if ( aAtomLayer[0][3]-listLayer[-1][-1][-1][3] > dThreshold):
                    #if ( len(listLayer[-1]) != nUnitLayer   ):
                        #print("Specific unit layer is %n atom layer but %n is found by threshold." % (nUnitLayer,len(listLayer[-1])))
                    listLayer.append([aAtomLayer])
                else:
                    listLayer[-1].append(aAtomLayer)
                      
        return listLayer              
                
    def CreateSlab(self,dThreshold = -1.0, nLatt = 7 , nVacuum = 5, nFix = 3, nShift = 0):
        '''
        Slice a n-layer slab from bulk cell with v layer vacuum along c axis. The atoms will be placed in the middle of the cell.
        A layer ( Unit layer ) is a combination of continuious atom layer. An atom layer is all atom in same c internal coordinate.
        Unit layer doesn't require to be same between each two.
        For example, MoS2 is a layer structure, and a layer contains S-Mo-S, 3 atom layers
        In AB stack, A is both an atom layer and an unit layer. B is another.
        repeat 5 times will create a A B A B A slab. Shift = 1 create B A B A B slab.

        :param nUnitLayer: the mono atom layer count in one unit layer in this manipulation. ( not used currently)
        :param nLatt: the layer counts of unit layer.
        :param nVacuum: the layer of vacuum.
        :param nFix: the layer of atom never move ; the fixed indexs will be returned as return value
        :param nShift: the start layer, counting from  
        :param dThreshold: The minimum distance between 2 layer. Auto detect if < 0. For example, all atom in range of 4 bohr is treated as same atom layer ( van der Waals ). For a ionic crystal, it should less than 0.1. 
        '''
        result = copy.deepcopy(self)
         
        #get atom in in layer listLayer->AtomLayer -> Atom
        listLayer = self.GetLayerAtom(dThreshold)
        
        print("Program detected layer:")
        for i,aLayer in enumerate(listLayer):
            stInfo = "Layer %d: " % (i+1)
            for aAtomLayer in aLayer:
                for aAtom in aAtomLayer:
                    stInfo += aAtom[0] + " "
            print(stInfo)
        

        
        #Calculate the maxium supercell along c axis contained in the slab; rounded celing
        nCmax = (nLatt + len(listLayer) -1) / len(listLayer) - 1
        #record cell vector to move atom
        listVector = result.ConventionalCellVector
        #print(listVector)
        #Expand to minium supercell with 1 more layer vacuum to add residual atoms
        result.ExpandCell("001",nCmax,nVacuum + 1)
        #print(result.fLatticeLength)
        
        #Those step is useless for atom should be sorted in slab
        result.listAtom = [] 
        
        
        #Add Atoms
        i0 = 0
        nAtomCount = 0
        listFix = []
        while ( i0 < nLatt):
            i = i0 + nShift # layer number is relative to i
            i2 = i % len(listLayer)
            nC = i / len(listLayer)
            bFix = ( i0 >= (nLatt-nFix)/2  and i0 < (nLatt+nFix)/2  ) # fix count is relative to i0
            for aAtomLayer in listLayer[i2]:
                for aAtom in aAtomLayer:
                    aAtom2 = copy.deepcopy(aAtom)
                    for j in range(0,3):
                        aAtom2[j+1] += listVector[2][j] * nC
                    result.AddAtom([aAtom2], "bohr", "conv")
                    if ( bFix):
                        listFix.append(nAtomCount)
                    nAtomCount += 1
            i0 += 1
            
        #print(result.listAtom)
        
        #Move atoms to centre of c axis
        #Due to a slab model can never be a complex lattice, primitive cell and c axis can be directly used
        #Only do when Vacuum layer > 0
        if ( nVacuum > 0):
            listC = [x[3] for x in result.listAtom]
            fMove = 0.5 - (max(listC)+min(listC))/2
            for aAtom in result.listAtom:
                aAtom[3] += fMove
        
        return result,listFix
            
    def CreateSurfCell(self,stRefAxis="001",dRot=0,mRot=[[1,0],[0,1]]):
        '''
        Slice a cell that have a plance perpendicular to specifc stRefAxis ( for surface model )
        The input miller index should be simlified, like 002 may cause unexpected result
        '''
        
        vRef = [int(stRefAxis[0]),int(stRefAxis[1]),int(stRefAxis[2])]
            
        #Alogrithm error
        #Get new axis
        #listAns = f_SolveCross(vRef[0],vRef[1],vRef[2])
        #if there is [1,x,0] and [y,1,0], keep it if it is among the lowest
        #nLowest = listAns[0][0]
        #ans = None
        #for aAns in listAns:
        #    if ( aAns[0] == nLowest ):
        #        if ( aAns[1][0] == 1 and aAns[1][4] == 1 ):
        #            ans = aAns[1]
        #    else:
        #        break
        #
        #if ( ans == None):
        #    ans = listAns[0][1]
        
        AxisInt = f_GetSurfaceFromMiller(vRef)
        print(AxisInt)
        
        #AxisInt = [[],[],[]] # new axis ,bases are old axis
        #AxisInt[0] = ans[0:3]
        #AxisInt[1] = ans[3:6]
        #AxisInt[2] = copy.deepcopy(vRef)

        #Rotate
        #dRot = dRot * math.pi / 180
        #mRot = f_RotateMatrix(dRot,AxisInt[2])
        #AxisInt[0] = lu.f_Matrix_transpose(lu.f_Matrix_dot(mRot,lu.f_Matrix_transpose([AxisInt[0]])))[0]
        #AxisInt[1] = lu.f_Matrix_transpose(lu.f_Matrix_dot(mRot,lu.f_Matrix_transpose([AxisInt[1]])))[0]
        AxisIntNew0 = lu.f_List_Op_List(lu.f_List_Op_Scalar(AxisInt[0],"*",mRot[0][0]),"+", lu.f_List_Op_Scalar(AxisInt[1],"*",mRot[0][1]))
        AxisIntNew1 = lu.f_List_Op_List(lu.f_List_Op_Scalar(AxisInt[0],"*",mRot[1][0]),"+", lu.f_List_Op_Scalar(AxisInt[1],"*",mRot[1][1]))
        AxisInt[0] = AxisIntNew0
        AxisInt[1] = AxisIntNew1
        #print(AxisInt)

        
        #Create cartesian axis
        OriginalVector = self.ConventionalCellVector
        #OriginalVector = self.PrimitiveCellVector
        listVector = []
        for aAxis in AxisInt:
            Vector = [0,0,0]
            for i in range(0,3):
                for j in range(0,3):
                    Vector[j] += aAxis[i] * OriginalVector[i][j]
            listVector.append(Vector)
        
        #print(OriginalVector)
        #print(listVector)
        
        #Get possible cell enlarge times to slice a supercell from it
        #Detect 8 vertex of the supercell
        #nMax = 1
        listVertex = []
        for a1 in range(0,2):
            for a2 in range(0,2):
                for a3 in range(0,2):
                    list_a = [[a1],[a2],[a3]]
                    listVertex.append(lu.f_Matrix_transpose(lu.f_Matrix_dot(lu.f_Matrix_transpose(AxisInt),list_a))[0])
        print(listVertex)
        arMax = [0,0,0,0,0,0]
        for aVertex in listVertex:
            for i in range(0,3):
                if ( aVertex[i] < arMax[i*2]):
                    arMax[i*2] = aVertex[i]
                elif ( aVertex[i] > arMax[i*2+1]):
                    arMax[i*2+1] = aVertex[i]
                    
        arMax = [int(x) for x in arMax]

        #print(listVertex)
        print(arMax)
        
        listAtom = self.GetAtomList("bohr","conv",arMax)
        print(listAtom)
        result = copy.deepcopy(self)
        result.listAtom = []
        result.ReadFromPrimitiveCellVector(listVector, 14)
        #convert all lattice under non-standard vectors to primitive coords ( always standards )
        for aAtom in listAtom:
            result.listAtom.append([aAtom[0]] + lu.f_Matrix_transpose(lu.f_Matrix_dot( lu.f_List_inv3(lu.f_Matrix_transpose(listVector)) , lu.f_Matrix_transpose([aAtom[1:4]])))[0])
            #print(result.listAtom[-1])
            result.CheckAtom()
        #self.AddAtom([ [x[0]] + lu.f_Matrix_transpose(lu.f_Matrix_dot(mRot,lu.f_Matrix_transpose([x[1:4]])))[0] for x in listAtom] , "bohr", "conv")
              
        #print(result.listAtom)
        return result
            
    def CheckAtom(self):
        '''
        Check whether an Atom is < 0 or > 1, or there are duplicated atoms
        '''
        #remove precision problem
        self.CorrectError()
        # check 0 < x < 1
        for aAtom in self.listAtom:
            for i in range(1,4):
                while ( aAtom[i] >= 1 ):
                    aAtom[i] -= 1
                while ( aAtom[i] < 0 ):
                    aAtom[i] += 1
        #check duplicate
        i = 0
        while ( i < len(self.listAtom)):
            j = i + 1
            while ( j < len(self.listAtom)):
                aAtom = self.listAtom[i]
                aAtomRef = self.listAtom[j]
                bSame = True
                for k in range(1,4):
                    if ( abs(aAtomRef[k] - aAtom[k]) > Lattice.dErrMax):
                        bSame = False
                        break
                if ( bSame):
                    #print("%d atom is duplicated with %d" % (j,i))
                    #print(aAtom,aAtomRef)
                    del self.listAtom[j]
                else:
                    j += 1
            i += 1       
                
def f_GetReciprocalLattice(matR):
    '''
    Get reciprocal lattice vectors from real lattice  or vice versa
    Recprocal lattice is just inversion of real lattice * 2\pi, indeed.
    '''
    matK =  lu.f_Matrix_Op_Scalar(lu.f_Matrix_transpose(lu.f_List_inv3(matR)),"*",2*math.pi)
    return matK
    


class KPointsT():
    '''
    Represent a set of kpoints
    Contains k-point coordinate, k-point name
    '''

#format enum
    format_xyz = 1
    format_xyzn = 2
    format_xyzw = 3
    format_xyzwn = 4
    format_nxyz = 5
    format_nxyzw = 6

#Dictionary to create k-points list from specified number of high-symmetryk-points
#All lines between k-points are connected, but for 4,6,8, one duplicate path is choosed as it is impossible to draw all lines in one stroke !
    dicOneLine = {
            2:(1,2),
            3:(1,2,3,1),
            4:(1,2,3,4,1,3,2,4),
            5:(1,2,3,4,5,1,3,5,2,4,1),
            6:(1,2,3,4,5,6,1,3,5,2,4,6),
            7:(1,2,3,4,5,6,7,1,3,5,7,2,4,6,1,4,7,3,6,2,5,1),
            8:(1,2,3,4,5,6,7,8,1,3,5,7,1,4,6,8,2,4,7,2,5,8,3,6,1,5,8,4,7,3,2,6)
            }

    def __init__(self,stFile=None):
        self.listKPt= [] #The k-point list. Format: name,x,y,z,weight ( x,y,z maybe crystal, cartesian, tpiba )
        self.stMode = "crystal" # The k-point mode. Possible value include cartesian(cart), crystal, tpiba, tpibabc. Note, cartesian coordination are in unit bohr or bohr^{-1}. "ang" and "bohr" means Cartesian, but in given unit
        self.latt = None #The structure corresponding to this k-point list,used to convert in 2pi/a unit and crystal unit
        self.kFormat = KPointsT.format_xyz #default k-points format to read
        if ( stFile != None ):
            self.ReadFromFile(stFile)

    def __ReFormat__(self,kpt):
        '''
        make kpt format as namekx,ky,kz,weight

        Acceptable k-point format ( "," means space or tab )
           [name],kx,ky,kz,[weight]
        or kx,ky,kz,[name]
        or kx,ky,kz,weight,[name]
        
        All element can be either string or float
        '''
        if ( len(kpt) < 3 ):
            raise ValueError("Unexpected data format in k-point file")
        elif ( len(kpt) == 3 ): #must be 3 digit
            return [""] + [float(x) for x in kpt] + [1.0]
        else:
            # detect does any value start with character
            listNum = []
            listString = [] 
            for i in range(0,len(kpt)):
                if ( isinstance(kpt[i],str) and not kpt[i][0].isdigit() and not kpt[i][0] == "-"):
                    listString.append(kpt[i])
                else:
                    listNum.append(float(kpt[i]))

            listString = listString + listNum

            if ( len(listString) == 4 ): 
                if ( len(listNum) == 3):
                    # if there is not weight, add it
                    listString.append(1.0)
                else:
                    # there is no name, add it
                    listString = [""] + listString


            return listString

    def CreateFromPointLine(self,listK2,listPointsCount=None,nTotal=None,mLatt=None,nDivisor=None):
        '''
        Read k-point from x,y,z + point along the line format
        kpt must be in formar [x,y,z], [x,y,z,name] or [name,x,y,z]
        The mode is not changed during this process

        :param listK2: the list of k-points in KPointsT readable format
        :param listPointsCount: the number of points along each line
        :param nTotal: the total number of points. The lattice must be specified if use this option, and k-points will be spread homogeneously along lines
        :param mLatt: the lattice vectors for k-points, necessary if k-points is in internal coordinate
        :param nDivisor: if this number is set, then number of points along each line must be a divisor/factor of this number. This function is especially useful for WIEN2k. Set to 0 means the program will automatically choose a number.
        :return: the divisor
        '''
        def FindNear(x,ar):
            '''
            Find nearest number of x in list y ( y in ascending order )
            '''
            for i,y in enumerate(ar):
                if ( x < y):
                    if ( i == 0):
                        return y
                    else:
                        if ( x*x < y*ar[i-1]):
                            return ar[i-1]
                        else:
                            return y
            return ar[-1] #return final one if not found

        self.listKPt = []
        listK = [self.__ReFormat__(x) for x in listK2]

#Use or create number of points along each lines
        if ( listPointsCount == None):
            if ( nTotal == None):
                raise ValueError("Must specify number of k-points along each lines or the total number of k-points")
            if ( self.stMode != "cart" and self.stMode != "cartesian" and mLatt  == None):
                raise ValueError("Must specify lattice vectors when create k-points list from internal coordinate")
#Create ca
            stModeOld = self.stMode
#Use self as temp
            self.listKPt = copy.deepcopy(listK)
            self.ConvertUnit("cart",mLatt) 
#Calculate distance
            listDist = []
            for i in range(0,len(listK2)-1):
                listDist.append(lu.f_List_norm3(lu.f_List_Op_List(self.listKPt[i][1:4],"-",self.listKPt[i+1][1:4])))
#Calculate number of k-points
            dDist = sum(listDist)
            listPointsCount = [ int(round( (nTotal-1) * x / dDist)) for x in listDist]
#Restore self
            self.stMode = stModeOld

#Modify number of k-points to fit divisor. Note the difference between k-points coordinate must be considered too
        if ( nDivisor != None):
            nDivisorOld = nDivisor
            if ( nDivisor == 0):
                nDivisor = 2**4*3**4*5**2 # The guess value for most possible 
            listPointsCount2 = []
            for i,nCount in enumerate(listPointsCount):
#Calculate fractional coordinate in output unit
                vDiff = lu.f_List_Op_List(listK[i][1:4],"-",listK[i+1][1:4])
#Get n/m fractional representaion
                vDiff = [ f_FindFraction(x) for x in vDiff ]

#Check the divisor for each dimension
                for n1,n2 in vDiff:
                    if ( nDivisor % n2 != 0):
                        n3 = nDivisor
                        nDivisor = f_FindLCM(nDivisor,n2)
                        print("%i is not divided by k-points coordinate, increased to %i" % ( n3,nDivisor))
#Find the GCD for divisor for 3-D fractional
                n3 = f_FindGCD(nDivisor/vDiff[0][1],f_FindGCD(nDivisor/vDiff[1][1],nDivisor/vDiff[2][1]))

#Calculate all divisor
                listDivisor = filter(lambda x: n3/x*x==n3,range(1,n3+1))
                listPointsCount2.append(FindNear(listPointsCount[i],listDivisor))
            #Search for one line end

            print("Number of k-points is modified to %i according to idv" % (sum(listPointsCount2)+len(listPointsCount2)+1) )
            print("Equal distance distribution:")
            print(listPointsCount)
            print("Current distribution:")
            print(listPointsCount2)
            listPointsCount = listPointsCount2
            

        self.listKPt = [] 
        if ( len(listK) - len(listPointsCount) != 1):
            print("Error: Number of k-points must be number of lines + 1")
            return

        for i in range(0,len(listPointsCount)):
            l = listPointsCount[i]
            self.listKPt.append(listK[i])
            for j in range(1,l):
                kpt = lu.f_List_Op_List(lu.f_List_Op_Scalar(listK[i][1:4],"*",1.0*(l-j)/l),"+",lu.f_List_Op_Scalar(listK[i+1][1:4],"*",1.0*j/l))
                self.listKPt.append(self.__ReFormat__(kpt))

        #add the last value
        self.listKPt.append(listK[-1])
        return nDivisor

    def CreateKPointsFromSymmetry(self,nGroup):
        '''
        Create high-symmetry k-points list from spacegroup symmetry
        The lattice must be set before invoke this function

        :param nGroup: the 230 space group index (1-230)
        '''
        def ParseDiv(st):
            '''
            Return numerical for "n/m"
            '''
            list1 = st.split('/')
            if ( len(list1) == 2):
                return float(list1[0])/float(list1[1])
            else:
                return float(st)

        if ( self.latt == None):
            raise ValueError("Lattice must be specified before create k-points in band structure")
        listSpec = None
#Detect UAC : whether this space group have A / C axis variants depending on the angles
        bHasUAC = kpt_spec.listSpaceGroupKPointsUAC[nGroup-1] != None
        nUAC = 0
#Get k-points set
        if ( bHasUAC):
            if ( abs(self.latt.fLatticeAngle[1]-90) > 0.01):
                nUAC += 2
            if ( abs(self.latt.fLatticeAngle[2]-90) > 0.01):
                nUAC += 3
            if ( nUAC != 2 and nUAC != 3):
                raise ValueError("The cell shape is not consistent with space group")
#Unique axis C
            if ( nUAC == 3):
                listSpec = kpt_spec.listSpaceGroupKPointsUAC[nGroup-1]
        
        if ( listSpec == None):
#For non-UAC cases, it is possible to find multiple k-points definition due to different a,b,c relationship
            (a,b,c) = self.latt.fLatticeLength
            for kpt_cond in kpt_spec.listSpaceGroupKPoints[nGroup-1]:
                if ( kpt_cond[0](a,b,c)):
                    listSpec = kpt_cond[1]
                    break
            if ( listSpec == None):
                raise ValueError("The cell shape is not consistent with space group")

        #Parse k-points data
        #Always use primitive cell
        listK = [ [x[0]]+[ParseDiv(y) for y in x[1]] + [1.0] for x in listSpec]
#Try to add line for each two k-points
        listK2 = [ listK[x-1] for x in KPointsT.dicOneLine[len(listK)]]

        #print(listK2)
        return listK2
        #self.CreateFromPointLine(listK2,nTotal=nTotal,mLatt=self.latt.PrimitiveCellVector)
     
    def ReadFromList(self,list_kp):
        '''
        Create object from simple k-point list in format 
        [ [x1,y1,z1] ...]
        or 
        [ [x1,y1,z1,w1] ...]
        or
        [ [name,x1,y1,z1,w1] ...]
        The unit is not concerned in this case
        '''
        n_item = len(list_kp[0])
        if (n_item == 3):
            self.listKPt = [ [""] + x + [1.0] for x in list_kp]
        elif (n_item == 4):
            self.listKPt = [ [""] + x for x in list_kp]
        elif (n_item == 5):
            self.listKPt = list_kp
        else:
            raise TypeError("Unknown format")

    def ReadFromFile(self,stFileName):
        '''
        Read k-point from a file in raw format
        format as first line is k-point data mode; if it is not "crystal"/"cartesian"/"tpiba"/"piba", then it is set to crystal and this line is used as k-point
        each k-point format ( "," means space or tab )
        * [name],kx,ky,kz,[weight]
        * kx,ky,kz,[name]
        * kx,ky,kz,weight,[name]
        name must start with character and no space in it; if it start with number, then it is treated as kx!
        '''
        f = open(stFileName,'r')
        
        stLine = f.readline().strip()
        ar_stLine = f.readlines()
        f.close()

        if ( not stLine in ("crystal","tpiba","cartesian")):#Not a format data
            ar_stLine = [stLine] + ar_stLine
        else:
            self.stMode = stLine

        for stLine in ar_stLine:
            arLine  = stLine.strip().split()
            if ( len(arLine) == 0 ): #Empty line
                continue
            self.listKPt.append(self.__ReFormat__(arLine))
        pass

    def WriteToFile(self,stFileName,bWriteName=False,bWriteWeight=False):
        '''
        Write k-point list to a file with/without

        :param bUseName: write name or not, default no
        '''
        f = open(stFileName,'w')
        if ( bWriteName ):
            for kpt in self.listKPt:
                f.write(" ".join(["%10s" % kpt[0]]+["%12.7f" % x for x in kpt[1:4]]))
                f.write('\n')
        else:
            for kpt in self.listKPt:
                f.write(" ".join(["%12.7f" % x for x in kpt[1:4]]))
                f.write('\n')
   
    def GetTurningPoint(self):
        '''
        Return index of turning points in the k-point path
        Start point and end point are not included
        Index start from 0

        Warning: sometimes a band structure may contains duplicated points (like band-mode in VASP), which is bad so we skip such k-points
        '''
        list_turnpt = [] 
        for i in range(1,len(self.listKPt)-1):
#Check if two points are the same, if it is, we always mark it is
            if ( all([x==y for x,y in zip(self.listKPt[i-1][1:4], self.listKPt[i][1:4])]) or 
                    all([x==y for x,y in zip(self.listKPt[i][1:4], self.listKPt[i+1][1:4])])):
                list_turnpt.append(i)
                continue

            dAngle= lu.f_List_angle3(lu.f_List_Op_List(self.listKPt[i-1][1:4],"-",self.listKPt[i][1:4]),lu.f_List_Op_List(self.listKPt[i][1:4],"-",self.listKPt[i+1][1:4]))
            if ( dAngle > math.pi/3):
                list_turnpt.append(i)
        return list_turnpt


    def ConvertUnit(self,stModeNew,mLatt=None):
        '''
        Convert k-points coordinate between crystal,cartesian and tpiba unit based on structure
        Note this conversion mostly only appeared in orthogonal lattice, as in practice, non-orthogonal lattice never use tpiba unit.
        Equation : matLattVec . vecCrystal = vecCartesian

        :param stModeNew: the new k-point unit
        :param mLatt: the real-space lattice vector that k-point is based on; note for some program it is primitive lattice ( like VASP and YAehmop ), and others are conventional lattice (Wien2K cubic). If None is used, we use conventional lattice for tpiabc (which belongs to Wien2K cubic) and primitive for others. 
        :return: whether the conversion completed
        '''
        from constants import Ang2Bohr
        #Run even if two modes are same to verify its 
        #if ( stModeNew == self.stMode):
        #    return

        #Just * or / if ang -> bohr and bohr -> ang
#Note this is reciprocal lattice so reveresed
        dScale = 1
        if (stModeNew == "ang"):
            dScale *= Ang2Bohr
        if (self.stMode == "ang"):
            dScale /= Ang2Bohr

        b_auto = mLatt is None
        if (b_auto):

            if (self.latt is None):
                print("Unable to convert k-point unit as lattice information unavailable.")
                return False

            if (self.stMode == "tpibabc"):
                mLatt = self.latt.ConventionalCellVector
            else:
                mLatt = self.latt.PrimitiveCellVector
#tpiba unit
        if (self.stMode in ["tpiba","tpibabc"] or stModeNew in ["tpiba","tpibabc"]):
            mLattConv = self.latt.ConventionalCellVector
            listUnit = [ 2*math.pi/lu.f_List_norm3(x) for x in mLattConv ]
            dTpibaUnit = listUnit[0]
#       dTpibaUnit = 2*math.pi/lu.f_List_norm3(mLattConv[0])
#       listUnit = [dTpibaUnit,2*math.pi/lu.f_List_norm3(mLattConv[1]), 2*math.pi/lu.f_List_norm3(mLattConv[2])]
#       print(listUnit)

        if (self.stMode in ["crystal","cart","cartesian","bohr","ang"]  or stModeNew in ["crystal","cart","cartesian","bohr","ang"] ):
            matR = f_GetReciprocalLattice(lu.f_Matrix_transpose(mLatt))
            matRi = lu.f_List_inv3(matR)
#           print("matR",matR)
#           print("matRi",matRi)

        #Convert current k-point in any unit to cartesian
        mKPt = []
        for kpt in self.listKPt:
            mKPt.append(kpt[1:4])

        mKPt = lu.f_Matrix_transpose(mKPt)

        if ( self.stMode == "crystal"):
            mKPt = lu.f_Matrix_dot(matR,mKPt)
        elif ( self.stMode == "tpiba"):
            mKPt = lu.f_Matrix_Op_Scalar(mKPt,"*",dTpibaUnit)
        elif ( self.stMode == "tpibabc"):
            for i in range(0,len(mKPt[0])):
                for j in range(0,3):
                    mKPt[j][i] *= listUnit[j]

        #Convert cartesian to any unit

        #Recontruct lattice
        if (b_auto):
            if (stModeNew == "tpibabc"):
                mLatt = self.latt.ConventionalCellVector
            else:
                mLatt = self.latt.PrimitiveCellVector
            matR = f_GetReciprocalLattice(lu.f_Matrix_transpose(mLatt))
            matRi = lu.f_List_inv3(matR)
            dTpibaUnit = 2*math.pi/lu.f_List_norm3(mLatt[0])
            listUnit = [dTpibaUnit,2*math.pi/lu.f_List_norm3(mLatt[1]), 2*math.pi/lu.f_List_norm3(mLatt[2])]
#           print("matR",matR)
#           print("matRi",matRi)

        if ( stModeNew == "tpiba" ):
            mKPt = lu.f_Matrix_Op_Scalar(matR,"/",dTpibaUnit)
        elif ( stModeNew == "crystal"):
            mKPt = lu.f_Matrix_dot(matRi,mKPt)
        elif ( stModeNew == "tpibabc"):
            for i in range(0,len(mKPt[0])):
                for j in range(0,3):
                    mKPt[j][i] /= listUnit[j]


        mKPt = lu.f_Matrix_transpose(mKPt)

#Convert ang / bohr
        for i,kpt in enumerate(mKPt):
            self.listKPt[i][1:4] = [x*dScale for x in kpt]


        self.stMode = stModeNew
        return True

    def get_k_xaxis(self):
        '''
        Get a list of points on one axis for band structure plotting

        If two points are two far away, then they are treated as seperated, the distance between them is eliminated

        :return:list_x and list_name for x-axis and k-point names, list_break for the index of k-points where a new kpath starts
        '''
        list_dist2 = []
        for i, kpt in enumerate(self.listKPt):
            if (i == len(self.listKPt) -1):
                break
            list_dist2.append(lu.f_List_norm(lu.f_List_Op_List(kpt[1:4],"-",self.listKPt[i+1][1:4])))

        list_break = [0]
        for i in range(1, len(list_dist2)-1):
            s1 = list_dist2[i-1] + list_dist2[i+1]
            if (list_dist2[i] > s1 * 2 ): #Too long, treated as seperated
                list_dist2[i] = s1 * 2
                list_break.append(i+1)

        list_dist2 = [0] + list_dist2

        list_dist_set = []

        list_name = []
        list_dist = [] #x-axis in plot rely on k-k distance

        for i, dist in enumerate(list_dist2):
            if (i == 0):
                list_dist.append(0)
            else:
                list_dist.append(list_dist[-1] + dist)
            kpt = self.listKPt[i]
            if (kpt[0] != ""):
                name = kpt[0]
                if ( name[0] == '\\' ): #delete first backslash as neither gnuplot nor xmgrace use it. SIESTA use "\" for Latex. 
                    name = name[1:]
                list_name.append([i,name,list_dist[-1]])


#       for i,kpt in enumerate(self.listKPt):
#           if (i == 0):
#               list_dist.append(0.0)
#           else:
#               list_dist.append(list_dist[-1] + lu.f_List_norm(lu.f_List_Op_List(kpt[1:4],"-",self.listKPt[i-1][1:4])))
#           if (kpt[0] != ""):
#               name = kpt[0]
#               if ( name[0] == '\\' ): #delete first backslash as neither gnuplot nor xmgrace use it. SIESTA use "\" for Latex. 
#                   name = name[1:]
#               list_name.append([i,name,list_dist[-1]])

        return list_dist, list_name, list_break
      
