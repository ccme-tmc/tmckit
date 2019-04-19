#!/usr/bin/env python
import sys,os,shutil,re,math,subprocess,time
from chem_utils import *
from constants import *
from struct_utils import * 
from latt_utils import *
from list_utils import f_Matrix_Op_Scalar,f_List_Op_Scalar
import math

import copy
from common_caseutil import Lattice,f_GetLibDataPath,deprecated,KPointsT,f_env_RunMpirunCommand,NamelistIO,FormatError
from dos_utils import OrbitalNameT
from band_utils import BandsT
from py_xmlio import XMLSerilizer
import xml.etree.ElementTree as ET
from common_caseutil import f_Data_Write


mod_name = "qesp_utils"

# default pseudopotential tag for each element
qe_default_psp = "-pbe-tm"
kgrid_cutoff = 30.0   # the default k-mesh on i-direction is estimated by int(kgrid_cutoff/latt[i])+1
  
def f_QE_get(file,tag,unit=None):
  """ 
  f_QE_get: get the value of a paramter as represented by 'tag'
  """
  ifile= open(file,'r') 
  lines = ifile.readlines()
  nw = len(tag)
  for line in lines: 
    if line.strip()[0:nw] == tag:
      ltmp = line[nw:]

def f_qe_get_out(out_file,tag):
    pass  
  

def f_QE_Reset(inp,name,val):
  """
  Reset the value of a paramter in a PWscf input file  
  """ 
  sname = "f_QE_Reset"

  ofile = open(inp,'r')
  ofile.close()

def f_QE_Read_Struct_Out(file,iop_out=0):
  """ 
  Read struct from QE input or output file
    iop_out = 0 -- return lattice vectors and basis (internal coordinates) 
            = 1 -- return Cartesian coordinates in units of angstrom  
  """

  ifile= open(file,'r')

  found_none = False
    
  ibrav = -1
  natom = 0 
  lat_vec = None
  while (1):
    line = ifile.readline()
    if not line: 
      found_none = True
      break

    # ibrav
    if ibrav == -1 and "bravais-lattice index" in line:
      info = line.split()
      ibrav = int(info[3])    
      print "bravais-lattice index (ibrav) = %2d "%(ibrav)
  
    # alat 
    if "lattice parameter" in line:
      data = line.split()
      alat = float(data[4])
      print "lattice parameter (alat) = %12.6f "%(alat)
    
    # natom 
    if natom == 0 and "number of atoms/cell" in line:
      info = line.split()
      natom = int(info[4]) 
      print "the number of atoms/cell (natom) =%4d"%(natom) 

    if lat_vec == None and "crystal axes:" in line:
      print "Lattice Vectors for the initial structure:"
      lat_vec = []
      for i in range(3):
        data = ifile.readline().split()
        (x,y,z) = ( float(data[3]), float(data[4]), float(data[5]) )
        print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,x,y,z)
        lat_vec.append( [x,y,z] )
     
    if "Begin final coordinates" in line: break 

  if found_none: 
    print "ERROR: no info on the optimized structure has been found!"
    sys.exit(1)      

  print ""
  line = ifile.readline()
  if "new unit-cell volume" in line:
    print "Read output from a 'vc-relax' calculation"
    irelax = 1 
    data = line.split()
    vol_new = float(data[4]) 
    print "  new volume = %12.6f Bohr^3"%(vol_new)

    # get lattice vectors 
    line = ifile.readline()  # skip a blank line 

    line = ifile.readline()  # "CELL_PARAMETERS ..." 
    print "  lattice vectors in terms of the lattice parameter(alat=%12.6f):"%(alat) 
    lat_vec_old = lat_vec[:] 
    lat_vec = [] 
    for i in range(3):
      data = ifile.readline().split()
      (x,y,z) = (float(data[0]), float(data[1]), float(data[2]) )
      lat_vec.append([x,y,z]) 
      print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,x,y,z)

  else:
    irelax = 0 
    print "Read output from a 'relax' calculation"
  
  line = ifile.readline()  # skip the blank line 
  line = ifile.readline().strip()  #"ATOMIC_POSITIONS (...)"

  if   "alat"     in line:  mol_type = 'alat'
  elif "bohr"     in line:  mol_type = 'bohr' 
  elif "angstrom" in line:  mol_type = 'angstrom'
  elif "crystal"  in line:  mol_type = 'crystal' 
  mol = []   
  print line 
  for ia in range(natom):
    line = ifile.readline()
    data = line.split()
    atom = data[0]
    (x,y,z) = (float(data[1]), float(data[2]), float(data[3])) 
    print "%-6s %12.6f %12.6f %12.6f"%(atom,x,y,z)
    mol.append([atom,[x,y,z]])

  if mol_type == "alat" and irelax ==1 :
    # get the ratio between new and old lattice parameters
    alen_old = 0.0
    alen = 0.0
    for i in range(3):
      alen_old += lat_vec_old[0][i]**2
      alen += lat_vec[0][i]**2

    alat_old = alat
    alat = alat_old * sqrt(alen/alen_old)

    print "\n  new lattice vectors in terms of new lattice parameter (alat=%12.6f):"%(alat)
    for i in range(3):
      for j in range(3): lat_vec[i][j] *= alat_old/alat

      print "a(%1d)=(%12.6f,%12.6f,%12.6f)"%(i+1,lat_vec[i][0],lat_vec[i][1],lat_vec[i][2])

    print "\natomic positions in terms of new lattice parameter:"
    for ia in range(natom):
      for i in range(3):
         mol[ia][1][i] *= alat_old/alat
      atom = mol[ia][0]
      (x,y,z) = (mol[ia][1][0],mol[ia][1][1],mol[ia][1][2])
      print "%-6s %12.6f %12.6f %12.6f"%(atom,x,y,z)
      
  ifile.close()

  for iv in range(3):
    for i in range(3):
      lat_vec[iv][i] *= alat

  if mol_type == 'alat':
    for ia in range(natom):
      for i in range(3):
         mol[ia][1][i] *= alat

    basis = f_Latt_Cart_to_Crys(mol,lat_vec) 

  # convert units of lat_vec to Angstrom 
  for iv in range(3):
    for i in range(3):
      lat_vec[iv][i] *= Bohr2Ang

  if iop_out == 0: 
    return basis,lat_vec,ibrav 
  else:
    return mol
  

def f_QE_Write_Struct(prefix,mol,latt=None,ibrav=0):
  print """ calling f_Write_Struct_qesp: 
  Using the input molecular (and lattice constants and type ) to generate 
  an input file for PWscf with a few default options.
  "mol" delivers the molecular coordinates in units of angstrom
  !!!  the default options is most likely not appropriate for your systems !!! 
  """
  sname = "f_QE_Write_Struct"

  print "Create QE input for given structure"

  name = prefix.strip()+"-pw.in"
  ofile = open(name,'w')
  print "Create new PWscf inpup file: " + name 

  # check how many types of atoms are present in the molecule 
  nat = len(mol)
  ntyp,species,sp_index = f_Check_Species(mol) 

  outdir = './out-'+prefix.strip()
  wfcdir = os.getenv("ESPRESSO_TMPDIR",outdir) 
  pspdir = os.getenv("ESPRESSO_PSEUDO",'./') 
 
  # default input for control block 
  ofile.write("""&control
  calculation = 'scf'
  restart_mode = 'from_scratch'
  tstress = .true.
  tprnfor = .true.
""")
  ofile.write("  outdir = '"+ outdir+"'\n") 
  ofile.write("  wfcdir = '"+ wfcdir+"'\n") 
  ofile.write("  pseudo_dir = '"+ pspdir+"'\n") 
  ofile.write("  prefix = '"+ prefix+"'\n") 
  ofile.write("/\n\n")

  ofile.write("\n&system\n")
  ofile.write("  nat = %d, ntyp = %d\n"%(nat,ntyp)) 

  xyz_shift = [ 0.0, 0.0, 0.0 ]  

  if latt is None : # finite systems 

    sc = f_Struct_Size(mol)
    ofile.write("  ibrav = %d\n"%(8))
    ofile.write("  celldm(1) = %12.6f\n"%(sc[0]*Ang2Bohr))
    ofile.write("  celldm(2) = %12.6f\n"%(sc[1]/sc[0]))
    ofile.write("  celldm(3) = %12.6f\n"%(sc[2]/sc[0]))

  else:           # periodic systems 
    print "  Lattice Constants:"
    print "\t    a =%12.6f,    b =%12.6f,     c=%12.6f"%(latt[0],latt[1],latt[2])
    print "\talpha =%12.6f, beta =%12.6f, gamma=%12.6f"%(latt[3],latt[4],latt[5])
    a = latt[0]
    boa=latt[1]/latt[0]
    coa=latt[2]/latt[0]
    cos_bc= cos(latt[3]*Deg2Rad)
    cos_ac= cos(latt[4]*Deg2Rad) 
    cos_ab= cos(latt[5]*Deg2Rad) 
    
    ofile.write("  ibrav = %d\n"%(ibrav))
    ofile.write(  "  celldm(1) = %12.6f\n"%(a*Ang2Bohr))
    if ibrav == 4 or ibrav == 6 or ibrav == 7:  # hex or trigonal P
      ofile.write("  celldm(3) = %12.6f\n"%(coa))
    elif ibrav == 5:  # triangle R
      ofile.write("  celldm(4) = %12.6f\n"%(cos_bc))
    else:
      ofile.write("  celldm(2) = %12.6f\n  celldm(3) = %12.6f\n"%(boa,coa))
      if ibrav == 12 or ibrav == 13:  # Monoclinic
        ofile.write("  celldm(4) = %12.6f\n"%(cos_ab))
      if ibrav == 14: 
        ofile.write("  celldm(4) = %12.6f\n  celldm(5) = %12.6f\n  celldm(6) = %12.6f\n"%(cos_bc,cos_ac,cos_ab))

  # default input for the block "system"
  ofile.write("  ecutwfc = 20.0 \n/\n") 

  # block electrons 
  ofile.write("""&electrons
  diagonalization='cg'
  mixing_mode = 'plain'
  mixing_beta = 0.7
  conv_thr =  1.0d-8
/
""")

  # block ions  
  ofile.write("""&ions
  ion_dynamics = 'bfgs'
  ion_positions = 'default'
/
""")

  # block cell
  ofile.write("""&cell

/
""")
  
  # card ATOMIC_SPECIES 
  ofile.write("ATOMIC_SPECIES\n")
  for i in range(ntyp):
    atom= species[i]
    mass = f_Element_Symbol_to_Mass(atom) 
    psp = atom + qe_default_psp.strip() + ".UPF"
    ofile.write("  %-4s %6.2f %s\n"%(atom,mass,psp)) 
 
  # card ATOMIC_POSITIONS

  if latt is None: 
    ofile.write("ATOMIC_POSITIONS angstrom\n")
  else:
    ofile.write("ATOMIC_POSITIONS crystal\n")

  for i in range(nat):
    atom=mol[i][0]
    x = mol[i][1][0]
    y = mol[i][1][1]
    z = mol[i][1][2] 
    ofile.write("  %-4s %12.6f %12.6f %12.6f\n"%(atom,x,y,z))

  # card K_POINTS
  ofile.write("K_POINTS automatic\n")
  if latt is None:
    ofile.write(" 1 1 1 0 0 0 \n")
  else:
    # estimate the number of k-points needed 
    nk=[1,1,1]
    for i in range(3):
      nk[i] = int(kgrid_cutoff/latt[i])+1
    ofile.write("%2d %2d %2d  1  1  1 \n"%(nk[0],nk[1],nk[2]))
  
  ofile.close()  

def f_QE_DefaultRunCommand(stCase,stExec="pw.x"):
    '''
    Create a command to run a QE program with specific input, output and error files.
    :param stCase: the input filename exclude .in
    :param stExec: the program to run, default pw.x
    '''
    return "%s <%s.in >%s.out 2>%s.err" % (stExec,stCase,stCase,stCase)

def f_QE_Run(stInp,stOut,stExec="pw.x",np=1):
    '''
    Create a command to run a QE program with specific input, output and error files.
    :param stCase: the input filename exclude .in
    :param stExec: the program to run, default pw.x
    '''
    qejob= "%s <%s >&%s" % (stExec,stInp,stOut)
    f_env_RunMpirunCommand(qejob, np)

def QESP_GetVersion(object):
    '''
    Get version number of Quantum-Espresso
    Note some input files are different in QE4 and QE5, like dos.x and projwfc.x
    To set version number used in tmckit, set enviroment variable tmckit_qe_version=4 or 5
    By default it is 4
    '''
    nDefault = 4
#Read enviroment variable
    if ( os.environ.has_key("tmckit_qe_version")):
        return int(os.environ["tmckit_qe_version"])
    else:
#If there is none, then try to detect and set it
#Try run pw.x
        err = subprocess.call(["which","pw.x"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if ( err != 0):#No pw.x found, return default 4
            return nDefault
        else:
#Try read from pw.x output
            process = subprocess.Popen(["pw.x"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            time.sleep(2)
            process.terminate()
            time.sleep(1)
            output,outerr = process.communicate()
            err = process.poll()
            if ( err == None):
                raise RuntimeError,"Cannot execute pw.x correctly!"
            else:
                nIndex = output.index("v.")
                result = int(output[nIndex+2])
                print("Detected Quantum-Espresso version : %i" % result)
                os.environ["tmckit_qe_version"] = str(result)
                return result
    pass


class QESP_ld1_input(NamelistIO):
    def __init__(self,stFileName = None):
        super(QESP_ld1_input,self).__init__()
        self.title = "H"
        self.zed = 1.0
        self.rel = 1
        self.config = ""
        self.iswitch=3
        self.dft="LDA"
        self.lloc = 2
        self.pseudotype=3
        self.filepseudo_pw = ''
        self.author = "TMC"
        self.PPCard = []
        if ( stFileName != None):
            self.ReadFromFile(stFileName)
    
    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        list_stLine = fIn.readlines()
        i = 0
        stCurrentPart = ""
        while ( i < len(list_stLine)):
            stLine = list_stLine[i].strip()
            if ( len(stLine) == 0): # skip empty line
                i += 1
                continue
                        
            if( stLine[0] == "/"):#List part
                if ( stCurrentPart == "inputp"): # after inputp part
                    i += 1
                    while ( list_stLine[i].strip() == ""):
                        i += 1
                    stCount = int(list_stLine[i].strip())
                    i += 1
                    for j in range(i,i+stCount):
                        orb = list_stLine[j].split()
                        orb = [orb[0],int(orb[1]),int(orb[2]),float(orb[3]),float(orb[4]),float(orb[5]),float(orb[6])]
                        if ( len(orb) == 8 ): # for rel != 1
                            orb.append(float(orb[7]))
                        else:
                            orb.append(0.0)
                        self.PPCard.append(orb)
                    break
                else:
                    i += 1
                    continue
                
            
            # normal part

                
            if ( stLine[0] == '&'): #mark line, skip
                if ( stLine[0] == '&'):
                    stCurrentPart = stLine[1:].strip()
                    self.dicExistName[stCurrentPart]= []
                i += 1
                continue
            # comma split
            if ( stLine.find(",") != -1): # split ','
                del list_stLine[i]
                for stSplit in stLine.strip().split(","):
                    list_stLine.insert(i,stSplit)
                continue
            stCurrentName,var_tmp = self.__GetNameValue__(stLine)
            if ( not stCurrentName in  self.dicExistName[stCurrentPart]):
                self.dicExistName[stCurrentPart].append(stCurrentName)
            i+=1
        fIn.close()
        
    def WriteToFile(self,stFileName):
        fOut = open(stFileName,'w')
        
        self.__WriteNameValue__(fOut,"input",[])
        if ( "inputp" in self.dicExistName.keys()):
            self.__WriteNameValue__(fOut,"inputp",[])
        if ( "test" in self.dicExistName.keys()):
            self.__WriteNameValue__(fOut,"test",[])
       #self.__WriteNameValue__(fOut,"input",["title","zed","rel","config","iswitch","dft"])
        #self.__WriteNameValue__(fOut,"inputp",["lloc","pseudotype","file_pseudopw","author"])
        #PP Generation Card
        if ( len(self.PPCard) > 0):
            fOut.write(str(len(self.PPCard)))
            if ( self.rel != 1 ):
                fOut.write("%6.2f")
            fOut.write("\n")
            
            for PPOrbital in self.PPCard:
                fOut.write("%s  %d  %d%6.2f%6.2f%6.2f%6.2f" % tuple(PPOrbital[:-1]))
                if ( self.rel != 1 ):
                    fOut.write("%6.2f" % PPOrbital[-1] )
                fOut.write("\n")
    #fOut.write(stLine+"\n")
            
        fOut.close()

class QESP_ld1_conf(object):
    '''
    The single configuration output part of ld1.x
    '''
    def __init__(self):
        self.dEtotAE = 0.0
        self.dEtotPS = 0.0
        self.dEtotAbsoluteDiff = 0.0 # the absolute difference of total energy between AE and PS. Its absolute value is not important; but important when compare with other configuration
        self.listAEOrbital = []
        self.listPPDiff = []
        self.listPPTest = [] #only the test with highest cutoff Bessel functions is record
        self.listConfOrb = [] # configuration occupation

       
    def Read(self,list_stLine,nStart=-1):

        global step,step_prev

        def StepNew(stepNew):
            global step,step_prev
            step_prev = step
            step = stepNew
            #print(step)

        i  = nStart
        
        listC = []
        step_none = 0 #not in usable part
        step_ae = 1
        step_pp = 2
        step_ppb = 3
        step_list = 4
        step_list2 = 5
        step_list_pp = 6
        step_pp_ene = 8
        step_list_ppb = 7
        step = step_none
        step_prev = step_none
    
        dicEmptyLineEnd = {step_list:step_ae,step_list2:step_none,step_list_pp:step_pp,step_list_ppb:step_ppb}# end these part if read empty line; if end, goto its value 

        stPart = ""
        while (i < len(list_stLine)-1):
            i += 1
            stLine = list_stLine[i]
            if ( stLine.strip() == ""):
                if ( step in dicEmptyLineEnd.keys()):
                    StepNew(dicEmptyLineEnd[step])
                else:
                    continue
            if ( step == step_none ):
                if ( "---" in stLine and not "End" in stLine):
                    if ( "All-electron" in stLine):
                        StepNew(step_ae)
                        listC = self.listAEOrbital
                    elif ( "Generating" in stLine):# not used here
                        i += 1
                    elif ( "Testing the pseudopotential" in stLine):
                        StepNew(step_pp)
                        listC = self.listPPDiff
                    elif ( "Bessel" in stLine):
                        StepNew(step_ppb)
                        listC = self.listPPTest
                    else:
                        raise FormatError,"Unknow content %s of ld1.x output" % stLine.strip()
            elif ( step == step_ae ):
                if ( "n l" in stLine):
                    StepNew(step_list)
                    listC = self.listAEOrbital
                elif ( "normal" in stLine):
                    StepNew(step_list2)
                    listC = self.listAEOrbital
                    i += 1

            elif ( step == step_list): #read energy
                listC.append([int(stLine[5]),int(stLine[7]),stLine[13:15],float(stLine[18:22]),float(stLine[56:69])])

            elif ( step == step_list2): #read corresponding r
                stName = stLine[7:9]
                if ( stName != stLine[10:12] ): #different orbital overlap, neglect
                    continue
                listTmp = None
                for Orbital in listC:
                    if ( Orbital[2] == stName ):
                        listTmp = Orbital
                        break
                if ( listTmp == None):
                    raise ValueError,"%s does not exist in orbital energy list"
                listTmp += [float(stLine[33:41]),float(stLine[49:59]),float(stLine[70:78])]

            elif ( step == step_pp):
                if ( "n l" in stLine):
                    StepNew(step_list_pp)
                    listC = self.listPPDiff
                elif ( "Etot =" in stLine):
                    self.dEtotAE = float(stLine[11:27])
                elif ( "Etotps =" in stLine):
                    self.dEtotPS = float(stLine[13:27])
                    self.dEtotAbsoluteDiff = self.dEtotPS-self.dEtotAE
                elif ( "End" in stLine):
                    StepNew(step_none)

            elif ( step == step_list_pp):
                listC.append([int(stLine[5]),int(stLine[7]),stLine[13:15],float(stLine[20:24]),float(stLine[29:42]),float(stLine[42:57]),float(stLine[58:72])])


            elif ( step == step_ppb ):
                if ( "End" in stLine):
                    break #end of read
                if ( "Cutoff" in stLine):
                    dCut = float(stLine[18:25])
                    self.listPPTest.append([dCut])
                    listC = self.listPPTest[-1]
                    i += 1
                    StepNew(step_list_ppb)

            elif ( step == step_list_ppb): 
                listC.append([int(stLine[9]),float(stLine[19:28]),float(stLine[32:41]),float(stLine[45:54])])

        #create occupation description
        self.listConfOrb = []
        for Orbital in self.listPPDiff:
            self.listConfOrb.append([Orbital[2].lower(),Orbital[3]])
        return i

class QESP_ld1_output(object):
    '''
    The output of ld1.x 
    Note: the test part may contain a lot of seperate configuration test result
    '''
    def __init__(self,stFileName=None):
        self.listConf = []
        self.bOK = False 
        self.stErrorPos = ""
        self.stError = ""
        if ( stFileName != None):
            self.ReadFromFile(stFileName)

    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        list_stLine = fIn.readlines()
        fIn.close()

        #Detect error
        self.bOK = True
        for i in range(len(list_stLine)-1,0,-1):
            if ( "%%%" in list_stLine[i]):
                self.bOK = False
                self.stError =  list_stLine[i-1].strip()
                self.stErrorPos = list_stLine[i-2].strip()
                return
            elif ( "warning:" in list_stLine[i]):
                self.bOK = False
                self.stError =  list_stLine[i].strip()
                self.stErrorPos = list_stLine[i-1].strip() 
            elif ( "Errors" in list_stLine[i]):
                self.bOK = False
                self.stError = list_stLine[i].strip()
                self.stErrorPos = list_stLine[i-1].strip() 

        i = -1
        stPart = ""
        step_none = 0
        step_e = 1
        step_r = 2
        step_list = 3
        step_list2 = 4
        step = step_none
        
        self.listConf = []

        while ( i < len(list_stLine)-1):
            i += 1 
            stLine = list_stLine[i]
            if ( "---" in stLine ):
                #print("Start Conf Part: %d" % i)
                i -= 1
                Conf = QESP_ld1_conf()
                i = Conf.Read(list_stLine,i)
                self.listConf.append(Conf)
    


class QESP_projwfc_input(NamelistIO):
    '''
    The input file of projwfc.x
    '''
    def __init__(self,stFileName=None):
        '''
        :param stFileName:: the file to read when initialize.
        '''
        super(QESP_projwfc_input,self).__init__()

#Detect QE version
        nVersion = QESP_GetVersion()
        stTemplate ="template-qe/projwfc.in" 
        if ( nVersion == 4):
            stTemplate = "template-qe/projwfc-v4.in"
        
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),stTemplate))

        if ( stFileName != None):
            self.ReadFromFile(stFileName)
    
    #def WriteToFile(self,stFileName):
    #    self.__DetectFileExist__(stFileName)
    #    fOut = open(stFileName,'w')
    #    self.__WriteNameValue__(fOut,"inputpp",[])

class QESP_dos_input(NamelistIO):
    '''
    The input file of dos.x
    '''
    def __init__(self,stFileName=None):
        '''
        :param stFileName:: the file to read when initialize.
        '''
        super(QESP_dos_input,self).__init__()

        nVersion = QESP_GetVersion()
        stTemplate ="template-qe/dos.in" 
        if ( nVersion == 4):
            stTemplate = "template-qe/dos-v4.in"
        
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),stTemplate))
        if ( stFileName != None):
            self.ReadFromFile(stFileName)

class QESP_bands_input(NamelistIO):
    '''
    The input file of bands.x
    '''
    def __init__(self,stFileName=None):
        '''
        :param stFileName:: the file to read when initialize.
        '''
        super(QESP_dos_input,self).__init__()

        self.ReadFromFile(os.path.join(os.path.dirname(__file__),"template-qe/bands.in"))
        if ( stFileName != None):
            self.ReadFromFile(stFileName)


class QESP_PW_input(NamelistIO):
    def __init__(self,stFileName=""):
        '''
        Create a Quantum-Espresso Input file object
        :param stFileName:: The file to read when initialize. If not set, use default value.
        '''
        super(QESP_PW_input,self).__init__()

        self.listAtomSpecies = []
        self.stAtomPositionMode = "alat"
        self.listInnerCoord = [] #[7], elements are AtomName,x,y,z,xfix,yfix,zfix
        self.listAtomForce = []
        self.listKPt = []
        self.listKPtGrid = [0,0,0,0,0,0]
        self.stKPtMode = "automatic"
        self.aCell = Lattice() #IMPORTANT NOTE: the coordinate conversion is based on Lattice(), __GetCell__() must be used before any attempt to convert coordinate format
        self.cell_parameters = None #! 3 cell vectors, only used when ibrav=0, this is always stored as alat
        self.unit_cell_parameters = None #! Unit of cell vectors
                
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),"template-qe/pw.in"))
        #clear default read
        
        if ( stFileName != ""):
            self.dicExistName = {}
            self.ReadFromFile(stFileName)
            
    def __refresh_celldm__(self):
        '''
        clean celldm that not used ( 1 is always used )
        '''
#Don't use celldm with ibrav==0
        if (self.ibrav == 0):
            return
        list_used_celldm = [\
                            [],\
                            [],\
                            [],\
                            [],\
                            [3],\
                            [4],\
                            [3],\
                            [3],\
                            [2,3],\
                            [2,3],\
                            [2,3],\
                            [2,3],\
                            [2,3,4],\
                            [2,3,4],\
                            [2,3,4,5,6]\
                            ]
        
        celldm_new = [-1 for x in range(0,100)]
        celldm_new[1] = self.celldm[1]
        for nIndex in list_used_celldm[self.ibrav]:
            celldm_new[nIndex] = self.celldm[nIndex]
        self.celldm = celldm_new

    def __write_cards__(self,fOut,stLattMode,stCellParaMode):
        '''
        Write cards into the file
        :param fOut: the opened file
        :param stLattMode: specify how to store atom information ( default: the form stored )
        :param stCellParaMode: specify how to store CELL_PARAMETERS (default: the form stored)
        '''
        #CELL_PARAMETERS only used when ibrav=0
        if (self.ibrav == 0):
            cell2 = self.cell_parameters
            if (stCellParaMode != self.unit_cell_parameters):
#Convert unit
                if (self.unit_cell_parameters == "alat"):
                    cell2 = f_Matrix_Op_Scalar(cell2,"*",self.celldm[1])
                elif (self.unit_cell_parameters == "angstrom"):
                    cell2 = f_Matrix_Op_Scalar(cell2,"*",Ang2Bohr)

                if (stCellParaMode == "alat"):
                    cell2 = f_Matrix_Op_Scalar(cell2,"/",self.celldm[1])
                elif (stCellParaMode == "angstrom"):
                    cell2 = f_Matrix_Op_Scalar(cell2,"/",Ang2Bohr)

            fOut.write("CELL_PARAMETERS %s\n" % stCellParaMode)
            for vec in cell2:
                fOut.write("%16.8f %16.8f %16.8f\n" % tuple(vec))
            fOut.write("\n")

        fOut.write("ATOMIC_SPECIES\n")
        for aAtom in self.listAtomSpecies:
            fOut.write(aAtom.ToString()+"\n")
        fOut.write("\n")
        
        #atom position part
        if ( stLattMode == None ):
            fOut.write("ATOMIC_POSITIONS %s\n" % self.stAtomPositionMode)
            for aCoord in self.listInnerCoord:
                fOut.write("%4s %16.8f %16.8f %16.8f %d %d %d\n" % tuple(aCoord) )
        else:
            fOut.write("ATOMIC_POSITIONS %s\n" % stLattMode)
            #Convert atom format
            if ( stLattMode == "crystal"):
                listCoord = [aAtom+[1,1,1] for aAtom in self.aCell.listAtom]
            else:
                listCoord = self.aCell.listAtomCartesian
                if ( stLattMode == "alat"):
                    listCoord = [ [aAtom[0]]+[x/self.celldm[1] for x in aAtom[1:4]] + [1,1,1] for aAtom in listCoord]
                elif ( stLattMode == "angstrom"):
                    listCoord = [ [aAtom[0]]+[x*Bohr2Ang for x in aAtom[1:4]] + [1,1,1] for aAtom in listCoord]
            for i in range(0,len(listCoord)):
                listCoord[i][4:7] = self.listInnerCoord[i][4:7]
            for aCoord in listCoord :
                fOut.write("%4s %16.8f %16.8f %16.8f %d %d %d\n" % tuple(aCoord) )
        fOut.write("\n")
        
        fOut.write("K_POINTS %s\n" % self.stKPtMode)
        if ( self.stKPtMode == "automatic" ):
            fOut.write("%d %d %d %d %d %d\n" % tuple(self.listKPtGrid))
        elif ( self.stKPtMode == "gamma"):
            pass
        else:
            fOut.write("%d\n" % (len(self.listKPt)))
            for aKPt in self.listKPt:
                fOut.write("%11.8f %11.8f %11.8f %4.2f\n" % tuple(aKPt))
        
        fOut.write("\n")
    
    def WriteToFile(self,stFileName,stLattMode=None,stCellParaMode=None):
        '''
        Write information into pw.x input file
        :param stFileName: The file to write into ( automatically overwrite)
        :param stLattMode: specify how to store atom information ( default: the form stored )
        :param stCellParaMode: specify how to store CELL_PARAMETERS (default: the form stored)
        '''
        self.GetCell() # Refresh cell for later use
        
        self.__DetectFileExist__(stFileName)
        
        #clean celldm
        self.__refresh_celldm__()     
        
        fOut = open(stFileName,'w')
        #self.__WriteNameValue__(fOut,"control",["calculation","restart_mode","tstress","tprnfor","outdir","prefix","pseudo_dir"])
        self.__WriteNameValue__(fOut,"control",["calculation","restart_mode","tstress","tprnfor","outdir","prefix","wf_collect","verbosity"])
        #self.__WriteNameValue__(fOut,"system", ["nat","ntyp","ibrav","celldm","ecutwfc","input_dft","smearing","degauss","occupations"])
        self.__WriteNameValue__(fOut,"system", ["nat","ntyp","ibrav","celldm","ecutwfc"])
        self.__WriteNameValue__(fOut, "electrons",["diagonalization","mixing_mode","mixing_beta","conv_thr"])
        self.__WriteNameValue__(fOut,"ions",["ion_dynamics","ion_positions"])
        
        self.__WriteNameValue__(fOut,"cell",[])
        #fOut.write("&cell\n\n\n/\n\n")

       #Write cards 
        self.__write_cards__(fOut,stLattMode,stCellParaMode)

#       print("Struct write into %s" % stFileName)
        fOut.close()
        
    @deprecated    
    def __GetQESPCartesianLattice__(self,ibrav,celldm):
        listLength = [celldm[0],celldm[0],celldm[0]]
        listAngle = [math.pi/2,math.pi/2,math.pi/2]
        if ( ibrav == 0):
            return None
        elif ( ibrav <= 3): # default : cubic, no change
            pass
        elif ( ibrav == 4 or ibrav == 6 or ibrav == 7):
            listLength[2] *= celldm[2]
            listAngle[0] = [math.pi*2/3]
        elif ( ibrav == 5):
            listAngle = [math.acos(celldm[3]) for x in range(0,3)]
        elif ( ibrav <= 11): # ortho
            listLength[1] *= celldm[1]
            listLength[2] *= celldm[2]
        elif ( ibrav <= 13): #monoclinic
            listLength[1] *= celldm[1]
            listLength[2] *= celldm[2]
            listAngle[0] = math.acos(celldm[3])
        elif ( ibrav == 14): # triclinic
            listLength[1] *= celldm[1]
            listLength[2] *= celldm[2]
            listAngle = [math.acos(x) for x in celldm[3:6]]                       
        
        
    def GetCell(self):
        '''
        Get cell information from input file 
        '''
        if ( not hasattr(self,"ibrav")):
            return None;
        
        self.aCell = Lattice()
        a = self.celldm[1]
        if (self.ibrav != 0): #Read celldm
            #Align celldm as a,b,c,alpha,beta,gamma
            if ( self.ibrav ==12 or self.ibrav == 13):
                self.celldm[6] = self.celldm[4]
            self.aCell.ReadFromABC(self.ibrav, [a,self.celldm[2]*a,self.celldm[3]*a], [math.acos(x) if x <> -1 else 0 for x in self.celldm[4:7]],unit_angle = "rad")
        #print(self.listInnerCoord)
        else: #Read vectors
            cell2 = self.cell_parameters
            if (self.unit_cell_parameters == "alat"):
                cell2 = f_Matrix_Op_Scalar(cell2,"*",self.celldm[1])
            elif (self.unit_cell_parameters == "angstrom"):
                cell2 = f_Matrix_Op_Scalar(cell2,"*",Ang2Bohr)

            self.aCell.ReadFromRawCellVector(cell2)

        self.aCell.AddAtom([Atom[0:4] for Atom in self.listInnerCoord],self.stAtomPositionMode,a=a)
        return self.aCell
    
    def ReadFromCell(self,Latt):
        '''
        Read Cell Information into a pw.x input file
        Ref : doc/INPUT_PW.html
        '''
        self.aCell = copy.deepcopy(Latt)
        
        self.celldm = [-1 for x in range(0,100)]
        self.celldm[1] = Latt.fLatticeLength[0]
        if (Latt.IsRaw):
            self.ibrav = 0
            self.cell_parameters = copy.deepcopy(Latt.PrimitiveCellVector)
            self.unit_cell_parameters = "bohr"
        else:
            self.ibrav = Latt.nLatticeType
            arLength = [x  for x in Latt.fLatticeLength]
            arCosAngle = [math.cos(x /180.0 * math.pi) for x in Latt.fLatticeAngle ]
            
            if ( self.ibrav == 4 or self.ibrav == 6 or self.ibrav == 7): 
                self.celldm[3] = arLength[2]/self.celldm[1]
            elif ( self.ibrav == 5):
                self.celldm[4] = arCosAngle[0]
            elif ( self.ibrav >= 8):
                self.celldm[2] = arLength[1]/self.celldm[1]
                self.celldm[3] = arLength[2]/self.celldm[1]
                if ( self.ibrav >= 12):
                    self.celldm[4] = arCosAngle[2]
                    if ( self.ibrav == 14):
                        self.celldm[4] = arCosAngle[0]
                        self.celldm[5] = arCosAngle[1]
                        self.celldm[6] = arCosAngle[2]
        #Read Atom ( stored as alat form )
        self.nat = len(Latt.listAtom)
        listAtomSpecies = []
        self.listInnerCoord = []
        
        self.stAtomPositionMode = "alat"
        for aAtom in Latt.listAtomCartesian:
            self.listInnerCoord.append([aAtom[0]] + [x/self.celldm[1] for x in aAtom[1:4]]+ [1,1,1])
            #Count atom species
            bExist = False
            for aKind in listAtomSpecies:
                if ( aKind.stName == aAtom[0]):
                    bExist = True
                    break
            if ( not bExist ):
            #Check if already has this kind of atom before 
                aAtomOrg = None
                for aKind in self.listAtomSpecies:
                    if ( aKind.stName == aAtom[0]):
                        aAtomOrg = copy.deepcopy(aKind)
                        break
                if ( aAtomOrg != None):
                    listAtomSpecies.append(aAtomOrg)
                else:
                    aAtomNew = QESP_Atom()
                    aAtomNew.stName = aAtom[0]
                    aAtomNew.stPP =  aAtom[0]+".pbe-rrkjus.UPF"
                    listAtomSpecies.append(aAtomNew)
        self.listAtomSpecies = listAtomSpecies
        self.ntyp = len(listAtomSpecies)
        
    def ReadFromFile(self,stFileName):
        '''
        Read in content from input file
        Record all name readed in list to write later.
        '''
        fIn = open(stFileName,'r')
        list_stLine = fIn.readlines()
        self.list_stLine = list_stLine
        i = 0
        stCurrentPart = None
        while ( i<len(list_stLine) ):
            stLine = list_stLine[i].strip()
            #print(i,stLine)      
            #end part

            #Read cards
            if( stLine.strip() == "ATOMIC_SPECIES"):
                #Read atom species
                self.listAtomSpecies = []
                i += 1
                for j in xrange(i,i+self.ntyp):
                    self.listAtomSpecies.append(QESP_Atom(list_stLine[j]))

                i += self.ntyp
                stLine = list_stLine[i]
                continue

            elif (stLine.find("ATOMIC_POSITIONS") != -1):
                #Read atom coordinates
                arLine = stLine.strip().split()
                self.stAtomPositionMode = "alat" if len(arLine)==1 else arLine[1].replace('{','').replace('}','')
                self.listInnerCoord = []
                i += 1
                for j in xrange(i,i+self.nat):
                    stLine = list_stLine[j]
                    arLine = stLine.strip().split()
                    if ( len(arLine) == 4): #no fixed position information
                        arLine = arLine + ['1','1','1']
                    listCoord = [arLine[0]] + [float(x) for x in arLine[1:4]] + [int(x) for x in arLine[4:]]
                    self.listInnerCoord.append(listCoord)

                i += self.nat
                stLine = list_stLine[i]

                #print(i,stLine)                
                continue
            elif (stLine.find("K_POINTS") != -1):
                #K point
                arLine = stLine.strip().lower().split()
                self.listKPt = []
                self.listKPtGrid = [0,0,0,0,0,0]
                self.stKPtMode = "tpiba" if len(arLine)==1 else arLine[1].replace('{','').replace('}','')
                i+=1
                if (  self.stKPtMode == "tpiba" or  self.stKPtMode == "crystal" or  self.stKPtMode == "tpiba_b" or  self.stKPtMode == "crystal_b"):
                    nCount  = int(list_stLine[i].strip())
                    i += 1
                    for j in range(i,i+nCount):
                        self.listKPt.append([float(x) for x in list_stLine[j].strip().split()])
                elif( arLine[1] == "automatic"):
                    self.listKPtGrid = [int(x) for x in list_stLine[i].strip().split()]
                elif ( arLine[1] == "gamma"):
                    pass
                continue
            elif (stLine.find("CELL_PARAMETERS") != -1):
                i += 1
                if (self.ibrav != 0): #Skip
                    i += 3
                    continue

                ar = stLine.strip().lower().split()
                if (len(ar) == 1):
                    self.unit_cell_parameters = "bohr"
                else:
                    self.unit_cell_parameters = ar[-1]
                self.cell_parameters = []
                for j in xrange(3):
                    self.cell_parameters.append([float(x) for x in list_stLine[i+j].strip().split()])
                
                i += 3
                continue
            
            # normal part
            i = self.__ProcessLine__(i,stLine)

        fIn.close()
        self.list_stLine = []
        #print(self.dicExistName)
        
class QESP_PW_output:
    '''
    Class for reading output of pw.x
    Please verify bOK, bFinalReached and bNotConverge to see whether the result is converged
    Possible information of end with error results are also read
    '''
    def __init__(self,stFileName = None,bNoSym = False):
        self.stMode = "" # the calculation mode, currently useless
        self.nProcess = 1; #the Process Count in the calculation
        #self.nWALLTime = 0; # the real time used ( in second )         
        #self.nCPUTime = 0; # the CPU time used  ( in second )
        self.stWALLTime = "" # the real time used
        self.stCPUTime = "" # the CPU time used
        self.arEnergy= [] # all the energy in the relax process; final one is final result
        self.arLattice = [] # all the lattice in the relax process; note final one is final result
        self.arBand = [] # all band information in the relax process ; final one is final result
        self.arBandKPt = [] # all k-points of band information in the relax process ; final one is final result. Unit is cartesian  
        self.arForce = [] # all the force in the relax process
        self.arkVector = [] # the first reciporcal vectors
        self.bEnd = False # if this calc ends normally
        self.bFinalReached = False # if this calc reach end structure in optimization
        self.bNotConverge = False # if scf calculation does not converge, but ends normally
        self.bOK = False
        self.nElectron = -1 # number of electrons
        self.nKSStates = -1 #number of states ( bands)
        self.lattInit = None # the initial lattice ( depreceted)
        self.bNoSym = bNoSym # Symmetry is not used in this calculation ( Triclinc cell will be used in all vc-relax results)
        self.fermi = None # Fermi Energy
        self.vbm = None # VBM Energy
        self.b_metal = None # Whether it is a metal
        self.n_spin = 1 # Number of spin
        if ( stFileName != None):
            self.ReadFromFile(stFileName)

    def FirstEnergy(self):
        '''
        Energy of first structure in relaxation
        '''
        if ( len(self.arEnergy) > 0):
            return self.arEnergy[0]
        else:
            return None

    @property
    def FirstLattice(self):
        '''
        The first structure (or the only one if not relaxed)
        '''
        return self.arLattice[0]
    
    def FinalEnergy(self):
        '''
        Energy of last structure in relaxation
        '''
        if ( self.bOK):
            return self.arEnergy[-1]
        else:
            return None
    
    def FinalLattice(self):
        '''
        Last Structure in relaxation
        '''
        if ( self.bFinalReached):
            return self.arLattice[-1]
        else:
            return None
    
    def __TimeToSecond__(self,stTime):
        '''
        Convert a string like in pwscf output like xd ah bm cs to second
        '''
        dicTime = {"d":3600*24,'h':3600,'m':60,'s':1}
        re_time = re.compile("(.*)(\S)(.*)(\S)")
        aMatch = re_time.match(stTime)
        dTime = float(aMatch.group(1))*dicTime(aMatch.group(2))+float(aMatch.group(3))*dicTime(aMatch.group(4))
        return dTime    
    
    def ReadFromFile(self,stFileName):
        #step mark
        step_process = 0 # read process
        step_init_latt = 1 # read initial lattice
        step_scf = 2 # read scf energy
        step_latt = 3 # read new lattice
        
        err = 0
        
        fIn = open(stFileName,"r")
        list_stLine = fIn.readlines()
        
        #Detect whether it ended normally
        self.bEnd = False
        
        if ( len(list_stLine) < 50 ): # treat output < 50 line as not converg
            return
        
#Read QE version
        version_qe = None
        for line in list_stLine[:100]:
            if ("Program" in line):
                version_qe = int(line.split()[2][2:][0])
                
        for i in range ( len(list_stLine)-1, len(list_stLine)-20,-1):
            stLine = list_stLine[i]
            if ( stLine.find("JOB DONE") != -1 ): # find end tag
                self.bEnd = True
            if ( self.bEnd): # find time
                if ( stLine.find("PWSCF") != -1 ):
                    re_time  = re.compile(".*PWSCF.*:(.*)CPU(.*)WALL.*")
                    aMatch = re_time.match(stLine)
                    #self.nCPUTime = self.__TimeToSecond__(aMatch.group(2))
                    #self.nWALLTime = self.__TimeToSecond__(aMatch.group(3))
                    self.stCPUTime = aMatch.group(1).strip()
                    self.stWALLTime = aMatch.group(2).strip()
                    break
        
        if ( not self.bEnd):
            print (stFileName + " does not terminate normally, please check your result.")
            
        #Detect whether final relax structure is obtained
        self.bFinalReached = False
        for i in range(len(list_stLine)-1,0,-1):
            stLine = list_stLine[i]
            if ( stLine == "Begin final coordinates\n"):
                self.bFinalReached = True
                break
        if ( not self.bEnd):
            if ( self.bFinalReached ):
                print("But final structure has been found. Reading Continue")
            else:
                print("Reading stop.")
                err =  1

        if ( self.bEnd or self.bFinalReached): # detect whether SCF converge
            self.bNotConverge = False
            for i in range ( len(list_stLine)-100, len(list_stLine)-1):
                stLine = list_stLine[i]
                if ( stLine.find("convergence NOT") != -1):
                    print("SCF does not converge. Please check you result.")
                    self.bNotConverge = True
                    break
            
        #Judge whether calculation useful
        if ( self.bEnd and not self.bNotConverge):
            self.bOK = True
        else:
            self.bOK = False
        
        #Try read fermi energy or VBM energy
        stFermi = commands.getoutput("grep -a Fermi %s" % stFileName)
        if ( len(stFermi.strip()) != 0):
#Fermi energy
            if ("spin" in stFermi):
                self.fermi = float(stFermi[42:49]) #Only spin up one is read
                self.b_metal = None #We don't know if it is metal or not
            else:
                self.fermi = float(stFermi[26:36])
                self.b_metal = True
        else:
#VBM and CBM
            stFermi = commands.getoutput("grep -a \"highest occupied\" %s" % stFileName)
            if ( len(stFermi.strip()) != 0):
                l1 = stFermi.split()
#Find the number in the line because QE4 and QE5 use different format so no fixed position
                for st2 in l1:
                    if (st2[-1].isdigit()):
                        self.vbm = float(st2)
                        break
                self.b_metal = False

        
        nStep = step_process
        i = -1
        dA = 0 # A axis length
        nAtom = 0 # Atom Count
        while (i < len(list_stLine)-1 ):
            i += 1
            stLine = list_stLine[i]
            if ( nStep == step_scf):# Find SCF end
                if ( "End of" in stLine and "calculation" in stLine):
#Read band structure
                    if (self.n_spin == 1):
                        i += 2
                        if ("SPIN" in list_stLine[i]):#Correct to from 1 to 2
                            self.n_spin = 2
                            i += 3
                    else: #Skip  --- part
                        i += 5
                    stLine = list_stLine[i]
                    list_kp = []
                    list_band = []
                    list_band1 = []
                    while ( "band" in stLine ):
#Note qe's tpiba are always based on initial cell 
#Convert it to cartesian
                        list_kp.append([float(stLine[13:20])*2*pi/dA,float(stLine[20:27])*2*pi/dA,float(stLine[27:34])*2*pi/dA])
                        i += 2
                        stLine = list_stLine[i]
                        list_band.append([])
                        while ( len(stLine.strip()) != 0):
                            list_band[-1] += [float(x) for x in stLine.split()]
                            i += 1
                            stLine = list_stLine[i]
                        i += 1
                        stLine = list_stLine[i]
#Maybe occupation number is printed, if so, then skip this part
                        if ( "occupation" in stLine):
                            while( len(stLine.strip()) != 0):
                                i += 1
                                stLine = list_stLine[i]
#Read to next line
                            i += 1
                            stLine = list_stLine[i]
#Spin down
                        if (self.n_spin == 2 and "SPIN DOWN" in stLine):
                            list_band1 = list_band
                            list_band = []
                            list_kp = [] #We directly clean this as k-points are the same
                            i += 3
                            stLine = list_stLine[i]


                    if (self.n_spin == 2):
                        self.arBand.append([list_band1, list_band])
                    else:
                        self.arBand.append(list_band)
#                   print("Get Kp", len(list_kp))
                    self.arBandKPt.append(list_kp)

#Read total energy information
                if ( stLine[0] == "!"):  #End 
                    self.arEnergy.append(float(stLine[32:50]))
                    #find "convergence"
                    while ( not "convergence" in stLine):
                        i += 1
                        stLine = list_stLine[i]
                    #read force 
                    if ( "Force" in list_stLine[i+2]):
                        i += 2
                        listForce = []
                        while ( not "atom" in stLine):
                            i += 1
                            stLine = list_stLine[i]
                        for j in range(0,nAtom):
                            #print(stLine)
                            listForce.append([float(x) for x in stLine.split()[-3:]])
                            i += 1
                            stLine = list_stLine[i]
                        #dispersion force ( +D correction)
                        if ( "Dispersion" in list_stLine[i+5]):
                            i += 5
                            for j in range(0,nAtom):
                                i += 1
                                stLine = list_stLine[i]                                
                                arDisForce = [float(x) for x in stLine.split()[-3:]] 
                                #listForce[j][0] += arDisForce[0]
                                #listForce[j][1] += arDisForce[1]
                                #listForce[j][2] += arDisForce[2]
                        self.arForce.append(listForce)
                    nStep = step_latt
                    continue
            elif ( nStep == step_latt):
                if ( stLine.find("Final energy") != -1 ): # Final reached
                    self.arEnergy.append(float(stLine[22:41]))
                    continue            
                if ( stLine[0:4] == "CELL" or stLine[0:4] == "ATOM" ): # Start Cell parameter or Atom parameter ( for vc-relax and relax, respectively )
                    Latt =  copy.deepcopy(self.lattInit)
                    if ( stLine[0:4] == "CELL"): #always alat 
                        i += 1
                        arCell = list_stLine[i:i+3]
                        mCellVector = [ [float(x)*dA for x in stLine.strip().split()] for stLine in arCell]
                        if ( self.bNoSym ):
                            Latt.nLatticeType = 14
                            Latt.RawCellVector = mCellVector
                        Latt.ReadFromPrimitiveCellVector(mCellVector)
                        dANew = Latt.fLatticeLength[0]
                        i += 4
                        stLine = list_stLine[i]
                    else:
                        dANew = dA # relax doesn't change A axis
                    stAtomMode = stLine[18:22]
                    if (stAtomMode == "angs"):
                        stAtomMode = "ang"
                    stLatt = "raw" if self.bNoSym else "prim"
                    for j in range(0,nAtom):
                        i += 1
                        stLine = list_stLine[i]
                        arLine = stLine.strip().split()
                        if ( stAtomMode == "alat"):
                            arLine = [arLine[0]] + [float(x)*dA for x in arLine[1:4]]
                            stAtomModeNew = "bohr"
                        else:
                            arLine = [arLine[0]] + [float(x) for x in arLine[1:4]]
                            stAtomModeNew = stAtomMode
                        Latt.AddAtom([arLine], stAtomModeNew,stLatt)
                        
                    Latt.CorrectError()
                    self.arLattice.append(Latt)
                    nStep = step_scf

                continue
            
            elif ( nStep == step_process):
                if ( stLine.find("processors") != -1):
                    self.nProcess = int(stLine[40:45]) # Read process
                    nStep = step_init_latt
                continue
            elif (nStep == step_init_latt):
                if ( stLine.find("bravais-lattice") != -1):
                    nBravisType = int(stLine[34:45])
                    
                    i += 3 # Read atom/cell 
                    stLine = list_stLine[i]
                    nAtom = int(stLine[34:45])                    
                    
                    i += 2 # Read number of electrons
                    stLine = list_stLine[i]
                    self.nElectron = int(float(stLine[34:45]))
                    if ("up" in stLine):
                        self.n_spin = 2
                    
                    i += 1 # Read states
                    stLine = list_stLine[i]
                    self.nKSStates = int(float(stLine[34:45]))                    
                    
                    celldm = [0,0,0,0,0,0,0] 
                    while ( True):
                        i += 1
                        stLine = list_stLine[i]
                        if ( stLine.find("celldm") != -1):
                            celldm[1] = float(stLine[15:26])
                            celldm[2] = float(stLine[38:49])
                            celldm[3] = float(stLine[61:72])
                            i += 1
                            stLine = list_stLine[i]
                            celldm[4] = math.acos(float(stLine[15:26]))
                            celldm[5] = math.acos(float(stLine[38:49]))
                            celldm[6] = math.acos(float(stLine[61:72]))
                            #pw.x put angle ab at celldm(4) if it is the only angle
#Swap it!
                            if (nBravisType == 12 or nBravisType == 13):
                                celldm[6] = celldm[4]
                                celldm[4] = 0.0
                            dA = celldm[1] # A axis  
                            self.lattInit = Lattice()
                            self.lattInit.ReadFromABC(nBravisType, [celldm[1]] + [dA * x for x in celldm[2:4]], celldm[4:7],unit_angle="rad")
                            break
                        else:
                            continue
                    while ( stLine.find("reciprocal") == -1):
                        i += 1
                        stLine = list_stLine[i]
                    for j in range(1,4):
                        self.arkVector.append([float(x) for x in list_stLine[i+j][24:54].split()])
                    i += 4
                    #read initial atom position
                    while ( not "site" in stLine ):
                        i += 1
                        stLine = list_stLine[i]
                    Latt =  copy.deepcopy(self.lattInit)
                    i += 1
                    for j in range(0,nAtom):
                        stLine = list_stLine[i+j]
                        #print([stLine[17:23].strip(),float(stLine[38:49])*dA,float(stLine[50:61])*dA,float(stLine[62:73])*dA])
                        Latt.AddAtom([[stLine[17:23].strip(),float(stLine[38:49])*dA,float(stLine[50:61])*dA,float(stLine[62:73])*dA]])
                    self.arLattice.append(Latt)                         
                    nStep = step_scf
                continue    
        fIn.close()
        return err            
            
    def ShowSummary(self):
        if ( self.bEnd ):
            print("Result ends normally.")
        elif ( self.bFinalReached):
            print("Result ends abnormally but final structure obtained.")
        else:
            print("Result ends abnormally.")
            return
        
        print("Final SCF energy: %f Ry" % self.arEnergy[-1])
        if ( len(self.arLattice) != 0):
            print("Final Structure:")
            self.arLattice[-1].ShowSummary()
            
class QESP_Atom(object):
    '''
    Represent a atom species definition in QESP
    '''
    def __init__(self,stLine=""):
        self.stName = "none"
        self.dWeight = 1
        self.stPP = "none"
        if ( stLine != ""):
            self.Read(stLine)
    
    def Read(self,stLine):
        arLine = stLine.strip().split()
        self.stName = arLine[0]
        self.dWeight = float(arLine[1])
        self.stPP = arLine[2]
    
    def ToString(self):
        return ("%s %f %s" % (self.stName,self.dWeight,self.stPP))

class QESP_ph_input(NamelistIO):
    '''
    Quantum-Espresso ph.x input file
    '''
    def __init__(self,stFileName=""):
        '''
        Create a Quantum-Espresso ph.x input file object
        :param stFileName:: The file to read when initialize. If not set, use default value.
        '''
        self.stTitle = None
        self.listKPtGrid = []

        self.dicExistName = {}
        self.listExistName = []
               
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),"template-qe/ph.in"))
        if ( stFileName != ""):
            self.dicExistName = {}
            self.listExistName = []
            self.ReadFromFile(stFileName)
 
    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        self.list_stLine = fIn.readlines()
        self.stTitle = self.list_stLine[0].strip()
        i = 1
        self.stCurrentPart = ""
#read inputph
        while ( i < len(self.list_stLine)):
            stLine = self.list_stLine[i].strip()
            i = self.__ProcessLine__(i, stLine)
            if ( self.stCurrentPart == "inputph" and  stLine == "/"):
                break
        while ( self.list_stLine[i].strip() == ""):
            i += 1
        self.listKPtGrid = [int(x) for x in self.list_stLine[i].strip().split()]
        #clean
        self.list_stLine = []

    def WriteToFile(self,stFileName):
        self.__DetectFileExist__(stFileName)

        fOut = open(stFileName,'w')
        fOut.write(self.stTitle+"\n")
        for stPart in self.listExistName:
            self.__WriteNameValue__(fOut,stPart,[])
        fOut.write(" ".join([str(x) for x in self.listKPtGrid]))
        fOut.write("\n")

class QESP_matdyn_input(NamelistIO):
    '''
    matdyn.x input file
    '''
    def __init__(self,stFileName=""):
        '''
        Create a Quantum-Espresso ph.x input file object
        :param stFileName:: The file to read when initialize. If not set, use default value.
        '''
        self.fildyn = None
        self.listKPt = []

        self.dicExistName = {}
        self.listExistName = []
               
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),"template-qe/matdyn.in"))
       
        if ( stFileName != ""):
            self.dicExistName = {}
            self.listExistName = []
            self.ReadFromFile(stFileName)

    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        self.list_stLine = fIn.readlines()
        i = 0
        self.stCurrentPart = ""
        while ( i < len(self.list_stLine)):
            stLine = self.list_stLine[i].strip()
            i = self.__ProcessLine__(i, stLine)
            if ( self.stCurrentPart == "input" and  stLine == "/"):
                break
        while ( self.list_stLine[i].strip() == ""):
            i += 1
        nCount = int(self.list_stLine[i].strip())
        for j in range(i+1,nCount):
            self.listKPt.append( [ float(x) for x in self.list_stLine[i]])
        #clean
        self.list_stLine = []

    def WriteToFile(self,stFileName):
        fOut = open(stFileName,'w')
        self.__DetectFileExist__(stFileName)

        fOut = open(stFileName,'w')
        for stPart in self.listExistName:
            self.__WriteNameValue__(fOut,stPart,[])
        fOut.write(str(len(self.listKPt))+"\n")
        for KPt in self.listKPt:
            fOut.write(" ".join([str(x) for x in KPt]) + "\n")
        fOut.write("\n")


def f_QE_LoadOptimizedCell(stInput,stOutput,stNew,bNoSym=False):
    '''
    After a pw.x vc-relax calculation, load a optimized cell from output file, 
    and create a new input file from old input file parameters and optimized cell
    Also set new input file calculation = "scf"
    '''
    aIn = QESP_PW_input(stInput)
    aOut = QESP_PW_output(stOutput,bNoSym)
    if ( aIn.calculation != "vc-relax" and aIn.calculation != "relax"):
        print("Input file is not a relax/vc-relax calculation, no new input will be generated.")
        return 1
    
    if ( not aOut.bFinalReached):
        print("The calculation does not get an optimized cell, no new input will be generated.")
        return 2
    
    aIn.ReadFromCell(aOut.FinalLattice())
    aIn.calculation = "scf"
    aIn.WriteToFile(stNew,"alat")
    return 0

def f_QE_ReadBand(stFileName,bLast=False):
    '''
    Read first or last band sturcture from output file
    This is a light-weight version of QESP_ReadBand
    '''
    fIn = open(stFileName)
    list_stLine = fIn.readlines()
    i = 0
    list_kp= []
    list_band = []
    
    step_start = 0
    step_mid = 1
    step_end = 2
    
    step = step_start
    
    if ( bLast ): #Find the last one
        for i in range(len(list_stLine)-1,0,-1):
            if ( "End of self" in stLine or "End of band" in stLine):
                break
        if ( i == 0 ):
            raise ValueError,"Cannot find band information in file"
        i -= 2

    while True:
        stLine = list_stLine[i]
        if ( step == step_start ):
            if ( not "End of " in stLine ):
                i += 1
                continue
            step = step_mid
            i += 2
        elif ( step == step_mid):
            if ( not "band" in stLine):
                i += 1
                if ( "Writing output" in stLine):
                    break
            else:
                list_kp.append([float(stLine[13:20]),float(stLine[20:27]),float(stLine[27:34])])
                i += 2
                stLine = list_stLine[i]
                list_band.append([])
                while ( len(stLine.strip()) != 0):
                    list_band[-1] += [float(x) for x in stLine.split()]
                    i += 1
                    stLine = list_stLine[i]
                
    fIn.close()
    return list_band,list_kp

def QESP_ReadBand(stFileName):
    '''
    Read the lastest band information from pw.x output file
    @todo spin is always 1
    '''
    qOut = QESP_PW_output(stFileName)
    aKPt = KPointsT()
    aKPt.ReadFromList(qOut.arBandKPt[-1])
    aKPt.stMode = "cart"
    aKPt.latt = qOut.FinalLattice()
    if (aKPt.latt is None):
        aKPt.latt = qOut.FirstLattice
#Convert to crystal
    aKPt.ConvertUnit("crystal",aKPt.latt.PrimitiveCellVector)
    #print(qOut.FinalLattice())
#    dFermi = qOut.fermi

#   return aKPt,qOut.arBand[-1],dFermi,qOut.nElectron,1
    return BandsT(aKPt,qOut.arBand[-1],qOut.vbm,qOut.fermi,qOut.nElectron,qOut.n_spin,qOut.b_metal)

def QESP_ReadPDOS(stPrefix,stFolder=None):
    '''
    Read the PDOS data from projwfc.x result
    :param stPrefix: the prefix of output file
    :param stFolder: the folder contained the results. If not set, use current directory
    '''
    dicL = {"s":0,"p":1,"d":2,"f":3} #Order of l
    dicM = [ [0],[0,1,-1],[0,1,-1,2,-2]] #Order of m in different l, in real spherical
    Regex= re.compile(pattern=re.escape(stPrefix)+"\\.pdos_(.*)") # filename regex
    reAO = re.compile("atm#(.*)\\((.*)\\)_wfc#(.*)\\((.*)\\)") # atm#1(M)_wfc#1(o) regex
    list_name = []
    list_pdos = []
    list_energy = []
    if ( stFolder == None):
        stFolder = os.getcwd()
    for stFileName in os.listdir(stFolder):
        aMatch = Regex.match( stFileName )
        if ( aMatch != None):
            stPara = aMatch.group(1)
            #read file
            fIn = open(stFileName,'r')
            list_stLine = fIn.readlines()
            fIn.close()
            #read name
            if ( stPara == "tot"): #Total DOS
                name = OrbitalNameT()
                name.total = True
                list_name.append(name)
            else: #Projected DOS
                aMatch = reAO.match(stPara)
                if ( aMatch == None):
                    raise NameError,"File name %s seems a part of PDOS but cannot be recognized" % stFileName
                name = OrbitalNameT()
                name.species = aMatch.group(2).strip()
                name.i = int(aMatch.group(1)) 
                name.l = dicL[aMatch.group(4)]
                name.n = int(aMatch.group(3)) #This is index of orbital, not true n. Deal with later

                #read each line
                #
                for i in range(0,name.l*2+1):
                    name2 = copy.copy(name)
                    name2.m = dicM[name.l][i]
                    list_name.append(name2)
            
            #read energy and pdos
            list_pdos_part = zip(*[ [float(y) for y in x.split()] for x in list_stLine[1:]])
            #Get energy
            list_energy = list_pdos_part[0]
            #Add to list
            #LDOS is not used
            if ( stPara == "tot"): #We choose Total DOS (index=2 is summation of projected DOS)
                list_pdos.extend([list_pdos_part[1]])
            else: #choose PDOS
                list_pdos.extend(list_pdos_part[2:])

#Deal with quantum number n
    list_name2 = []
    for name in list_name:
        bFound = False
        for group in list_name2:
            for name2 in group:
                if ( cmp(
                    [name.i,name.l,name.m,name.z,name.spin],
                    [name2.i,name2.l,name2.m,name2.z,name2.spin],
                    )==0):
                    group.append(name)
                    break
            if ( bFound):
                break
        if ( not bFound):
            list_name2.append([name])
#Deal with multiple n orbital
    for group in list_name2:
        for i,name in enumerate(group):
            name.n = i+1

    return list_name,list_energy,list_pdos,None


class QESP_plotband_input(object):
    '''
    Input file of plotband.x
    '''
    def __init__(self,stFileName=None):
        self.ReadFromFile(os.path.join(os.path.dirname(__file__),"template-qe/plotband.in"))
        if ( stFileName!= None):
            self.ReadFromFile(stFileName)
    
    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        list_stLine = [x.strip() for x in fIn.readlines()]
        fIn.close()
        self.stFile = list_stLine[0]
        (self.fMinEnergy,self.fMaxEnergy) = [float(x) for x in list_stLine[1].split()]
        self.stXM = list_stLine[2]
        self.stPS = list_stLine[3]
        self.fFermi = float(list_stLine[4])
        (self.fInterval,self.fFermi2) = [float(x) for x in list_stLine[5].split()]
        if ( self.fFermi != self.fFermi2):
            print("\33[36mWarning: Fermi level is not consistent in input file! \33[m")
    
    def WriteToFile(self,stFileName):
        fOut = open(stFileName,'w')
        fOut.write(self.stFile)
        fOut.write("\n")
        fOut.write("%4.2f %4.2f\n" % ( self.fMinEnergy,self.fMaxEnergy))
        fOut.write(self.stXM)
        fOut.write("\n")
        fOut.write(self.stPS)
        fOut.write("\n")
        fOut.write("%4.2f\n" % self.fFermi)
        fOut.write("%4.2f %4.2f\n" % ( self.fInterval,self.fFermi))


def f_QE_ReadFermi(stFileName=None):
    '''
    Read fermi or VBM energy from pq_dos.fermi file under current directory or specific pw.x output
    :param stFileName: if given, read from this file as pw.x output
    '''
    stFermiFile="pq_dos.fermi"
    fFermi1 = None
    if ( stFileName == None):
        if ( os.path.exists(stFermiFile)):
            fIn = open(stFermiFile)
            fFermi1 = float(fIn.readline().strip())
            fIn.close()
            print("Fermi level read from pq_dos.fermi: %f" % fFermi1)
    else:
        qOut = QESP_PW_output(stFileName)
        fFermi1 = qOut.fermi

    return fFermi1

def f_QE_ExtendBand(qIn,stFileNameOut=None,fExtendBand=-1):
    '''
    Extend number of bands for DOS/Band calculation
    :param qIn: QESP_PW_input which will be modified
    :param stFileNameOut: pw.x output filename  which is used as reference
    :param fExtendBand: at least fExtendBand * occupied states will be used 
    '''
    if ( stFileNameOut != None):
        qOut = QESP_PW_output(stFileNameOut)        
        if ( fExtendBand != -1):
            dPercent = fExtendBand - 1
            nbnd = int(qOut.nElectron / 2 * dPercent + 1)
            if ( nbnd < 4):
                nbnd= 4
            nbnd += qOut.nElectron / 2
            print "the number of electrons=",qOut.nElectron
            #if ( nbnd > qOut.nElectron):
            #    nbnd = qOut.nElectron
            qIn.AddValue("system", "nbnd", nbnd)
            print("Number of bands in calculation is extended to %d (%3.0f%% of occupied states)" % (nbnd,nbnd*200.0/qOut.nElectron))
        else:            
            print("Use number of bands in scf calculation: %d (%3.0f%% of occupied states)" % (qOut.nKSStates,qOut.nKSStates*200.0/qOut.nElectron) )
            if ( qOut.nKSStates <= qOut.nElectron / 2 ):
                print("Warning: condunction band / unoccupied states are not included in this calculation!")            
    else:
        print("Warning: no pw.x output is specified. Use default number of bands, this may cause conduction bands lost!")
    
    return qIn


dicBandKPt = {
              'simple_cubic':[
                    ['R',[0.5,0.5,0.5],10],
                    ['lambda',[0.25,0.25,0.25],10],
                    ['Gamma',[0,0,0],5],
                    ['Delta',[0.25,0,0],5],
                    ['X',[0.5,0,0],5],
                    ['Z',[0.5,0.25,0],5],
                    ['M',[0.5,0.5,0],5],
                    ['Sigma',[0.25,0.25,0],5],
                    ['Gamma',[0,0,0],5]
               ],
              "fcc":[ 
                     ['W',[1,0.5,0],20],
                     ['L',[0.5,0.5,0.5],10],
                     ['Lambda',[0.25,0.25,0.25],10],
                     ['Gamma',[0,0,0],20],
                     ['Delta',[0.5,0,0],20],
                     ['X',[1,0,0],10],
                     ['Z',[1,0.25,0],10],
                     ['W',[1,0.5,0],10],
                     ['Kappa',[0.75,0.75,0],10]
                     ],
              "hcp":[
                     ['Gamma',[0,0,0],10],
                     ['Sigma',[0.25,0,0],10],
                     ['M',[0.5,0,0],10],
                     ['K',[1.0/3.0,1.0/3.0,0],10],
                     ['lambda',[1.0/6.0,1.0/6.0,0],10],
                     ['Gamma',[0,0,0],5],
                     ['Delta',[0,0,0.25],5],
                     ['A',[0,0,0.5],5]
                     ],
              "bcc":[
                     ['Gamma',[0,0,0],20],
                     ['Delta',[0.5,0,0],20],
                     ['H',[1,0,0],20],
                     ['N',[0.5,0.5,0],10],
                     ['Sigma',[0.25,0.25,0],10],
                     ['Gamma',[0,0,0],10],
                     ['Lambda',[0.25,0.25,0.25],10],
                     ['P',[0.5,0.5,0.5],10]
                     ]
              }

def f_QE_GetSpecialKPt(stLattice):
    '''
    Get special k-points list of a lattice name ( simple_cubic,fcc,bcc or hcp). If not recognized, use simple_cubic.
    '''
    if ( not stLattice in dicBandKPt.keys()):
        print("Warning: Unknown lattice type for band plot, simple_cubic will be used.")
        stLattice = "simple_cubic"
    return dicBandKPt[stLattice]

def f_QE_DetectSpeicalKPt(ibrav):
    '''
    Determine the special k-points in band structure plot
    '''
    dic_kpt = {2:"fcc",3:"bcc",4:"hcp"}
    st_kpSetUse = "simple_cubic"
    if ( dic_kpt.has_key(ibrav)):
        st_kpSetUse = dic_kpt[ibrav] 
    print("Auto-detected lattice symmetry:%s" % st_kpSetUse)
    list_kp = dicBandKPt[st_kpSetUse]
    return  list_kp

def f_QE_ReadAverage(stAvgFile):
    '''
    Read data from average.x output
    Note:for convinience, an additional line will be added as last line axis value + interval,while value is first line ( due to the system is periodic )
    '''
    f = open(stAvgFile)
    list_stLine = f.readlines()
    f.close()
    i = 0 
    while i < len(list_stLine):
        if ( "0.0000" in list_stLine[i]):
            break
        i += 1

    if ( i > 50 ):# not found
        raise ValueError,"Charge at 0.0000 position not found!"
    data = [] 
    while ( i < len(list_stLine)):
        data.append(  [ float(x) for x in list_stLine[i].split() ]  )
        i += 1
    #data = [ [ float(x) for x in y.split()] for y in list_stLine[i:]] # too slow !
    data.append(copy.deepcopy(data[0]))
    data[-1][0] = data[1][0] + data[-2][0]


    #rewrite to verify
    #f = open("r.txt","w")
    #for line in data:
    #    f.write("%14.9f,%14.9f,%14.9f\n" % tuple(line))
    #f.close() 

    return data

def f_QE_GetFileType(stFileName):
    '''
    Get file type ( input or output, of which Quantum-Espresso program
    '''
    f = open(stFileName,'r')
    bIn = False
    while ( True ):
        stLine = f.readline()
        if ( "&" in stLine ):
            bIn = True
            break
        elif ( "starts on" in stLine):
            bIn = False
            break
    f.close()
    
    if ( bIn ):
        return "QESP_pw_input"
    else:
        return "QESP_pw_output"

def f_QE_ReadLattice(stStructFile):
    '''
    Read struct from Quantum-Espresso input or output file ( only the lastest structure is used )
    Automatically detect file type
    :param stStructFile: the filename
    '''
#detect file type
    stType = f_QE_GetFileType(stStructFile)
    if ( stType == "QESP_pw_input" ):
        return QESP_PW_input(stStructFile).GetCell()
    else:
        return QESP_PW_output(stStructFile).arLattice[-1]

class PPHeader(object):
    def __init__(self):
        pass

class PPNonlocal(object):
    '''
    The index of beta functions starts from 1
    '''
    prefix_func = "PP_BETA"
    def __init__(self):
        pass

    def get_beta(self):
        list1 = []
        a = dir(self)
        n = 0
        for key in a:
            val = getattr(self,key)
            if (val is not None):
                if (PPNonlocal.prefix_func in key):
                    list1.append(val)
        self.list_beta = list1
        return list1

class PPMesh(object):
    def __init__(self):
        pass

class PsWaveFuncNL(object):
    def __init__(self):
        pass

class PsWaveFunc(object):
    def __init__(self):
        pass

class PPPswfc(object):
    def __init__(self):
        pass

class UPFv2(object):
    '''
    This class represents a upf file, format v2 in QE
    '''
    def __init__(self,filename=None):
        if (filename is not None):
            self.read(filename)

    @classmethod
    def load_upf(self,filename):
        '''
        Read from a upf file
        '''
        xs = XMLSerilizer(globals(),"upf.xsd")
        obj = xs.Deserilize(filename)
#Read additional informations
        obj.__post_xml__()

        return obj

    def __post_xml__(self):
        '''
        Parse the informations from string to arrays
        '''
        self.PP_MESH.PP_RAB = [float(x) for x in self.PP_MESH.PP_RAB.split()]
        self.PP_MESH.PP_R = [float(x) for x in self.PP_MESH.PP_R.split()]
        self.PP_NONLOCAL.PP_DIJ = [float(x) for x in self.PP_NONLOCAL.PP_DIJ.split()]
        for nl in self.PP_NONLOCAL.get_beta():
            nl.__base__ = [float(x) for x in nl.__base__.split()]
        self.PP_NONLOCAL.get_beta()

def qesp_read_latt(filename="data-file.xml"):
    '''
    Read Quantum Espresso cell from data-file.xml
    '''
    for event, elem in ET.iterparse(filename,events=('end',)):
        if (elem.tag == "CELL"):
            node_cell = elem
        elif (elem.tag == "IONS"):
            node_ions = elem
            break

    node_vec = node_cell.find("DIRECT_LATTICE_VECTORS")
    m_vec = [[float(x) for x in node_vec.find(key).text.split()] for key in ["a1", "a2", "a3"]]
    
    list_species = [x.text for x in node_ions.findall('././ATOM_TYPE')]
    list_atoms = [[x.attrib["SPECIES"].strip()] + [float(y) for y in x.attrib["tau"].split()] for
            x in node_ions if "ATOM." in x.tag]

    unit = node_ions.find("./UNITS_FOR_ATOMIC_POSITIONS").attrib["UNITS"]

    latt = Lattice()
    latt.ReadFromRawCellVector(m_vec)
    if (unit == "Bohr"):
        latt.AddAtom(list_atoms,unit="bohr")
    else:
        raise ValueError("Unknown unit %s" % unit)

    return latt


def qesp_read_band(filename="data-file.xml"):
    '''
    Read Quantum Espresso band
    if the filename ends with xml, read as format data-file.xml
        else use QESP_ReadBand to read pw output
    '''
    if (not filename.endswith("xml")):
        return QESP_ReadBand(filename)

    dir_rel = os.path.dirname(filename)

    latt = qesp_read_latt(filename)

#Read cell and eigenvalues
    for event, elem in ET.iterparse(filename,events=('end',)):
        if (elem.tag == "CELL"):
            node_cell = elem
        elif (elem.tag == "BAND_STRUCTURE_INFO"):
            node_info = elem
        elif (elem.tag == "EIGENVALUES"):
            node_eigenvalues = elem
            break
    
    d_electron = float(node_info.find("./NUMBER_OF_ELECTRONS").text)
    n_electron = int(d_electron)
    b_electron_fraction = abs(d_electron != n_electron) > 1e-7

    num_spin = int(node_info.find("./NUMBER_OF_SPIN_COMPONENTS").text)
    eig_fermi = float(node_info.find("./FERMI_ENERGY").text) * Ha2eV

    list_kpt = []
    list_band = [[], []]
    list_occ = [[], []]
    for kpt_data in node_eigenvalues:
        list_kpt.append([float(x) for x in (
            kpt_data.find("./K-POINT_COORDS").text.split() + [kpt_data.find("./WEIGHT").text])])
        for ix_spin, child in enumerate(filter(lambda x:"DATAFILE" in x.tag, kpt_data)):
            node2 = ET.parse(os.path.join(dir_rel, child.attrib["iotk_link"]))
            unit_name = node2.find("./UNITS_FOR_ENERGIES").attrib["UNITS"]
            if (unit_name == "Hartree"):
                unit = Ha2eV
            else:
                raise ValueError("Unknown unit %s" % unit_name)
            list_band[ix_spin].append([float(x) * unit for x in node2.find("EIGENVALUES").text.split()])
            list_occ[ix_spin].append([float(x) for x in node2.find("OCCUPATIONS").text.split()])

    if (num_spin == 1):
        list_band = list_band[0]
        list_occ = list_occ[0]


    aKPT = KPointsT()
    aKPT.ReadFromList(list_kpt)
    aKPT.stMode = {"2 pi / a":"tpiba"}[node_info.find("UNITS_FOR_K-POINTS").attrib["UNITS"]]
    aKPT.latt = latt
    aKPT.ConvertUnit("cart")

    band = BandsT(aKPT,list_band,None,None,n_electron,num_spin, list_occ=list_occ)

    #Calculate VBM
#Only do if the number of electrons is an integer
    if (b_electron_fraction):
        band.fermi = eig_fermi
    else:
        if (not band.guess_vbm()):
            band.fermi = eig_fermi

    return band


def qesp_read_vacuum_level(filename, pos_vac):
    '''
    Read the vacuum level in eV from average.x output

    :param pos_vac: the position of the middle of vacuum, in crystal coordinates
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if ("Reading data" in line):
            ix_start = i
        elif ("AVERAGE  " in line):
            ix_end = i
    ix = int(ix_start + (ix_end - ix_start) * pos_vac)
    return float(lines[ix].split()[1]) * Ry2eV

