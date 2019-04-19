#!/usr/bin/env

from common_caseutil import KPointsT
from io_utils import io_grep_lines
from constants import Ha2eV
from band_utils import BandsT
import re,copy

class abi_BandStructure:
    '''
    Class to represent information of one eigenvalue
    Mostly stored GW corrected information
    '''

    def __init__(self):
        self.colName = [] #Column name
        self.listKPt = [] #k-points coordinate, unit unknown
        self.listBand = [] #Band information of each eigenvalue*kpoints. Unlike BandsT, each element is a list, stored eigenvalues and corrections
        self.dFermi = None #Fermi energy

#For iterator
        self.ix_n = 0
        self.ix_k = 0

    def __getitem__(self, x):
        '''
        Return specific columns of eigenvalue information
        :param x: tuple of the index of k-points, the band index and the information name. For perturbation GW, use "E(N)_pert"; for DFT eigenvalue, use "E_lda"
        '''
        (k, n, name) = x
        if (not name in self.colName):
            raise ValueError("Specified name %s does not exist in band structure information")
        i = self.colName.index(name)
        return self.listBand[k][n][i]

    def next(self):
        if (len(self.listBand) == 0):
            raise StopIteration()
        if (len(self.listBand[0]) == 0):
            raise StopIteration()
        if (self.ix_k == len(self.listBand)):
            self.ix_k = 0
            self.ix_n = 0
            raise StopIteration()
        val = self.listBand[self.ix_k][self.ix_n]
        if (self.ix_n == len(self.listBand[self.ix_k])-1):
            self.ix_k += 1
            self.ix_n = 0
        else:
            self.ix_n += 1
        return val
    
    def __iter__(self):
        return self

    def get_list_of_name(self,name):
        '''
        Get a 2-D list for specific column name
        '''
        if (not name in self.colName):
            raise IndexError("Name %s not found!" % name)
        ix = self.colName.index(name)
        return [x[ix] for x in self]


    def ReadFromFile(self, f, stLine1=None):
        '''
        Read band structure information from abinit output file
        :param f: the file handler to be read from
        :param stLine1: the preceding line before file pointer; This is useful if first line we need is already skipped
        '''
#Read first line
        if (stLine1 != None):
            stLine = stLine1
        else:
            stLine = f.readline()

        aKPt = None
        aBand = None
#Read file
        while stLine != "":
            if ("New Fermi" in stLine): #Read Fermi energy and break
                self.dFermi = float(stLine.split()[-2])
                break
            if (stLine.startswith("==")): #End of a dataset without Fermi level inside
                break

            if (aKPt != None):
                if (stLine.strip() == "" or stLine[1] != " "):#End of eigenvalues in one k-point, empty lines or non-number lines
                    self.listKPt.append(aKPt)
                    self.listBand.append(aBand)
                    #Clean temporary variables
                    aKPt = None
                    aBand = None
                else: #Read eigenvalues
#Skip unused second line in AC method, in case where 2-4 columns are all zero
                   # l1 = [float(x) for x in stLine.split()]
                    l1 = stLine.split()
#This line is dangerous as some numbers may collide together..
                    b_err = 0
                    for st in l1:
                        if (len(st) > 8):#Maybe we should reread
                            b_err = 1
                            break
                    if (b_err==1):#This error maybe occured only in AC?5+8*9
                        l1 = [float(stLine[:5])]
                        for i in range(9):
                            l1.append(float(stLine[5+i*8:5+i*8+8]))
                    else:
                        l1 = [float(x) for x in l1]

                    if (l1[3] == 0.0 and l1[1] == 0.0 and l1[2] == 0.0):
                        pass
                    else:
                        aBand.append(l1)

            if (" k =" in stLine): #Read k-points
                aKPt = [float(x) for x in stLine.split()[-3:]]
                aBand = []
                if (len(self.colName) == 0): #If column names are not read, then read it
                    self.colName = f.readline().strip().split()
                else: #Just skip one line
                    f.readline()
            stLine = f.readline()


def abi_ReadBand_Output_GW(stFileName):
    '''
    Read GW band structure from output file
    :return: a list of abi_BandStructure object, one object for one DATASET
    '''
    f = open(stFileName)
    listResult = []
    stLine = " "

    while stLine != "": #Read to end of file
        stLine = f.readline()
        if (stLine.startswith(" k =")):
            ab = abi_BandStructure()
            ab.ReadFromFile(f, stLine1=stLine)
            listResult.append(ab)

    f.close()

    return listResult


def abi_ReadBand_GW(stFileName):
    '''
    Read band structure from _GW file
    '''
    f = open(stFileName)
    nKPt, nSpin = [int(x) for x in f.readline().split()]
    listKPt = []
    listBand = []
    for i in xrange(0, nKPt):
        listKPt.append([float(x) for x in f.readline().split()])
        n = int(f.readline())
        listB = []
        for j in xrange(0, n):
            listB.append(float(f.readline().split()[1]))
        listBand.append(listB)

    aKPt = KPointsT()
    aKPt.ReadFromList(listKPt)

#   return aKPt, listBand, None, None, nSpin
    return BandsT(aKPt,listBand,None,None,None,nSpin)


def abi_ReadBand_EIG(stFileName):
    '''
    Read band structure from _EIG file
    '''
    f = open(stFileName)
    stLine = f.readline()
    if ("hartree" in stLine): # "(hartree) ="
        energy_unit = Ha2eV
    else:
        energy_unit = 1.0 # energy unit; if it is hartree then convert to eV
    if ("Fermi" in stLine):# Find the first equal
        ixEqual = stLine.index("=")
        dFermi = float(stLine[ixEqual + 1:ixEqual + 10])
        stLine = f.readline()
    else:
        dFermi = None
    nKPt = int(stLine[32:36])
    listKPt = []
    listBand = []
    listB1 = []
#    for i in xrange(0,nKPt):
    while True:
        stLine = f.readline()
        if not stLine:
            break
        stLine = stLine[:-1] #Remove \cr for splitting
        if ("kpt" in  stLine):
            listKPt.append([float(stLine[41:49]), float(stLine[49:57]), float(stLine[58:65]), float(stLine[26:35])])
            if (len(listB1) > 0):
                listBand.append(listB1)
            listB1 = []
        else:
#To avoid two maximum length float recognized as one
#           listB1 += [float(x) if not "*" in x else 0.0 * energy_unit for x in stLine.split()]
            n = 10
            listB1 += [float(x) if not "*" in x else 0.0 for x in [stLine[i:i+n] for i in xrange(0,len(stLine),n)]]
    #Find band
    listBand.append(listB1)

    aKPt = KPointsT()
    aKPt.ReadFromList(listKPt)

#   return aKPt, listBand, dFermi, None, 1
    return BandsT(aKPt,listBand,dFermi,None,None,1)


def abi_GetElectron(stFileName):
    '''
    Get electron count from log file
    '''
    stEle = io_grep_lines(stFileName, "nelect", 1, 2)
    if (stEle == ""):
        raise ValueError("Cannot find nelect in log file %s" % stFileName)
    dEle = float(stEle)
    if (abs(int(dEle) - dEle) < 0.01):
        return int(dEle)
    else:
        raise ValueError("Electron number %f is not an integer!" % dEle)

def abi_ReadFermi(filename):
    '''
    Get the last fermi energy from log file in unit eV
    '''
    line = io_grep_lines(filename, "Fermi energy", -1, 7)
    if (line == ""):
        raise ValueError("Cannot find nelect in log file %s" % stFileName)
    else:
        return float(line)*Ha2eV

def abi_ReadBand(stFileName, stLog):
    '''
    Read band structure from _GW or _EIG files, automatically detected by filename
    :param stFileName: the filename, must be end with _GW or _EIG, which is used to detected file type
    :param stLog: the log file of the calculation. The sole purpose is to read number of electrons.
    '''
    if ("_GW" in stFileName):
#       aKPt, listBand, dFermi, nElectron, nSpin = abi_ReadBand_GW(stFileName)
        band = abi_ReadBand_GW(stFileName)
    elif ("_EIG" in stFileName):
#       aKPt, listBand, dFermi, nElectron, nSpin = abi_ReadBand_EIG(stFileName)
        band = abi_ReadBand_EIG(stFileName)
    else:
        raise ValueError("Filename %s is not a valid abinit eigenvalue output file" % stFileName)
#   nElectron = abi_GetElectron(stLog)
#   return aKPt, listBand, dFermi, nElectron, nSpin
    band.num_electron = abi_GetElectron(stLog)
    if (band.vbm is None):
        band.vbm = abi_ReadFermi(stLog)
    vbm_old = band.vbm
    band.vbm = None
    if (band.guess_vbm() == False):
        band.fermi = vbm_old
    return band

class abi_input_var():
    '''
    A class represents a variable in an ABINIT input file
    maybe n*x or array
    '''
    tab_name = 16 #The length of variable name used in output

    def __init__(self,f=None):
        self.name = ""
        self.set1 = " "
        self.set2 = " "
        self.value = []
        if (f is not None):
            self.read(f)

    def __read_nocomment__(self,f):
        '''
        Read the file and return the first nonempty line (stripped) without comments
        return an empty string if EOF reached
        '''
        line = f.readline()
        while (line != ""):
            ix_comment = line.find("#")
            if ( ix_comment != -1):
                line = line[:ix_comment]
            line = line.strip()
            if (line != ""):
                return line
            line = f.readline()
        return ""#Indicate 

    def read(self,f):
        '''
        Read a variable from an input file object
        The pointer in the file is shifted back to the end of this variable,
        which makes it possible to read next one
        :return: -1 for EOF reached and 0 else
        '''
        while (True):
            ix_last = f.tell()
            line = self.__read_nocomment__(f)
            if (line == ""):
                f.seek(ix_last)
                return -1
#Get the name
            if (self.name == ""):
                ix = line.find(" ")
                if (ix != -1):
                    self.name = line[:ix]
                    line = line[ix+1:].strip()
                else:
                    self.name = line
                    line = ""#This empty makes values and name in different lines
#Get the No of set
                i = re.search(r"[\d\?\+:]",self.name)
                if (i is not None):
                    i = i.start()
                    if (i == len(self.name)-1): #single loop
                        self.set1  = self.name[-1]
                        self.name  = self.name[:-1]
                    elif (i == len(self.name)-2):
                        self.set1  = self.name[-2]
                        self.set2  = self.name[-1]
                        self.name  = self.name[:-2]
            else:
                if (not line[0].isdigit()):#New variable
                    f.seek(ix_last)
                    return 0
            self.value.append(line)
        pass

    def write(self,f):
        '''
        Write this variable to a new file
        '''
#Align dataset number
        st_format = "%%%is%%1s%%1s  %%s\n" % (abi_input_var.tab_name)
        f.write(st_format % (self.name,self.set1,self.set2,self.value[0]))
        for line in self.value[1:]:
            f.write(st_format % (" "," "," ",line))


class abi_input():
    '''
    A class represents a abinit input file (version 7)
    Some conventions are assumed
    '''
    dic_group = {} #Group defined in ABINIT documentation
    filename_vardef = "template-abinit/var.def"

    @classmethod
    def __read_groupinfo__(cls):
        '''
        Read ABINIT variables groups
        '''
        f = open(os.path.join(os.path.dirname(__file__),abi_input.filename_vardef))
        lines = f.readlines()
        f.close()
        for i in xrange(len(lines-1)/2):
            vars = lines[i*2+2].strip().split()
            group = lines[i*2+1].strip()
            for var in vars:
                abi_input.dic_group[var] = group

        return

    def __init__(self,filename=None):
        self.list_var = []
        self.dic_var = {}
        if (filename is not None):
            self.read(filename)

    def read(self,filename):
        '''
        Read ab ABINIT input file
        '''
        f = open(filename,'r')
        while(True):
            var1 = abi_input_var()
            err = var1.read(f)
            if (err == -1):
                break
            self.list_var.append(var1)
            self.dic_var[var1.name] = var1
        f.close()
        return 

    def write_simple(self,filename):
        '''
        Write ABINIT input file
        Variables are written one by one
        '''
        f = open(filename,'w')
        for var in self.list_var:
            var.write(f)
        f.close()
        return

    def write_group(self,filename):
        '''
        Write ABINIT input file
        Variables are written grouped and automatically commented
        Rules: Xij will be put together with same j and then same i
        Xij are always before Xi and X
        '''

        f = open(filename,'w')
        
        def order_set(s1,s2):
            if ( s1.isdigit() and s2.isdigit() ): #Digit first
                return int(s1)-int(s2)
            elif (s1.isdigit()):
                return -1
            elif (s2.isdigit()):
                return 1
            elif (s1 != " "): #Other :,?,+ first
                return -1
            elif (s2 != " "):
                return 1

            return 0

        def order_var(v1,v2):
            i = order_set(v1.set2,v2.set2)
            if (i == 0):
                i = order_set(v1.set1,v2.set1)
            return i

        list_var2 = copy.copy(self.list_var)
        #Write runtime information first (ndtset,udtset and varpa)
        set_pass = set()

        f.write("# General information\n")
        for i,var in enumerate(list_var2):
            group = abi_input.dic_group[var]
            if (var.name == "ndtset" or var.name == "udtset" or group == "varpar"):
                var.write(f)
                set_pass.add(i)
        f.write("\n\n")

        list_var2 = [x for i,x in enumerate(list_var2)]
#Write loop split information
#Second loop first
        list_var2.sort(order_var)
        if (list_var2[0].set2 != " "):
            f.write("# Loop part\n")

            f.write("\n\n")

        b_last2 = None
        b_last1 = None

#Common part
        f.write("# Common part\n")
        for var in list2:
            var.write(f)
        f.close()
        f.write("\n")
        
        return
        

def abi_input2wannier(filename):
    '''
    Create wannier90 calculation input file based on original abinit input file
    '''
    pass
