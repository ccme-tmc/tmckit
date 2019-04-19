#!/usr/bin/env python
#yaehmop package utility

import os,copy 
import operator
import itertools
import list_utils as lu
from common_caseutil import deprecated,KPointsT,Lattice 
from constants import Bohr2Ang

from band_utils import f_Band_ReadFromFile,f_Band_GetVBMFromElectron,f_Band_GetGap,BandsT

class yaeh_input():
    '''
    Input object of yaehmop
    Note this is also used as a library of parameters
    '''
    def __init__(self,stFileName,b_write_fix=False):
        '''
        :param b_write_fix: Write "fixed" informations to the input file. YAehmop does not support this.
        '''
        self.listPara = [] # EHT parameters
        self.stPart1 = [] #Content before EHT parameters
        self.stPart2 = [] #Content after EHT parameters
        self.ReadFromFile(stFileName)
#        self.listFix = [] #Elements names that are fixed and will not change during read/write
#        self.listFixDZ = [] #Elements with second-zeta of d fixed and will not chagne during read/write
        self.b_write_fix = b_write_fix
        

    def ReadFromFile(self,stFileName):
        f = open(stFileName)
        ar_stLine = f.readlines()
        f.close()

        #Determine parameters start position
        for i,stLine in enumerate(ar_stLine):
            if ( len(stLine.strip()) == 0):
                continue
            if ( stLine.strip()[0] == ";"):
                continue
            if ( "parameters" in stLine):
                break
        self.stPart1 = ar_stLine[:i+1]
        nStart = i+1
        self.listPara = []
#read split string
        for i in range(nStart,len(ar_stLine)):
            stLine = ar_stLine[i]
            if ( len(stLine.strip()) == 0):
                break
            self.listPara.append(yaeh_para(stLine))

        #If there still lines then add it
        if ( i < len(ar_stLine)-1):
            self.stPart2 = ar_stLine[i:]
        #print(listPara)

    def WriteToFile(self,stFileName):
        f = open(stFileName,'w')
        for stLine in self.stPart1:
            f.write(stLine)
        for para in self.listPara:
            #f.write(" ".join([str(x) for x in para]))
            #print(para.Write())
            f.write(para.Write(b_write_fix=self.b_write_fix))
            f.write("\n")
        for stLine in self.stPart2:
            f.write(stLine)        

    def write(self,stFileName):
        self.WriteToFile(stFileName)

    def GetPara(self):
        '''
        Return all variable orbtial parameters
        '''
        list2 = []
        for para in self.listPara:
            list2 = list2 + para.Save()
        return list2

    def GetParamType(self):
        '''
        Return all variable type (nuclear charge, zeta or energy)
        '''
        a = [x.GetParamType() for x in self.listPara]
        return list(itertools.chain(*[x.GetParamType() for x in self.listPara]))


    def get_para(self):
        return self.GetPara()

    def __getitem__(self,x):
        '''
        Return parameters of one specific element
        :return: yaeh_para class instance. If not found, then return None
        '''
        for para in self.listPara:
            if ( para.stElement == x):
                return para
#Not found, return None
        return None

    def __setitem__(self,key,x):
        '''
        Set parameters of on specific element
        If not exist then append
        '''
        nIndex = -1
        for i,para in enumerate(self.listPara):
            if ( para.stElement == key):
                nIndex = i
                break
        if ( nIndex == -1):#Append new parameters
            self.listPara.append(x)
        else:
            self.listPara[nIndex] = x

    def set_para(self,para):
        '''
        Set all variable parameters in the input
        '''
        self.ReloadPara(para)

    def load_para(self,data):
        '''
        Reload all orbital parameters

        :param data: a yaeh_input instance with all necessary parameters
        '''
        self.ReloadParaFromBase(data.listPara)

    def ReloadPara(self,listNew):
        '''
        Reload all orbital parameters in a input file from an array
        '''
        #print("New para length %i" % len(listNew))
        nStart = 0
        for para in self.listPara:
            #print("Para count: %i" % para.nVarCount)
            para.Load(listNew[nStart:nStart+para.nVarCount])
            nStart += para.nVarCount

    def ReloadParaFromBase(self,listNew):
        '''
        Reload all orbital parameters from parameter database
        '''
        for i in range(0,len(self.listPara)):
            bFound = False
            for para2 in listNew:
                if ( para2.stElement == self.listPara[i].stElement):
                    self.listPara[i] = copy.deepcopy(para2)
                    bFound = True
                    break
            if ( not bFound):
                raise ValueError,"Element %s not found in parameter database" % self.listPara[i].stElement

class yaeh_para:
    '''
    Extended Huckel Method parameters of one elements
    f orbital is not support now
    for spd:
    Ti 22 4 4   3.80905906  -8.73503789 4  2.8544081   -7.13526927 3  2.75288322 -14.9012884   -2.80545149   2.66434072   2.99185815 
    for spd + charge: 
    Ti 22 4 nucval 2.0 3.0 4   3.80905906  -8.73503789 4  2.8544081   -7.13526927 3  2.75288322 -14.9012884   -2.80545149   2.66434072   2.99185815 
    '''
    key_fix = "fix"
    num_param_no_d = 4 #: Number of variable parameters without d orbitals
    num_param_with_d = 9 #: Number of variable parameters with d orbitals
    num_param_charge = 2 #: Number of variable parameters related to charges
    num_param_max = 11 #: Number of all possible parameters
    num_item_charge = 3 #: Number of items in line 
    mode_charge_ratio = "nucratio"  #: two numbers are charge and ratio of cf/3c
    mode_charge_val = "nucval" #two number are cf charge and 3c charge
    param_type_nuclear = 1 #: Type of a given parameter is a nuclear charge
    param_type_energy = 2 #: Type of a given parameter is an energy
    param_type_exponent = 3 #: Type of a given parameter is an exponent
    param_type_coef = 4 #: Type of a given parameter is a coefficient of a contracted orbital

    def __init__(self,stInput=None):
        self.stElement=""
        self.nZ = 0
        self.nVal = 0
        self.s_n = 0
        self.s_zeta = 0
        self.s_energy = 0
        self.p_n = 0
        self.p_zeta = 0
        self.p_energy = 0
#if d_n ==0, then it is ignored in output
        self.d_n = 0
        self.d_zeta_1 = 0
        self.d_energy = 0
        self.d_coef_1 = 0
        self.d_zeta_2 = 0
        self.d_coef_2 = 0
#Additional part for nuclear attaction
#Marked by charge_mode
#       self.charge_pseudo_nuc = 0.0
#       self.charge_crystal_field = 0.0
#       self.charge_ratio_nuc = 0.0
        self.charge_mode = None
        self.charge1 = 0.0
        self.charge2 = 0.0

#Fix or not
#This array must be as long as maximum number of parameters
        self.listFixed = [False for x in xrange(yaeh_para.num_param_max)]  #Fixed parameters. d-orbital parameters are always fixed if d_n =0

        if ( stInput != None):
            self.Read(stInput)

    def Read(self,st):
        '''
        Read parameter from one line in input of YAehmop
        Anything after "fix" is treated as fit information
        '''
        ar = st.split()
        ix_fix = ar.index(yaeh_para.key_fix) if (yaeh_para.key_fix in ar) else -1
#How to detect if there is charge?
#1. number of fixed options
#2. whether float or int for the third atom
#Cross-checked
        b_contains_charge = "nuc" in ar[3]

        if (ix_fix != -1):
            ar2 = [int(x) for x in ar[ix_fix+1:]]
            ar = ar[:ix_fix]
            l1 = len(ar)
            l2 = len(ar2)
            if (l2 == yaeh_para.num_param_no_d or l2 == yaeh_para.num_param_with_d):
                b2 = False
            elif (l2 == yaeh_para.num_param_no_d + yaeh_para.num_param_charge or
                    l2 == yaeh_para.num_param_with_d + yaeh_para.num_param_charge):
                b2 = True
            else:
                raise ValueError("Number of options fixed must be 4 or 9, 6 or 11 (with charges)")
            if (b2 != b_contains_charge):
                raise ValueError("Charges are included only in one of values / fixed (%s / %s)" % (b2, b_contains_charge))

        self.stElement = ar[0]
        self.nZ = int(ar[1])
        self.nVal = int(ar[2])
#Charge
        if (b_contains_charge):
            n_item_charge = yaeh_para.num_item_charge
            n_fix_charge = yaeh_para.num_param_charge
            self.charge_mode = ar[3]
            self.charge1, self.charge2 = [float(x) for x in ar[4:6]]
        else:
            n_item_charge = 0
            n_fix_charge = 0
#Orbitals
        (self.s_n,self.s_zeta,self.s_energy,self.p_n,self.p_zeta,self.p_energy) = \
                [float(x) for x in ar[n_item_charge + 3:n_item_charge + 9]]
        self.s_n = int(self.s_n)
        self.p_n = int(self.p_n)
        if ( len(ar) > 9 + n_item_charge):
            (self.d_n,self.d_zeta_1,self.d_energy,self.d_coef_1,self.d_zeta_2,self.d_coef_2) = \
                [float(x) for x in ar[n_item_charge + 9:n_item_charge + 15]]
            self.d_n = int(self.d_n)
        else:#Fix all d orbitals
            self.listFixed[yaeh_para.num_param_no_d + n_fix_charge:yaeh_para.num_param_with_d + n_fix_charge] = \
                    [True for x in xrange(yaeh_para.num_param_with_d-yaeh_para.num_param_no_d)]

        if (ix_fix != -1):
            self.listFixed[:len(ar2)] = [ x==1 for x in ar2]

    def Write(self, b_write_fix = False):
        '''
        Create a line that can be written into YAehmop

        :param b_write_fix: Write informations of fixed
        '''
        st = "%3s %3i %3i " % (self.stElement,self.nZ,self.nVal)
        if (self.charge_mode is not None):
            st += "%s %f %f " % (self.charge_mode, self.charge1, self.charge2)
        
        st += "%i %f %f %i %f %f" % (self.s_n,self.s_zeta,self.s_energy,self.p_n,self.p_zeta,self.p_energy)
        if ( self.d_n > 0):
            st = st + " %i %f %f %f %f %f" % (self.d_n,self.d_zeta_1,self.d_energy,self.d_coef_1,self.d_zeta_2,self.d_coef_2)
        if (b_write_fix):
            n_param = yaeh_para.num_param_with_d if self.d_n > 0 else yaeh_para.num_param_no_d
            if (self.charge_mode is not None):
                n_param += yaeh_para.num_param_charge
            list2 = ["1" if x else "0" for x in self.listFixed]
            st = "%s fix %s" % (st, " ".join([x for x in list2[:n_param]]))
        return st

    @property
    def nVarCount(self):
        '''
        Return number of variable parameters, for sp is 4, for spd is 9
        Fixed is possbile
        '''
        return 9 - self.nFixCount + yaeh_para.num_param_charge if (self.charge_mode is not None) else 0

    @property    
    def nFixCount(self):
        '''
        Return the number of fixed parameters, maximumly 9
        '''
        return sum([1 for x in self.listFixed if x])#Count how many parameters fixed


    def FixAll(self):
        '''
        Set all parameters for this set are fixed
        '''
        self.listFixed = [True for x in xrange(9)]

    def FixDZ(self):
        '''
        Set all parameters for the second zeta of d-orbital useless
        '''
        self.listFixed[6:] = [True for x in xrange(3)]
        self.d_zeta_2 = 10.0
        self.d_coef_1 = 1.0
        self.d_coef_2 = 0.0

    def Load(self,ar):
        '''
        Load variable parameters from an array
        '''
        nSkip = self.nFixCount
        if (len(ar) + nSkip != yaeh_para.num_param_max):
            raise ValueError,"Number of parameters not consistent with orbital count"
        listOld = (self.charge1, self.charge2, 
                self.s_zeta,self.s_energy,self.p_zeta,self.p_energy,
                self.d_zeta_1,self.d_energy,self.d_coef_1,self.d_zeta_2,self.d_coef_2)
        
        listNew = []
        j = 0
        for i in xrange(yaeh_para.num_param_max):
            if (self.listFixed[i]):
#               print("Skip %i" % i)
                listNew.append(listOld[i])
            else:
#               print("Add %i" %i)
                listNew.append(ar[j])
                j += 1

        (self.charge1, self.charge2,
                self.s_zeta,self.s_energy,self.p_zeta,self.p_energy,
                self.d_zeta_1,self.d_energy,self.d_coef_1,self.d_zeta_2,self.d_coef_2) = listNew

    def Save(self):
        '''
        Save variable parameters to an array
        '''
        if (self.d_n == 0):
            result = [self.s_zeta,self.s_energy,self.p_zeta,self.p_energy]
        else:
            result = [self.s_zeta,self.s_energy,self.p_zeta,self.p_energy,self.d_zeta_1,self.d_energy,self.d_coef_1,self.d_zeta_2,self.d_coef_2]
        if (self.charge_mode is not None):
            result = [self.charge1, self.charge2] + result
        #Remove fixed parameters
        return [ x[0] for x in zip(result,self.listFixed) if not x[1]]
        

    def GetParamType(self):
        '''
        Get all types (nuclear charge, energy and exponent) in 
        '''
        if (self.d_n == 0):
            result = [yaeh_para.param_type_exponent, yaeh_para.param_type_energy] * 2
        else:
            result = [yaeh_para.param_type_exponent, yaeh_para.param_type_energy] * 3 + [yaeh_para.param_type_coef, yaeh_para.param_type_exponent, yaeh_para.param_type_energy] 
        if (self.charge_mode is not None):
            result = [yaeh_para.param_type_nuclear] * 2 + result
        #Remove fixed parameters
        return [ x[0] for x in zip(result,self.listFixed) if not x[1]]
    
def yaeh_ReadBand(stFileName):
    '''
    Read k-points and band from yaeh band result.
    Fermi energy is not calculated.
    @todo spin is always 1
    '''
    f = open(stFileName)
    ar_stLine =  f.readlines()

    list_kp = []
    list_band = []

    nbnd = int(ar_stLine[4].split()[0]) #Number of bands in each kpt
    nStart = 0 #Start position of band data
    for i in range(0,len(ar_stLine)):
        if ( "Begin band data" in ar_stLine[i]):
            nStart = i+1 #line count that not band data
            break

    for i in range( 0,(len(ar_stLine)-nStart)/(nbnd+1)):
        i2  = (nbnd+1)*i + nStart
        list_kp.append([float(x) for x in ar_stLine[i2][10:].split()])
        list_band.append([float(x) for x in ar_stLine[i2+1:i2+nbnd+1]])


    #write
    aKPT = KPointsT()
    aKPT.ReadFromList(list_kp)

#   return aKPT,list_band,0,0,0
    band = BandsT(aKPT,list_band)
    band.num_spin = 1
    return band

def yaeh_CreateInputFile(stFileName, latt, kpt=None, with_charge=False):
    '''
    Create geometry part of yaeh input file from Lattice object
    '''
    listAtom = latt.GetAtomList(latt="prim",unit="conv") #primitive lattice
    listAppear = [] #Record appeared atom to use "*" as notation

#Write atom fractional coordinate list
    f = open(stFileName,'w')
    f.write("\nCase Name\n\n; define the geometry, in crystallographic coordinates\nGeometry Crystallographic\n;this is the number of atoms\n%i\n"%(len(listAtom)+3))
    stAtomFormat = "%3i %2s %f %f %f\n" 
    i = 0
    atomStart = listAtom[0]
    for atom in listAtom:
        i += 1
        stElement = atom[0]
        if ( not atom[0] in listAppear):
            listAppear.append(atom[0])
            stElement = "*"
        f.write(stAtomFormat% (i,stElement,atom[1],atom[2],atom[3]))
#Three atoms to indicate cell vectors
#Rule: use primitive vectors here
#Use conventional vectors in crystal description
    vecPrim = lu.f_Matrix_transpose(lu.f_Matrix_dot(lu.f_List_inv3(lu.f_Matrix_transpose(latt.ConventionalCellVector)),lu.f_Matrix_transpose(latt.PrimitiveCellVector)))
    for atom2 in vecPrim:
        i += 1
        f.write(stAtomFormat % (i,listAppear[0],atomStart[1]+atom2[0],atomStart[2]+atom2[1],atomStart[3] + atom2[2]))

    f.write("\nparameters\n")
    for atom in listAppear:
        f.write(atom+" 0 0 0 0.0 0.0 0 0.0 0.0 \n")

    f.write('''
; definition of the lattice
lattice
; dimensionality
3
; number of overlaps in each direction
3 3 3
; definition of a (it goes from 1->3)
1 %i
; definition of b (it goes from 1->4)
1 %i
; definition of c (it goes from 1->5)
1 %i

Crystal Spec
; the lengths of the lattice vectors a, b, and c
''' % ( i-2,i-1,i )  )

#unit angstorm
    f.write("%f %f %f\n" % tuple( [x * Bohr2Ang for x in latt.fLatticeLength] ))
    f.write("; the crystallographic angles alpha, beta, and gamma.\n")
    f.write("%f %f %f\n" % tuple( latt.fLatticeAngle ) )

    f.write("\nNearest Neighbor Contact\n2.8\n\n")

    if ( kpt != None):
        f.write("Band\n1\n%i\n" % len(kpt.listKPt))
        kpt.ConvertUnit("crystal")
        if ( kpt.stMode != "crystal"):
            raise ValueError,"k-points used by YAEHMOP must be in fractional coordinate"
        for a in kpt.listKPt:
            stName = a[0]
            if ( stName ==  ""):
                stName = "0"
            f.write("%7s %11.7f %11.7f %11.7f\n" % ( stName,a[1],a[2],a[3]))
        f.write("\n")

    else:
        f.write("Band\n1\n111\n")
    
    f.write('''
RHO
7

''')
    if (with_charge):
        f.write('''
CUTOFFNUC
2.5

INCLUDE THREE CENTER
INCLUDE CRYSTAL FIELD
''')
    f.write('''
print
energy levels
end_print

;wave func transp
;charge mat transp
;distance matrix
;wave functions
;overlap matrix
;hamiltonian

;average properties

electrons
1

; These k points are not reliable for a real average properties
; calculation, but they allow us to look at the wavefunctions,
; charge matrix, etc. at the high symmetry points.
k points
0
;1
;0 0 0 1

''')

    f.close()

class PEMiniSingleInput:
    '''
    Input of single-compound parameter fit
    This class does not have initial value, so it must be deserilized from xml to use
    '''
    def __init__(self):
#Below part will not be serilized or de-
        ## @var paraBase 
        # Link to database of EH parameters as a list
        self.paraBase = None
        self.yIn = None
        self.bandRef = None
        self.dVBMRef = None
        self.dGapRef = None
        self.listBandWeight = None
#Below part is used to store result
        self.dGap = None
        self.dRMSD = None
        pass
    
    def Load(self):
        '''
        Load external content with given parameters
        '''
        self.bandRef = f_Band_ReadFromFile(self.stBandRefFile)
        self.yIn = yaeh_input(self.stInputFile)
        self.listBandWeight = [0.001 for x in range(0,self.nBandStartIndex)] + [ 1.0 for x in range(self.nBandStartIndex,self.nBandStartIndex+self.nBandFitCount)]
        #Additional weight for VBM and CBM if range we considered include VBM and CBM 
        if ( len(self.listBandWeight) >= self.nElectron /2 + 1):
            print("Unoccupiede states included in calculation, use different weight for HOMO and LUMO")
            if ( self.dBandHighWeight == None):
                self.dBandHighWeight = 10
            self.listBandWeight[self.nElectron/2-1] = self.dBandHighWeight
            self.listBandWeight[self.nElectron/2] = self.dBandHighWeight

        if ( self.nBandTotal == None): #Not specified
            #Additional very small weight to make result not fluctuate too much
            self.listBandWeight.append(0.001)
        else:
            #Add all additional band to make result stabilize
            self.listBandWeight = self.listBandWeight + [0.001 for x in range(0,self.nBandTotal-len(self.listBandWeight))]


        self.dVBMRef = f_Band_GetVBMFromElectron(self.bandRef,self.nElectron)
        self.dGapRef,dt,nt = f_Band_GetGap(self.bandRef,self.dVBMRef,False)

    def LoadBMBandWeight(self):
        '''
        '''
class EHParaConstrain:
    '''
    Constrain of parameters
    Possible contrains: Fix all parameters, only single-zeta for d orbital
    '''
    def __init__(self):
        #  @var stElement
        #  The element name
        self.stElement= "" #Whether to remove 
        self.bFixAll = False #
        self.bFixDZ = False  #


class PEMiniMultiInput:
    '''
    Input of multi-compound parameter fit
    '''
    def __init__(self):
        self.EHParaFileInput=""
        self.EHParaFileOutput=""
        self.listCompoundConfig=[]
        self.listEHParaConstrain=[]
#Non-serilize part
        self.yaehDatabase = None #A fake yaeh input file. Used to I/O all orbital parameters.

    def Load(self):
        '''
        Load all single-case data and cache them after sin
        '''
        stCwd = os.getcwd()
#Load EH Parameters
        self.yaehDatabase = yaeh_input(self.EHParaFileInput)
#Set data constrain
        if ( self.listEHParaConstrain != None):
            for constrain in self.listEHParaConstrain:
                aPara = self.yaehDatabase[constrain.stElement]
                if( constrain.bFixAll):
                    aPara.FixAll()
                if ( constrain.bFixDZ):
                    aPara.FixDZ()
        
#Deal with single calculation in each folder

        for sIn in self.listCompoundConfig:
#Load indivdual file
            stFolder = sIn.stCompound
            print("Set up %s ..."% stFolder)
            os.chdir(stFolder)

            sIn.Load()
            sIn.ParaBase = self.yaehDatabase.listPara #Link every single-compound calculation to this multi-c
#Clear output
            f = open("single.out","w")
            f.close()

            os.chdir(stCwd)

