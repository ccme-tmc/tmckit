#!/usr/bin/env python
'''@package w2k_caseutil 
Contains classes and functions for WIEN2k
'''

import copy
import shutil
from optparse import OptionParser
import subprocess
import re

from common_caseutil import *
from latt_utils import f_Latt_Get_ibrav
import chem_utils
from constants import Ry2eV
from list_utils import f_List_permute
import list_utils as lu
from io_utils import io_grep_lines
from band_utils import BandsT 
from dos_utils import OrbitalNameT

def ListCompare(x1,x2):
    return cmp(x1[1],x2[1])


def f_w2k_GetW2kDir():
    '''
    return Wien2K $WIENROOT Enviroment Variable
    '''
    return os.environ["WIENROOT"]
    
def f_w2k_GetW2kFilePath(stPath):
    '''
    return absolute path of a file under $WIENROOT directory
    '''
    return os.path.join(f_w2k_GetW2kDir(),stPath)


class Wien2K_Structure:
    def __init__(self,stFileName=""):
        self.stCaseName = ""  #From file name
        self.stTitle = "" #Title 
        self.stSpaceGroupSign = ""  # the internal character indicate the cell vectors, maybe P/F/I/...
        self.stSpaceGroupName = "" #the space group name
        self.stCalcMode = "RELA" #REAL or NREL, indicate real/complex calculation
        self.stUnit = "bohr" #unit, no effect
        self.LatticeParameter=[0,0,0]
        self.LatticeAngle=[0,0,0]
        self.listAtom = [] #list of non-equiv atoms
        self.listSymmetry = []
        if ( stFileName != ""):
            self.ReadFromFile(stFileName)
        
    def WriteToFile(self,stFileName):
        '''
        Write Wien2K struct to specific file
        '''
        fOut = open(stFileName,mode="w")
        
        fOut.write(self.stTitle)
        fOut.write("\n")
        
        fOut.write("%-4s" % self.stSpaceGroupSign)
        fOut.write("LATTICE,NONEQUIV.ATOMS:")
        fOut.write("%3d" % len(self.listAtom))
        fOut.write(self.stSpaceGroupName)
        fOut.write("\n")

        fOut.write("MODE OF CALC=")
        fOut.write(self.stCalcMode)
        fOut.write(" unit=")
        fOut.write(self.stUnit)
        fOut.write("                                                   \n")
        
        for i in range(0,3):
            fOut.write("%10.6f" % self.LatticeParameter[i])
        for i in range(0,3):
            fOut.write("%10.6f" % self.LatticeAngle[i])
        fOut.write("                   \n")
        
        for aAtom in self.listAtom:
            fOut.write(aAtom.ToString())
        
        fOut.write("%4d      NUMBER OF SYMMETRY OPERATIONS\n" % len(self.listSymmetry))    
        for i,aSym in enumerate(self.listSymmetry):
            for j in range(0,3):
                fOut.write("%2d%2d%2d%11.8f\n"% (aSym[0][j][0],aSym[0][j][1],aSym[0][j][2],aSym[1][j]))
            fOut.write("%8d\n" % (i+1))
                     
        fOut.close()
        print("Structure write into %s" % stFileName)
    
    def ReadFromFile(self,stFileName):
        '''
        Read Wien2K struct file
        '''
        self.stCaseName = os.path.splitext(os.path.basename(stFileName))[0]
        
        fIn = open(stFileName)
        self.stTitle = fIn.readline()[:-1]
        
        stLine = fIn.readline()[:-1]
        self.stSpaceGroupSign = stLine[0:4].strip()
        nAtomCount = int(stLine[27:30])
        self.stSpaceGroupName = stLine[30:]
        
        stLine=fIn.readline()
        self.stCalcMode = stLine[13:17]
        self.stUnit = stLine[23:28]
        
        stLine=fIn.readline()[:-1]
        for i in range(0,3):
            self.LatticeParameter[i] = float(stLine[i*10:(i+1)*10])
        for i in range(3,6):
            self.LatticeAngle[i-3] = float(stLine[i*10:(i+1)*10])
        
        #read atom
        for i in range(0,nAtomCount):
            aAtom = Wien2K_Atom()
            aAtom.ReadFromFile(fIn)
            self.listAtom.append(aAtom)
        
        #read symmetry (no parse!)
        stLine = fIn.readline()
        nSymCount = int(stLine[:4])
        self.listSymmetry = []
        for i in range(0,nSymCount):
            mSym =[] 
            mSymMove = []
            for j in range(0,3):
                stLine = fIn.readline()
                mSym.append([int(stLine[0:2]),int(stLine[2:4]),int(stLine[4:6])])
                mSymMove.append(float(stLine[6:]))
            self.listSymmetry.append([mSym,mSymMove])
            fIn.readline() # skip index
            
        fIn.close()

    def HasInversion(self):
        '''
        Return whether the structure has inversion symmetry
        This is used to judge real/complex version
        '''
        for mSym,mSymMove in self.listSymmetry:
            if ( mSym[0][0] == -1 and mSym[1][1] == -1 and mSym[2][2] == -1):
                return True
        return False
    
    def GetCellRaw(self):
        '''
        Get Raw cell information
        This function is especially useful as wien2k cell vector convention is different from Lattice class
        '''
        aCell = Lattice()
        mBravaisPrim = self.__GetPrimitiveVector__()
        mBravaisConv = self.__GetConvetionalCellVector__()
        #Convert R vectors only here
        aCell.ReadFromRawCellVector(mBravaisConv,mBravaisPrim)
        if ( self.stSpaceGroupSign == "R" ):
            stUseVec = "prim"
        else:
            stUseVec = "conv"
        for aAtom in self.listAtom:                
            for aCoord in aAtom.InnerCoord:
                aCell.AddAtom([[ chem_utils.f_Element_Z_to_Symbol(int(aAtom.fZ))] + aCoord],unit="crystal",latt=stUseVec)

        return aCell

    def GetCell(self):
        '''
        Get standard cell information
        May not give incorrect result, if so, use GetCellRaw instead
        '''
        #get cell type
        aCell = Lattice()
        
        arTest = [ ['X',[0,2,1]]  ]#dummy variable
        ibrav = f_Latt_Get_ibrav(self.stSpaceGroupSign,self.LatticeParameter+self.LatticeAngle,arTest)[0]
        
        nRot = arTest[0][1][0] # premute times
        
        #Get primitive cell vector
        mPrim = self.__GetPrimitiveVector__()
        #print(mPrim)
        #Permute a/b/c and x/y/z to make it the same as standard
        for i in range(0,3): #rotate
            mPrim[i] = f_List_permute(mPrim[i],nRot)
        mPrim = f_List_permute(mPrim,nRot)
        #print(mPrim)
        aCell.ReadFromPrimitiveCellVector(mPrim, ibrav)
        #print(aCell.fLatticeLength,aCell.fLatticeAngle)
        #print(aCell.PrimitiveCellVector)
        
        #Add Atom
        if ( self.stSpaceGroupSign == "R" ): #Primitive crystal coordinate
#R do not need to permute axis so we do not do it here
            for aAtom in self.listAtom:                
                for aCoord in aAtom.InnerCoord:
                    aCell.AddAtom([[ chem_utils.f_Element_Z_to_Symbol(int(aAtom.fZ))] + aCoord],"crystal")
        else:
#            mBravis = ConvertLatticeToCartesian(self.LatticeParameter,self.LatticeAngle)
            mBravais = self.__GetConvetionalCellVector__()
            mBravais = f_List_permute(mBravais,nRot)
            for i in range(0,3): #rotate
                mBravais[i] = f_List_permute(mBravais[i],nRot)

            for aAtom in self.listAtom: #Conventional Cell crystal coordinates
                for aCoord in aAtom.InnerCoord:
                    #convert to Cartesian
                    aCell.AddAtom([ [chem_utils.f_Element_Z_to_Symbol(int(aAtom.fZ))] + lu.f_Matrix_transpose(lu.f_Matrix_dot(lu.f_Matrix_transpose(mBravais), lu.f_Matrix_transpose([f_List_permute(aCoord,nRot)])))[0]],"bohr")
                    #Direct Add
                    #aCell.AddAtom( [[chem_utils.elements[int(aAtom.fZ)-1]] + f_List_permute(aCoord,nRot)],"crystal","normal")
        aCell.CorrectError()
        return aCell                    
    
    def __GetPrimitiveVector__(self):
        '''
        Get primitive cell vector followed by WIEN2K convention
        strip from WIEN2K source ( the described Bravais matrix in document v091 about P is totally wrong, reasonable in v11 document, but it seems that neither is the same as lapw1 source code!)
        '''
        (a,b,c)=self.LatticeParameter
        (aa,ab,ac)= [x /180.0 * math.pi for x in self.LatticeAngle]
        
        dicPrim = {}
        
#        dicPrim["P"]=[[a*math.sin(ac)*math.sin(ab),a*math.cos(ac)*math.sin(ab),a*math.cos(ab)],[0,b*math.sin(aa),b*math.cos(aa)],[0,0,c]]       
        dicPrim["P"]=[[a/math.sin(aa)*math.sqrt(1-math.cos(aa)**2-math.cos(ab)**2-math.cos(ac)**2+2*math.cos(aa)*math.cos(ab)*math.cos(ac)),a*(math.cos(ac)-math.cos(ab)*math.cos(aa))/math.sin(aa),a*math.cos(ab)],[0,b*math.sin(aa),b*math.cos(aa)],[0,0,c]]
        dicPrim["R"]=[[a/math.sqrt(3)/2,-a/2.0,c/3],[a/math.sqrt(3)/2.0,a/2.0,c/3.0],[-a/(math.sqrt(3)),0,c/3.0]]
        dicPrim["F"]=[[a/2, b/2, 0], [a/2, 0, c/2], [0, b/2, c/2]]
        #dicPrim["B"] = [[a/2, -b/2, c/2],[a/2, b/2, -c/2], [-a/2, b/2, c/2]]
#Warning : Wien2K (v091)  user manual says B should be as above, but in program nn.f it is actually the one below !!
#May meet problem in calculation must have right order of abc like band structure ( where k-path rely on abc axis order )
        dicPrim["B"] = [[-a/2, b/2, c/2],[a/2, -b/2, c/2], [a/2, b/2, -c/2]]

        dicPrim["CXY"] = [[a/2,-b/2,0],[a/2,b/2,0],[0,0,c]]
        dicPrim["CYZ"] = [[a,0,0],[0,-b/2,c/2],[0,b/2,c/2]]
        dicPrim["CXZ"] = [[a*math.sin(ac)/2,a*math.cos(ac)/2,-c/2],[0,b,0],[a*math.sin(ac)/2,a*math.cos(ac)/2,c/2]]
        dicPrim["H"] = [[math.sqrt(3)/2*a,-a/2,0],[0,a,0],[0,0,c]]
        
        
        for stSymbol,arPrim in dicPrim.iteritems():
            if (self.stSpaceGroupSign.find(stSymbol) != -1):
                mPrim = arPrim
                break
        
        if ( mPrim == None):
            print("Unrecgonized symmetry symbol! Please check your .struct files")
            return None
        
        return mPrim                

    def __GetConvetionalCellVector__(self):
        '''
        Get primitive cell vector followed by WIEN2K convention
        CXZ calculated by hand
        '''
        (a,b,c)=self.LatticeParameter
        (aa,ab,ac)= [x /180.0 * math.pi for x in self.LatticeAngle]
        
        mOrtho = [[a,0,0],[0,b,0],[0,0,c]]
        dic = {}
        dic["F"]= mOrtho   
        dic["B"] = mOrtho
        dic["CXY"] = mOrtho
        dic["CYZ"] = mOrtho
        dic["CXZ"] = [[a*math.sin(ac),a*math.cos(ac),0],[0,b,0],[0,0,c]]

        mConv = None
        for stSymbol,arCell in dic.iteritems():
            if (self.stSpaceGroupSign.find(stSymbol) != -1):
                mConv = arCell
                break
        if ( mConv == None):
#Not found, it means it is the same as primitive cell
            mConv = self.__GetPrimitiveVector__()

        if ( mConv == None):
            print("Unrecgonized symmetry symbol! Please check your .struct files")
            return None

        return mConv

        
    def GetPrimitiveCell(self):
        '''
        Get primitive cell from Wien2K struct
        '''
        result=copy.deepcopy(self)
        result.stSpaceGroupSign="P" # mark primitive
        result.stSpaceGroupName=" " # remove group name
        #mBravis=numpy.zeros((3,3))
        
        mPrim = self.__GetPrimitiveVector__()
        
        print("Convert coordinates")
        print("Wien2k Primitive Bravis Matrix")
        print(mPrim)
        
        mBravis = ConvertLatticeToCartesian(self.LatticeParameter,self.LatticeAngle)
        print("New Bravis Matrix")
        print(mBravis)
        result.LatticeParameter,result.LatticeAngle = ConvertCartesianToLattice(mPrim)
        
        if ( not ("R" in self.stSpaceGroupSign or "H" in self.stSpaceGroupSign or "P" in self.stSpaceGroupSign)): 
            #don't change if it is R ( internal coordinate already primitive cell though the abc is not) or P,H (nothing changed)
            print("Convert...")
            mPrimTI = lu.f_List_inv3(lu.f_Matrix_transpose(mPrim))
            mBravisT = lu.f_Matrix_transpose(mBravis)
            for i,aAtom in enumerate(self.listAtom):
                for j,aCoord in enumerate(aAtom.InnerCoord):
                    vCoord = lu.f_Matrix_transpose([aCoord])
                    vNew = lu.f_Matrix_dot(mBravisT,vCoord)
                    vNew = lu.f_Matrix_transpose(lu.f_Matrix_dot(mPrimTI,vNew))[0]
                    print(vNew)
                    print(result.listAtom[i].InnerCoord[j])
                    result.listAtom[i].InnerCoord[j] = [ (x+1) if (x < 0) else ( (x -1) if (x >= (1.0-0.001)) else x ) for x in vNew]
        
        result.listSymmetry = []
        
        return result
    
    def CheckDuplicate(self,bOutput = True):
        '''
        Check whether there are duplicate atoms ( note only same kind of atom is checked, overlap of different atom is not considered )
        Move atoms in the range of precision error of machine of some specific value to precision
        '''
        fLimit = 0.0001
        for aAtom in self.listAtom:
            i = 0
            while ( i < len(aAtom.InnerCoord)):
                listDelete = []
                for j in range(0,aAtom.InnerCoord):
                    if ( abs(aAtom.InnerCoord[i][1] - aAtom.InnerCoord[j][1]) < fLimit and \
                         abs(aAtom.InnerCoord[i][2] - aAtom.InnerCoord[j][2]) < fLimit and \
                         abs(aAtom.InnerCoord[i][3] - aAtom.InnerCoord[j][3]) < fLimit
                          ):
                        listDelete.append(j)
                listDelete.reverse()
                for j in listDelete:
                    if (bOutput):
                        print("Atom %s at %f %f %f is duplicated" % (aAtom.stName,\
                                                                     aAtom.InnerCoord[j][1],\
                                                                     aAtom.InnerCoord[j][2],\
                                                                     aAtom.InnerCoord[j][3]))
                    del aAtom[j]
        pass

        
class Wien2K_File_in0:
    def __init__(self,stFileName=""):
        self.stSwitch = "TOT"
        self.nPotentialIndex=5
        self.stPotential="XC_LDA"
        self.stPrintMode="NR2V"
        self.stHybridMode="    "
        self.stFFTMode="    "
        self.fLUSE=0
        self.arFFT=[-1,-1,-1]
        self.fFFTFactor=0
        self.nIField=0
        self.fEField=0
        self.nWField=0
        self.nVersion=13
        if ( stFileName != ""):
            self.ReadFromFile(stFileName)
        
    def WriteToFile(self,stFileName):
        '''
        Write in0 file
        '''
        fOut = open(stFileName,mode="w")
        
        fOut.write("%-6s" % self.stSwitch)
#This is different for version 14+
        if (self.nVersion < 14):
            fOut.write("%2d" % self.nPotentialIndex)
            fOut.write("    (5...CA-LDA, 13...PBE-GGA, 11...WC-GGA)\n")
        else:
            fOut.write(self.stPotential)
        
        fOut.write("%4s %4s %4s " % (self.stPrintMode,self.stHybridMode,self.stFFTMode ))
        if ( self.fLUSE != 0):
            fOut.write("%4f" % self.fLUSE)
        else:
            fOut.write("    ")
        if ( self.stFFTMode != "    "):
            fOut.write(" (R2V)\n")        
            fOut.write("%4d%4d%4d%8.2f    min IFFT-parameters, enhancement factor\n" %(self.arFFT[0],self.arFFT[1],self.arFFT[2],self.fFFTFactor))
        fOut.close()
        print("in0 write into %s" % stFileName)
        
    def ReadFromFile(self,stFileName):
        '''
        Read in0 file
        '''
        fIn = open(stFileName)
        arLine=ReadSplitLine(fIn)
        self.stSwitch=arLine[0]
        if (arLine[1][0].isdigit()):
            self.nPotentialIndex=int(arLine[1])
        else:
#Version 14+
            self.stPotential = arLine[1]
            self.nVersion = 14

        arLine=ReadSplitLine(fIn)
        #1part: NR2V only 3part: no EECE/HYBD  4part: have EECE 5part have L 
        self.stPrintMode=arLine[0]
        if ( len(arLine) == 3 ):
            self.stFFTMode=arLine[1]
        elif ( len(arLine) == 4):
            self.stHybridMode = arLine[1]
            self.stFFTMode=arLine[2]
        elif ( len(arLine) == 5):
            self.stHybridMode = arLine[1]
            self.stFFTMode=arLine[2]
            self.nLUSE=int(arLine[3])
            
        if ( self.stFFTMode == "IFFT"    ):
            arLine=ReadSplitLine(fIn)
            self.arFFT = [int(x) for x in arLine[0:3]]
            self.fFFTFactor = float(arLine[3])
        
        arLine = ReadSplitLine(fIn)
        if ( len(arLine) > 0 ):
            self.nIField = int(arLine[0])
            self.fEField = float(arLine[1])
            self.nWField = int(arLine[2])
            
        
class Wien2K_File_in1:
    def __init__(self,stFileName=""):
        self.stSwitch="WFFIL"
        self.dFermi = 0.5
        self.fRKMax = 7.00
        self.nLMax = 10
        self.nLnsMax = 4
        self.arAtom = []
        self.nKFileNumber = 4
        self.fEMin = -9.0
        self.fEMax = 2.0
        self.nBand = 32
        if ( stFileName != ""):
            self.ReadFromFile(stFileName)
        
    def ReadFromFile(self,stFileName):
        fIn = open(stFileName)
        
        #arLine = ReadSplitLine(fIn)
        #self.stSwitch = arLine[0]
        #if ( "EF"in arLine[1]):
        #    self.dFermi =  
        stLine = fIn.readline()
        arLine = [x.strip() for x in stLine.split()]
        self.stSwitch = arLine[0]
        i = stLine.find("EF")
        if ( i != -1):
            iend = stLine[i+3:].find("(") + i+3
            self.dFermi = float(stLine[i+3:iend])
        
        arLine = ReadSplitLine(fIn)
        self.fRKMax = float(arLine[0])
        self.nLMax = int(arLine[1])
        self.nLnsMax = int(arLine[2])
        
        while True:
            arLine = ReadSplitLine(fIn)
            if ( arLine[0] == "K-VECTORS"):
                break
            aAtom = Wien2K_File_in1_atom()
            aAtom.fE = float(arLine[0])
            aAtom.nLAPW = int(arLine[2])
            for i in range(0,int(arLine[1])):
                aOrbital = Wien2K_File_in1_orbital()
                arLine = ReadSplitLine(fIn)
                aOrbital.nL = int(arLine[0])
                aOrbital.fE = float(arLine[1])
                aOrbital.fDE = float(arLine[2])
                aOrbital.stSwitch=arLine[3]
                aOrbital.nLAPW = int(arLine[4])
                
                aAtom.arOrbital.append(aOrbital)
            
            self.arAtom.append(aAtom)
        
        self.nKFileNumber = int(arLine[2][-1])
        self.fEMin = float(arLine[3])
        self.fEMax = float(arLine[4])
        self.nBand = int(arLine[5])
        
    def WriteToFile(self,stFileName):
        fOut = open(stFileName,mode="w")
        
        fOut.write("%5s  EF= %7.4f    (WFPRI, SUPWF)\n"%(self.stSwitch,self.dFermi))
        
        fOut.write("%6.2f       %2d    %1d (R-MT*K-MAX; MAX L IN WF, V-NMT\n"%(self.fRKMax,self.nLMax,self.nLnsMax) )
        
        for aAtom in self.arAtom:
            fOut.write("%6.2f   %2d  %1d      (GLOBAL E-PARAMETER WITH n OTHER CHOICES, global APW/LAPW)\n" % (aAtom.fE,len(aAtom.arOrbital),aAtom.nLAPW ))
            for aOrbital in aAtom.arOrbital:
                fOut.write(" %1d%8.2f      %5.3f %4s %1d\n" % ( aOrbital.nL, aOrbital.fE,aOrbital.fDE,aOrbital.stSwitch,aOrbital.nLAPW ) )
        
        fOut.write("K-VECTORS FROM UNIT:%1d%7.1f%10.1f%6d   red   emin/emax/nband\n"%(self.nKFileNumber,self.fEMin,self.fEMax,self.nBand))
        print(".in1 write into %s" % stFileName)
                
        
        

class Wien2K_File_in1_atom:
    '''
    a non-eq atom in in1 file
    '''
    def __init__(self):
        self.fE = 0.0
        self.nLAPW = 0
        self.arOrbital=[]

class Wien2K_File_in1_orbital:
    '''
    a orbital in a atom in in1 file
    '''
    def __init__(self):
        self.nL = 0
        self.fE = 0.0
        self.fDE = 0.0
        self.stSwitch ="CONT"
        self.nLAPW = 0
        

class Wien2K_Atom:
    '''
    Represent a non-equivalent atom in case.struct
    '''
    def __init__(self):
        self.listIndex = [] # atom-index
        self.InnerCoord = [] #x ,y ,z
        self.nSplit = 1
        self.stName="H"
        self.nCharge=1
        self.nNPT = 781
        self.fR0 = 0.0001
        self.fRMT = 2.00
        self.fZ = 1.00
        self.mRot = [[0 for col in range(3)] for row in range(3)]
    
    def ReadFromFile(self,fIn):
        '''
        Read atom info from file stream
        '''
        stLine = fIn.readline()
        self.listIndex.append(int(stLine[4:8]))
        aAtomCoord = [float(stLine[12:22]),float(stLine[25:35]),float(stLine[38:48])]
        self.InnerCoord.append(aAtomCoord)
        
        stLine = fIn.readline()
        nAtomCount = int(stLine[15:17])
        self.nSplit = int(stLine[34:36])
        
        for i in range(0,nAtomCount-1):
            stLine = fIn.readline()
            self.listIndex.append(int(stLine[4:8]))
            aAtomCoord = [float(stLine[12:22]),float(stLine[25:35]),float(stLine[38:48])]
            self.InnerCoord.append(aAtomCoord)
            
        stLine=fIn.readline()
        self.stName = stLine[0:11]
        self.nNPT = int(stLine[15:20])
        self.fR0 = float(stLine[25:35])
        self.fRMT = float(stLine[40:50])
        self.fZ = float(stLine[55:65])
        
        
        for i in range(0,3):
            stLine = fIn.readline()
            self.mRot[i] = [float(stLine[20:30]),float(stLine[30:40]),float(stLine[40:50])]
            

        
        
    def ToString(self):
        result = ""
        aAtom = self.InnerCoord[0]
        #add first atom
        result = result+"ATOM" + "%4d" % self.listIndex[0] + ": X=" + "%10.8f" % aAtom[0] + " Y=" + "%10.8f" % aAtom[1] + " Z=" + "%10.8f" % aAtom[2] + "\n"
        #add multi and split
        result = result + "          MULT=" + "%2d" % len(self.InnerCoord) + "          ISPLIT=" + "%2d" % self.nSplit + "\n"
        #add other atom
        for i in range(1,len(self.InnerCoord)):
            aAtom = self.InnerCoord[i]
            result = result+"    " + "%4d" % self.listIndex[i] + ": X=" + "%10.8f" % aAtom[0] + " Y=" + "%10.8f" % aAtom[1] + " Z=" + "%10.8f" % aAtom[2] + "\n"
        #add summary line
        result = result + self.stName +  "NPT=" + "%5d" % self.nNPT + "  R0=" + "%10.8f" % self.fR0 + " RMT=" + "%10.5f" % self.fRMT + "   Z:" + "%5.2f" % self.fZ + "              \n"
        #add rot matrix
        result = result + "LOCAL ROT MATRIX:   " + "%10.7f" % self.mRot[0][0] + "%10.7f" % self.mRot[0][1] + "%10.7f" % self.mRot[0][2] + "\n"
        result = result + "                    " + "%10.7f" % self.mRot[1][0] + "%10.7f" % self.mRot[1][1] + "%10.7f" % self.mRot[1][2] + "\n"
        result = result + "                    " + "%10.7f" % self.mRot[2][0] + "%10.7f" % self.mRot[2][1] + "%10.7f" % self.mRot[2][2] + "\n"
        
        return result
        

class Wien2K_File_in2:
    def __init__(self,stFileName=""):
        self.stMode="TOT"
        self.stEECE=""
        self.fEMin=0
        self.fElectron=0
        self.fEValMin=0
        self.fEValGap=0
        self.stEFermiMode="TETRA"
        self.fBroadening=0
        self.list_stLM=[]
        self.fGMax=0
        self.stKListFile="FILE"
        if ( stFileName != ""):
            self.ReadFromFile(stFileName)
        pass
    
    def ReadFromFile(self,stFileName):
        fIn = open(stFileName,'r')
        list_stLine = fIn.readlines()
        arLine = list_stLine[0].split()
        self.stMode = arLine[0]
        if ( len(arLine) >2  ):
            self.stEECE = arLine[1]
        self.fEMin,self.fElectron,self.fEValMin,self.fEValGap  = [ float(x) for x in list_stLine[1].split()[:4]]
        arLine = list_stLine[2].split()
        self.stEFermiMode,self.fBroadening = arLine[0],float(arLine[1])
        i =3
        while ( list_stLine[i].find('GMAX') == -1 ):
            self.list_stLM.append(list_stLine[i])
            i = i +1
        self.fGMax = float(list_stLine[i].split()[0])
        self.stKListFile = list_stLine[i+1].split()[0]
        
    def WriteToFile(self,stFileName):
        fOut = open(stFileName,'w')
        fOut.write("%-5s%-5s      (TOT,FOR,QTL,EFG,FERMI)\n" % (self.stMode,self.stEECE))
        fOut.write("%10.1f%10.1f%5.2f%5.2f                EMIN, NE, ESEPERMIN, ESEPER0\n" % ( self.fEMin,self.fElectron,self.fEValMin,self.fEValGap ) )
        fOut.write("%-s5f10.5\n" % (self.stEFermiMode,self.fBroadening))
        for stLine in self.list_stLM:
            fOut.write(stLine)
            fOut.write("\n")
        fOut.write("%6.2f          GMAX\n" % self.fGMax)
        fOut.write("%-6s        FILE/NOFILE  write recprlist\n" % self.stKListFile)
        fOut.close()
        print(".in2 write into %s" % stFileName)
        
        
        

class Wien2K_Case:
    '''
    Represent a Case of Wien2K calculation.
    Defined by current directory and case name ( Casename can be different from folder name, unlike WIEN2K calculation )
    '''
#Spin-polarize related files ( which means one file will be two file with +"up" and +"dn", the new name in the keys)
    dic_file_sp = {
            "clmcor":"clmcor",
            "clmsc":"clmsc",
            "clmsum":"clm",
            "clmval":"clmval",
            "dmat":"dmat",
            "vector":"vector",
            "energy":"energy",
            "vorb":"vorb",
            "vsp":"vsp"
            }

    def __init__(self,stDirectoryInput="",stCaseNameInput="",bInit=True):
        if ( stCaseNameInput == ""):
            self.stCaseName = f_GetCaseName()
        else:
            self.stCaseName=stCaseNameInput
        self.aStructure = Wien2K_Structure()
        self.a_in0 = Wien2K_File_in0()
        self.a_in1 = Wien2K_File_in1()
        self.a_in2 = Wien2K_File_in2()
        self.a_dayfile = Wien2K_File_dayfile()
        self.a_scf= Wien2K_File_scf()
        if ( bInit):
            self.ReadFromFile(stDirectory=stDirectoryInput, stSuffix="", listPart=["struct","in0","in1","in2"])
            #self.RestoreCase(stDirectory=stDirectoryInput, bDefaultSave=False, listPart=["struct","in0","in1","in2"], bCopy=False)
            #self.aStructure = Wien2K_Structure(self.stCaseName+".struct")
            #self.a_in0 = Wien2K_File_in0(self.stCaseName+".in0")
            #self.a_in1 = Wien2K_File_in1(self.stCaseName+".in1")
            #self.a_in2 = Wien2K_File_in2(self.stCaseName+".in2")
        
        
    
    def ShowSummary(self):
        '''
        Output Wien2K case summary to screen.
        '''
        print("CaseName: "+self.stCaseName)
        print("RMT:")
        for aAtom in self.aStructure.listAtom:
            print(" Atom :%3s %4.2f" %(aAtom.stName, aAtom.fRMT))
        
        print("RKMAX: %4.2f" % self.a_in1.fRKMax)
        print("GMAX : %4.2f" % self.a_in2.fGMax)
        print("FFT  : %3d %3d %3d" % (self.a_in0.arFFT[0],self.a_in0.arFFT[1],self.a_in0.arFFT[2]))
        
        
        dicFunctional = {}
        dicFunctional[5]="LDA"
        dicFunctional[11]="WC06"
        dicFunctional[13]="PBE"
        dicFunctional[19]="PBEsol"
        dicFunctional[20]="AM05"
        
        print("Functional: %s" % dicFunctional[self.a_in0.nPotentialIndex])
        
        if ( os.path.exists(self.stCaseName+".dayfile")):
            self.ReadFromFile(listPart=["dayfile"])
            if ( self.a_dayfile.bSuccess):
                print("Status:OK")
            else:
                print("Status:Error")
                print("Crash part:%s" % self.a_dayfile.stCrash)
                print(self.a_dayfile.stError)
        
        if ( os.path.exists(self.stCaseName+".scf")):
            self.ReadFromFile(listPart=["scf"])
            
            
    def ReadFromFile(self,stDirectory="",stSuffix="",listPart=[]):
        '''
        Read object from Wien2K files in Folder 
        '''
        stPartName = os.path.join(stDirectory,self.stCaseName+stSuffix+".")
        if ( "struct" in listPart):
            self.aStructure.ReadFromFile(stPartName+"struct")
        if ( "in0" in listPart):
            self.a_in0.ReadFromFile(stPartName+"in0")
#Real or complex?
        if ( self.aStructure.HasInversion() ):
            stSuffix = ""
        else:
            stSuffix ="c"

        if ( "in1" in listPart):
            self.a_in1.ReadFromFile(stPartName+"in1"+stSuffix)
        if ( "in2" in listPart):
            self.a_in2.ReadFromFile(stPartName+"in2"+stSuffix)
        if ( "dayfile" in listPart):
            self.a_dayfile.ReadFromFile(stPartName+"dayfile")
        if ( "scf" in listPart):
            self.a_scf.ReadFromFile(stPartName+"scf", stPartName+"dayfile")
        return       
    
    
    def WriteToFile(self,stDirectory="",stSuffix="",listPart=[]):
        '''
        Write object to Wien2K file in a folder
        '''
        #Detect if directory exists
        if ( stDirectory != ""):
            if(not os.path.exists(stDirectory)):
                os.mkdir(stDirectory)
                
        stPartName = os.path.join(stDirectory,self.stCaseName+stSuffix+".")
        if ( "struct" in listPart):
            self.aStructure.WriteToFile(stPartName+"struct")
        if ( "in0" in listPart):
            self.a_in0.WriteToFile(stPartName+"in0")
#Real or complex?
        if ( self.aStructure.HasInversion() ):
            stSuffix = ""
        else:
            stSuffix ="c"
        if ( "in1" in listPart):
            self.a_in1.WriteToFile(stPartName+"in1"+stSuffix)
        if ( "in2" in listPart):
            self.a_in2.WriteToFile(stPartName+"in2"+stSuffix)
        return        
        
   
    def SaveCase(self,stDirectory="",stSuffix="",bDefaultSave=True,listPart=[]):
        '''
        Backup important files of case
        
        :param stDirectory: The directory to save files in
        :param stSuffix: The extra suffix to append in casename, like Graphite.scf to Graphite_xxx.scf
        :param bDefaultSave: Is save_lapw -a is used
        :param listPart: Save specific file ( in case.xxx format) in the list of xxx, it is must if bCopy=False
   
        '''
        #Detect if directory exists
        if ( stDirectory != ""):
            if(not os.path.exists(stDirectory)):
                os.mkdir(stDirectory)
        
        #Rebuild filename for spin-polarized
#@todo not implentend
       # listPart2 = []
       # for part in listPart:
       #     if (Wien2K_Case.dic_file_sp.has_key(x)):
       #         part2 = Wien2K_Case.dic_file_sp[x]
       #         listPart2.extend([
       #     else:
       #         listPart2.append(part)


        for stName in listPart:
            stOrgName = self.stCaseName+"."+stName
            stFullName = os.path.join(stDirectory,self.stCaseName+stSuffix+"."+stName)
            if ( os.path.exists(stFullName)):
                print("Error: File %s exists. %s Skipped." % (stFullName,stOrgName))
            else:
                print("File %s copied to %s" % (stOrgName,stFullName))
                shutil.copy(self.stCaseName+"."+stName,stFullName)
                    
        if ( bDefaultSave):
            commands.getoutput("save_lapw -a -d %s %s" %(stDirectory,self.stCaseName+stSuffix))
        
    def RestoreCase(self,stDirectory="",stSuffix="",bDefaultSave=True,listPart=[]):
        '''
        Restore important files of case

        # @param stDirectory The directory to save files in
        # @param stSuffix The extra suffix to append in casename, like Graphite.scf to Graphite_xxx.scf
        # @param bDefault Is save_lapw -a is used
        # @param listPart Don't use built-in save_lapw, but save specific file ( in case.xxx format) in the list of xxx
            
        ''' 
        for stName in listPart:
            stOrgName = self.stCaseName+"."+stName
            stFullName = os.path.join(stDirectory,self.stCaseName+stSuffix+"."+stName)
            if ( not os.path.exists(stFullName)):
                print("Error: File %s not exists. %s Skipped." % (stFullName,stOrgName))
            else:
                print("File %s copied to %s" % (stFullName,stOrgName))
                shutil.copy(stFullName,stOrgName)
        if(bDefaultSave):
            commands.getoutput("restore_lapw -f -d %s %s" %(stDirectory,self.stCaseName+stSuffix))

class Wien2K_File_dayfile:
    '''
    Wien2K dayfile parser
    '''
    def __init__(self,stFileName=""):
        self.stError = ""
        self.stCrash = ""
        self.nCycle = 0
        self.bSuccess = True
        self.nTime = 0
        if ( stFileName != ""):
            self.ReadFromFile(stFileName)
    
    def ReadFromFile(self,stFileName=""):
        if ( stFileName == "" or stFileName==None):
            stCaseName = f_GetCaseName()
            stFileName = stCaseName+".dayfile"
        fIn = open(stFileName,'r')
        list_stAll=fIn.readlines()
        if ( list_stAll[-1].find("stop error") == -1):
            self.bSuccess = True
        else:
            self.bSuccess = False
        nStart = 0
        nEnd = 0
        for stLine in list_stAll:
            if ( stLine.find("crashed") != -1):
                self.stCrash=stLine.split()[1]
            elif( stLine.find("error:") != -1):
                self.stError = stLine[:-1]
            elif( stLine.find("cycle") != -1):
                self.nCycle = int(stLine.split()[1])
    
    
class Wien2K_machines:
    '''
    Wien2K parallel execution file
    '''
    def __init__(self,stType="PBS",bInit=True):
        self.listHost = []
        self.listHostWeight =[]
        self.stType = stType
        if ( bInit ):
            self.ReadFromHost(stType)
    
    def ReadFromHost(self,stType="PBS"):
        '''
        Generate host list from type string
        '''
        if ( stType == None or stType == ""):
            self.stType = f_env_GetJobSystem()
            stType = self.stType
        
        if ( stType == "PBS"):
            fIn = open(os.environ["PBS_NODEFILE"],"r")
            ar_stHost = fIn.readlines()
            for stHost in ar_stHost:
                self.listHost.append(stHost.strip())
                self.listHostWeight.append(1)
        elif ( stType =="LSF"):
            stHost = open(os.environ["LSB_MCPU_HOSTS"])
            ar_stHost = stHost.split(" ")
            for i in range(0,len(ar_stHost)/2):
                self.listHost.append(ar_stHost[i*2])
                self.listHostWeight.append(int(ar_stHost[i*2+1]))
            
        
    
    def WriteToFile(self,stFileName=".machines",bAllOne=True,bOverwrite=False,bMPI=False):
        '''
        Generate .machines file
        
        :param bAllOne: Treat every process eqaully, and generate weight 1 for all single process(hostname)
        :param bOverwrite: whether overwrite existed files
        :param bMPI: use MPI parallel in .machines file
        '''
        if ( (bOverwrite == False) and os.path.exists(stFileName) ):
            print("Skip existed .machine files")
            return
        fOut = open(stFileName,"w")
        fOut.write("# This is a machines file from w2k_caseutil.py\n")
        fOut.write("#\n")

        
        for i in range(0,len(self.listHost)):
            if ( bAllOne ):
                for nCount in range(0,self.listHostWeight[i]):
                    fOut.write("1:%s\n" % self.listHost[i])
            else:
                fOut.write("%d:%s\n" % (self.listHostWeight[i],self.listHost[i]) )
        fOut.write("granularity:1\n")
        fOut.write("extrafine:1\n")
        
        if ( bMPI):#write MPI part
            fOut.write("lapw1:")
            for i in range(0,len(self.listHost)):
                fOut.write("%s:%d " % (self.listHost[i],self.listHostWeight[i]) )
        
        fOut.write("\n")
            
class Wien2K_File_scf:
    def __init__(self,stCaseName=""):
        self.bConverge = True
        self.fEnergy = 0
        self.fGap = 0
        self.list_fForce = []
        self.a_dayfile=Wien2K_File_dayfile()
        if( stCaseName != ""):
            self.ReadFromFile(stCaseName+".scf",stCaseName+".dayfile")
    
    def ReadFromFile(self,stFileName="",stLogFileName=""):
        '''
        Read Gap,Energy and Force from SCF file, see dayfile to verify if it is complete
        '''
        
        #if no file exist, treat as not converge
        #if not converge, exit directly
        if ( not os.path.exists(stFileName) ):
            print(stFileName+" not exsits!")
            self.bConverge = False
            return
        
        if ( os.path.exists(stLogFileName) ):
            self.a_dayfile.ReadFromFile(stLogFileName)
            self.bConverge = self.a_dayfile.bSuccess
            if ( self.bConverge == False):
                return
        
        fIn = open(stFileName,'r')
        list_stAll = fIn.readlines()
        #reverse search specific string
        bFindForce = False
        bFindEnergy = False
        bFindGap = False
        bFindFermi = False
        for i in range(len(list_stAll)-1,-1,-1):
            stLine = list_stAll[i]
            if ( (not bFindForce) and (stLine.find("TOTAL FORCE") != -1) ):
                j = i+1
                stLine = list_stAll[j]
                while ( stLine != "\n" and stLine !=""):
                    ar_stLine = stLine.split()
                    self.list_fForce.append([float(x) for x in ar_stLine[2:5]])
                    j = j + 1
                    if ( j >= len(list_stAll)):
                        break
                    stLine = list_stAll[j]
                bFindForce = True
            if ( (not bFindEnergy) and ( stLine.find(":ENE") != -1)):
                self.fEnergy = float(stLine.split()[-1])
                bFindEnergy = True
            if ((not bFindGap) and ( stLine.find(":GAP") != -1)):
                self.fGap = float(stLine.split()[5])
                bFindGap = True
            if ((not bFindFermi) and ( stLine.find(":FER") != -1)):
                self.fFermi = float(stLine.split()[-1])
                bFindFermi = True
            if ( bFindGap and bFindEnergy and bFindForce and bFindFermi):
                break
        fIn.close()

            
   
def w2k_CreateIn0():
    '''
    use dstart -fft to generate a new in0 with proper fft parameter
    '''      
    status,output=commands.getstatusoutput("x dstart -fft")
    if ( status != 0):
        print("Error when create case.in0")
        return
    
    stCaseName = f_GetCaseName()
    shutil.copy(stCaseName+".in0_std", stCaseName+".in0")
    
    pass


def w2k_CreateKPointList(listkp,bAddInv=True,bShift=False):
    '''
    use x kgen program to create .klist file for scf non-iteractive
    :param listkp: list, contains one integer as the k-points total count, or three integers as k-points in each dimension
    :param bAddInv: whether to add inversion symmetry in k-points if asked, default yes. 
    :param bShift: whether to shift k-points, default no.
    '''

    po = subprocess.Popen(["x","kgen"],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    while (po.poll() == None):            
        out = po.stdout.readline().strip().lower()
        if  ("add inversion" in out):
            po.stdin.write("1\n")
            po.stdin.flush()
        elif ( "0 allows to specify" in out):
            if ( len(listkp) == 1):
                po.stdin.write("%i\n" % listkp[0] )
            elif ( len(listkp) == 3):
                po.stdin.write("0\n")
            else:
                raise ValueError,"To create k-points 1 or 3 integers must be specified"
            po.stdin.flush()
        elif ("3 mesh" in out):
            po.stdin.write("%i\n%i\n%i\n"% tuple(listkp))
            po.stdin.flush()
        elif ("shift of k-mesh allowed" in out):
            if ( bShift):
                po.stdin.write("1\n")
            else:
                po.stdin.write("0\n")
            po.stdin.flush()
        elif ("short" in out):
            po.stdin.write("0\n")
        else:
            pass

    return 0



def w2k_GetKPointList(stStructName,stFileName):
    '''
    Read WIEN2k case.klist file to KPoint class
    '''
    aStruct = Wien2K_Structure(stStructName)

    f = open(stFileName,"r")
    ar_stLine = f.readlines()
    #Determine it is 10 character split or 5 character split WIEN2k kpoint file
    nSplit = 10
    if ( ar_stLine[0][14] != " " ): # not empty at 15 character means it is 5
        nSplit = 5
    arValue = []
    i = 0 

    while ( ar_stLine[i][0:3] != "END"):
        stLine = ar_stLine[i]
        nBase = int(stLine[nSplit*3+10:nSplit*4+10])
        arValue.append( [ stLine[0:10].strip(),int(stLine[nSplit*0+10:nSplit*1+10])*1.0/nBase,int(stLine[nSplit*1+10:nSplit*2+10])*1.0/nBase,int(stLine[nSplit*2+10:nSplit*3+10])*1.0/nBase,float(stLine[nSplit*4+10:nSplit*4+15] ) ])
        #print(arValue[-1])
        i += 1

    aKPT = KPointsT()
    aKPT.listKPt = arValue
    aKPT.latt = aStruct.GetCellRaw()
#   print(aKPT.latt.ConventionalCellVector)
#   print(aKPT.latt.PrimitiveCellVector)
    arTest = [ ['X',[0,1,2]]  ]#dummy variable
    ibrav = f_Latt_Get_ibrav(aStruct.stSpaceGroupSign,aStruct.LatticeParameter+aStruct.LatticeAngle,arTest)[0]

    #print(aKPT.latt.nLatticeType)
    if ( ibrav <=3 or ( ibrav >= 6 and ibrav <= 11 ) ):
        aKPT.stMode = "tpibabc"
    else:
        aKPT.stMode = "crystal"

    return aKPT

def w2k_WriteKPoints(KPt,stFileName,nDivisor=60):
    '''
    Write k-points object to specified file ( values are not checked ! ) in 10-5-5-5-5 format. 
    The k-points coordinate must be stored as integer1/integer2
    The integer2 is fixed to 60 by default
    If band plot is necessary, then the user must guarantee every point in one line have exactly same seperation after contract to integer!
    :param KPt: the KPointsT object
    :param stFileName: the file name to write in
    :param nDivisor: the idv in .klist file
    '''
    print("Warning: Please ensure your k-points coordinate unit follows wien2k convention!")
    f = open(stFileName,"w")

    for i,line in enumerate(KPt.listKPt):
        f.write("%-10s%5i%5i%5i%5i%5.1f" % tuple( [line[0]] + [int(round(x*nDivisor)) for x in line[1:4]] + [nDivisor,2.0] ))
        if ( i != 0):
            f.write("\n")
        else: #This energy is useless
            f.write("-8.00 8.00       \n")
    f.write("END\n")
    f.close()

def w2k_ReadFermi():
    '''
    Read fermi energy from Wien2K folder in eV unit.
    The order of files we look for the fermi energy: scf2,qtl,scf
    '''
    stCaseName = f_GetCaseName()
    dFermi = None
    list_suffix = [None, "scf", "scf2", "qtl"]
    for i in range(1,4):
        filename = "%s.%s" % (stCaseName, list_suffix[i])
        if (not os.path.exists(filename)):
            continue
        if ( i == 1):
            stLine = io_grep_lines(filename, ":FER  :",-1,-1)
        elif ( i == 2):
            stLine = io_grep_lines(filename, ":FER  :",1,-1)
        elif ( i == 3):
            stLine = io_grep_lines(filename,"FERMI ENERGY",1,-1)

        if ( stLine != None):
            dFermi = float(stLine)
            break

    if ( dFermi == None):
        raise ValueError,"Cannot find Fermi energy in current folder"

    return dFermi*Ry2eV

def w2k_ReadElectron():
    '''
    Read number of electron from Wien2K folder from scf file
    '''
    filename =f_GetCaseName()+".scf"
    if (not os.path.exists(filename)):
       raise IOError("Cannot find file %s" % filename)
    return  int(float(io_grep_lines(filename,"NOE  :",-1,7)))

def w2k_ReadBand(suffix_energyfile="-band"):
    '''
    Return k-point and band structure  WIEN2k (must be under case folder)
    Number of electrons is calculated from band structure and fermi level * 2, so be careful with it(Now read from scf file)
    @todo Spin is always set to 1
    @todo Fermi energy is always DFT
    '''
    stCaseName = f_GetCaseName()
    #Read kpoints
#No suffix: case.energy and case.klist
#Any suffix: case.energy-$suffix and case.klist_band
    b_noband = suffix_energyfile == "" 
    suffix_kpt = ".klist" if b_noband else ".klist_band"
    aKPT = w2k_GetKPointList(stCaseName+".struct",stCaseName+suffix_kpt)
    #Remove names if no band
    if (b_noband):
        for kpt in aKPT.listKPt:
            kpt[0] = ""
    #Read band structure
    filename_energy = stCaseName+".energy"+suffix_energyfile 
    f = open(filename_energy)
    ar_stLine = f.readlines()
    if (len(ar_stLine) < 2):
        raise IOError("The energy file %s is empty" % filename_energy)
    #read energy one by one
    nMinStates = -1
    listBand = []
#find first k-point , in which line there is a "E" for floating point
    for i in range(0,100):
        if ( "E" in ar_stLine[i] ) :
            break
        
    while ( i < len(ar_stLine)):
        stLine = ar_stLine[i]
        kpt = [ float(stLine[0:19]),float(stLine[19:38]),float(stLine[38:57])]
        n = int(ar_stLine[i].split()[-2]) #energy in this k-points
        #Determine the minimum number of states and ignore those band not complete
        if ( nMinStates == -1 or n < nMinStates ):
            nMinStates = n

        listBandAtK = []
        for j in range(0,n):
            listBandAtK.append(float(ar_stLine[i+1+j].split()[1]))
        listBand.append( [kpt,listBandAtK] )
        i += n + 1
    f.close()

    #Write band
    #aKPT.WriteToFile(os.path.join(stOutputFolder,"band.klist"))

    #Convert band unit from Ry to eV
    listBand2 = []
    for aBand in listBand:
        listBand2.append( [ x*Ry2eV for x in aBand[1][:nMinStates]])

    #Store in unit eV
    #f = open(os.path.join(stOutputFolder,"band_ref.dat"),'w')
    #for i in range(0,len(listBand)):
    #    f.write(" ".join(["%18.14f" % (x*Ry2eV) for x in listBand[i][1][:nMinStates]]))
    #    f.write("\n")

    f.close()

    #deal with fermi energy
    dFermi = w2k_ReadFermi()
#   dGap,dVBM2,nHOIndex = f_Band_GetGap(listBand2,dFermi,False)
#Read number of electron 
    num_electron = w2k_ReadElectron()

    band = BandsT(aKPT,listBand2,None,None,num_electron,1)
    if (not band.guess_vbm()):
        band.fermi = dFermi

    return band
#    return aKPT,listBand2,dVBM2,nHOIndex*2+2,1

def w2k_ReadPDOS(bCustomName=False):
    '''
    Read PDOS data from dos*eV in wien2k case directory
    The PDOS name can be infered from case.int/case.qtl/case.struct or directly read from dos data
    @todo deal with spin
    :param bCustomName: use name in dos*eV (which is specified by user), default is use system-detected names
    ''' 
    listFile = os.listdir(os.getcwd());
    stCaseName = f_GetCaseName()
    reg1 = re.compile(r"%s\.dos(.+)ev" % stCaseName)

    listFile2 = [None for x in xrange(0,99)]
#Search all dos*eV files and sort them
    for stFileName in listFile:
        aMatch = reg1.match(stFileName)
        if ( aMatch != None):
            listFile2[int(aMatch.group(1))-1] = stFileName
#Reading all files
#Note "# Energy" and first col should be deleted
    ar1 = []
    ar2 = []
    listColName = []
    listEnergy = None
    for stFileName in listFile2:
        if ( stFileName == None):
            break
        f = open(stFileName)
        ar2 = [x.strip() for x in f.readlines()]
        f.close()
#Read energy and VBM only once
        ar3 = [[float(y) for y in x.split()] for x in ar2[3:]]
        if (listEnergy is None):
            dVBM2 = float(ar2[1][4:14])
            listEnergy = [x[0] for x in ar3]
        ar3 = [x[1:] for x in ar3]
        listColName += ar2[2].split()[2:]
        ar1.append(ar3)

#Read DOS
    listDOS = list(zip(*[reduce(list.__add__,x) for x in zip(*ar1)]))

#Read Fermi and DOS from qtl
    f = open(stCaseName+".qtl")
    ar_stLine = f.readlines()
    f.close()
    dVBM = float(ar_stLine[2][56:])
    dVBM = dVBM2 #Always use fermi from .dos*ev
    nAtom = int(ar_stLine[3][35:38])
    listOrb = [x[31:].strip().split(',') for x in ar_stLine[4:4+nAtom]]
    listOrb = [["tot"]] + listOrb #Totol dos index 0

#Read name
    if (bCustomName): #Read from .dos*eV
        listName = [ OrbitalNameT(name=x) for x in listColName]
    else: #Read from case.int and case.struct
        f=  open(stCaseName+".int")
        ar_stLine = f.readlines()
        nCol = int(ar_stLine[2].split()[0])
        ar2 = [[ int(y) for y in x.split()[:2]] for x in ar_stLine[3:3+nCol]]
        #Read atom name
        aStruct = Wien2K_Structure(stCaseName+".struct")
    #Use atom name from struct file, only the first column
        listName = []
        for orb in ar2:
            ont = OrbitalNameT()
            if ( orb[0] == 0):
                ont.total = True
            else:
                stOrbName = listOrb[orb[0]][orb[1]-1]
                stAtom = aStruct.listAtom[orb[0]-1].stName.split()[0]
                if ( stOrbName.isdigit()): #If it is a number, then it is the quantum number l
                    ont.l = int(stOrbName)
                    ont.i = orb[0]
                    ont.species = stAtom
                else:#Custom name
                    if ( stOrbName == "tot"):#tot of atoms
                        ont.i = orb[0]
                        ont.species = stAtom
                    else:
                        ont.name = "%s(%s)" % ( stAtom,stOrbName)
            listName.append(ont)

    return listName,listEnergy,listDOS,dVBM

def w2k_read_qtl(filename=None):
    '''
    Read casename.qtl to a 2-D array [band,kpt] format
    :return: the first is list of components name, the second is a nested list in [band,kpt], each element is [energy,components1,2,3...]
    '''
    if ( filename is None or filename == "" ):
        casename = f_GetCaseName()
        filename = casename + ".qtl"

    f = open(filename)
    lines = f.readlines()
    f.close()

    list_comp = []
    list_orb = []
    ix_atom = 0
    ln_start = -1
#Read name
    for line in lines:
        ln_start += 1
        if ("JATOM" in line):
            ix_atom += 1
            list_orb_atom = line[31:].split(",")
            list_orb += [(ix_atom,x.strip()) for x in list_orb_atom]
        if ("BAND" in line):
            break
    list_orb.append((-1,"inter"))

    num_atom = ix_atom + 1
#Read components
    list_band = None
    i = ln_start-1
    while (i<len(lines)-1):
        i += 1
        line = lines[i]
        if ("BAND" in line):
            if (list_band is not None):
                list_comp.append(list_band)
            list_band = []
        else:
            list_kpt = []
            list_kpt.append(float(line[:10]))
            for j in range(0,num_atom):
                list_kpt += [float(x) for x in lines[i+j][14:].split()]
            list_band.append(list_kpt)
            i += num_atom - 1
    if (list_band is not None):
        list_comp.append(list_band)
    
    return list_orb,list_comp
        
def w2k_read_band_character(filename=None,filename_struct=None):
    '''
    Read band characters from qtl
    '''
    list_orb, list_comp = w2k_read_qtl(filename)
    struct = Wien2K_Structure(filename_struct if filename_struct is not None else f_GetCaseName() + ".struct")
#Convert into common format
    l_last = -1
    list_orb2 = []
    for orb in list_orb:
        orb1 = OrbitalNameT()
        if (orb[0] == -1): #Interstetial region
            orb1.name = "inter"
            continue
        orb1.species = struct.listAtom[orb[0] - 1].stName.split()[0]
        orb1.i = orb[0]
        if (orb[1].isdigit()):
            orb1.l = int(orb[1])
        elif (orb[1] == "tot"): #Summation of all l (not exactly same)
            pass
        else: #Custom m
            orb1.l = OrbitalNameT.dic_l_rev[orb[1][0].lower()]
            orb1.m_name = orb[1]
        list_orb2.append(orb1)

#Sort and remove the first element (energy) from list
    list_comp = [ [y[1:] for y in x] for x in zip(*list_comp)]
        
    return list_orb2, list_comp

def w2k_rescale_character(list_orb, list_comp):
    '''
    w2k projections are generally much smaller than the total value, like the interstitial region
    This function rescale all values to summation of all $l$ as 1
    '''
#Find the new normalization factors for l=0,1,2,3...
    list_include = [i for i,orb in enumerate(list_orb)
            if orb.l != OrbitalNameT.qn_none and orb.m == OrbitalNameT.qn_none]

    list_comp2 = []
    for list_kp in list_comp:
        l2 = []
        for list_cha in list_kp:
            norm = sum([list_cha[ix] for ix in list_include])
            l2.append([x/norm for x in list_cha])
        list_comp2.append(l2)
    
    return list_comp2


def Main(ArgList):
    usage = "usage: w2k_caseutil [-h]  [-s] [-p] [--jobsystem JOBSYSTEM]\n \
    this script can be called to show summary of WIEN2k calculation with -s\n \
    this scirpt can be called to generate .machines files for parallel execution with -p \
    --jobsystem can be used to specify Job system"
    parser = OptionParser(usage)
    #parser.add_option("-c",dest="stCompound",help="The compound name")
    parser.add_option("-t",dest="bTest",action="store_true",default=False,help="Test ONLY")
    parser.add_option("-s",dest="bShow",action="store_true",default=False,help="Show current Wien2K case calculation parameter (Ignore all other options if used)")
    parser.add_option("-v",dest="bVersion",action="store_true",default=False,help="Show python version")
    parser.add_option("-p","--para",dest="bPara",action="store_true",default=False,help="Generate .machines files from host string ( PBS or LSF )")
    parser.add_option("--jobsystem",dest="stJobSystem",help="Set Job Management System type (PBS or LSF)")
    (options,args) = parser.parse_args(ArgList)
    if ( len(args) != 1 ):
        parser.error("incorrect number of arguments.")
    
    if(options.bTest):
        #w2k_Modify_in1()
        #w2k_ConvertWien2KStructureToMadelung(options.stCompound)
        return
    
    #if (options.stCompound == None):
    #    options.stCompound = f_GetCaseName()
    
    if ( options.bVersion):
        print(sys.version)
        return
 
    if ( options.bPara ):
        if ( options.stJobSystem == None ):
            stJobSystem = f_env_GetJobSystem()
        else:
            stJobSystem = options.stJobSystem
        Wien2K_machines(stJobSystem).WriteToFile()
        return
    
    if(options.bShow):
        aCase = Wien2K_Case("",options.stCompound)
        aCase.ShowSummary()
        return
        
if __name__ == "__main__":
    #print(sys.version)
    Main(sys.argv)
