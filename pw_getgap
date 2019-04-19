#!/usr/bin/env python

import re
import sys
import os
import string
import copy
import re
import subprocess
from optparse import OptionParser
from common_caseutil import f_File_IsNullOrEmpty,f_ReadStdout,f_GetCaseName
from w2k_caseutil import f_w2k_GetW2kDir,f_w2k_GetW2kFilePath
from constants import Ry2eV

# Get Gap function
class kPointBand:
    def __init__(self,file=None): 
        '''
        Init
        '''
        self.list_fVector=[0,0,0]
        self.stName=""
        self.nBase=0
        self.nNumber=0
        self.fWeight=1.0
        self.list_fEnergy=[]
        self.bError = False
        if ( file ):
            self.bError = self.ReadFromFile(file)
    
    def GetkVectorString(self):
        return "%.3f  %.3f  %.3f" % (self.list_fVector[0],self.list_fVector[1],self.list_fVector[2]) 
    
    def ShiftEnergy(self, fFermiEnergy):
        '''
        Shift energy by fermi energy as reference 
        '''
        self.list_fEnergy = [x-fFermiEnergy for x in self.list_fEnergy]
        return self
    
    
    def ReadFromFile(self,file):
        '''
        Read info from case.energy file, convert 2Ry to eV automatically 
        '''
        stPattern=r"(.*E...) *(.*E...) *(.*E.\d\d)(.{12}) *([^ ]*) *([^ ]*) *([^ ]*)"
        stLine = file.readline()
        bError = True
        # read head line
        if ( stLine != "" ):
            match = re.match(stPattern,stLine)
            if ( match ):
                self.list_fVector[0] = float(match.group(1))
                self.list_fVector[1] = float(match.group(2))
                self.list_fVector[2] = float(match.group(3))
                self.stName = match.group(4)
                self.nBase = int(match.group(5))
                self.nNumber = int(match.group(6))
                self.fWeight = float(match.group(7))
                
                # read band energy lines
                stPattern=r" *([^ ]*) +([^ ]*) *"
                for i in range(0,self.nNumber):
                    stLine=file.readline()
                    match = re.match(stPattern,stLine)
                    if ( match ):
                        self.list_fEnergy.append(float(match.group(2))*Ry2eV)                       
                bError = False
        return bError


def f_GetNonequalAtomCount(stCasename):
    return int(f_ReadStdout("head -2 %s.struct |tail -1 |cut -c28-30" % stCasename))


def f_BandCheck(nIndex , nMaxIndex, stName):
    '''
    Check if nIndex is larger than Max index and show error message. Return True if check passed, vice versa.
    '''
    bPass = True
    if ( nIndex < 0 ):
        print("ERROR: "+stName+" < 1  !!!")
        print("  --- Check the Fermi energy")
        bPass = False
    if ( nIndex > nMaxIndex ):
        print("ERROR: "+ stName +" is larger than the number of available bands  !!!")
        print("  --- Check the Fermi energy")
        bPass = False
    
    return bPass


def f_BandAnalyse(list_kPointBandOrg,fFermiEnergyInRy):
    '''
    Analyse band structure in specific energy range , with given Fermi energy 
    '''
    # get min band count , min/max energy in kPBand
    list_kPointBand = copy.deepcopy(list_kPointBandOrg)
    fFermiEnergy = fFermiEnergyInRy*Ry2eV
   
    nBandCountMin = min([ x.nNumber for x in list_kPointBand ])
    fEnergyMin = min([min(x.list_fEnergy) for x in list_kPointBand])
    fEnergyMax = max([max(x.list_fEnergy) for x in list_kPointBand])
    
    print("----------------------------------------------------")
    print("         Band Structure Analysis                    ")
    print("----------------------------------------------------")
    print('Range of bands considered: 1 ' + str(nBandCountMin))
    
    if ( fEnergyMin > fFermiEnergy or fEnergyMax < fFermiEnergy ):
        print("WARNING from bandanaly: ")
        print(" - Fermi energy outside the energy range of bande!")
        print("  Fermi Energy(eV):   " + str(fFermiEnergy*Ry2eV))
        print(" - Maximal energy:    "+str(fEnergyMax))
        print(" - Minimal energy:    "+str(fEnergyMin))
        return
    
    # find HO and LU . Index is the position in array
    nHOIndex = -1
    nLUIndex = 1000
    
    # note: in Fortran array starts from 0 but in Python from 1
    for aBand in list_kPointBand:
        #Sk mean single k point
        nSkHOIndex = -1
        nSkLUIndex = 1000
        for i in range(0,nBandCountMin):
            if ( aBand.list_fEnergy[i] < fFermiEnergy ):
                nHOIndex = i if nHOIndex < i else nHOIndex
            else:
                nLUIndex = i if nLUIndex > i else nLUIndex
    
    # error checking
   
    if  not ( f_BandCheck(nHOIndex,nBandCountMin,"HO Band") and f_BandCheck(nLUIndex,nBandCountMin,"LU Band") ):
        return
    
    if ( nHOIndex > nLUIndex ):
        print("WARNING: Valence and Conductance bands overlap !!")
        print("Swap HO Band and LU Band...")
        nLUIndex,nHOIndex = nHOIndex,nLUIndex
    else:
# Sometimes when the Fermi energy is not calculated very accurately, judging
# whether it is conductance or valence band only according to its sign is not reliable
# test which dominate in k space and change the other
        if ( nHOIndex == nLUIndex ):
            print("WARNING: Valence and Conductance bands identical !!")
            nCount = 0
            for aBand in list_kPointBand:
                if ( aBand.list_fEnergy[nHOIndex] > fFermiEnergy ):
                    nCount += 1
            if (nCount >= len(list_kPointBand) /2):
                nHOIndex -= 1
            else:
                nLUIndex += 1
    
    print('Highest occupied band: '+str(nHOIndex+1))
    print('Lowest unoccupied band:'+str(nLUIndex+1))
    
    # from HO and LU , find energy and gap
    for aPointBand in list_kPointBand:
        aPointBand.ShiftEnergy(fFermiEnergy*Ry2eV)
    
    
    list_HOBand= [ x.list_fEnergy[nHOIndex] for x in list_kPointBand ]
    list_LUBand= [ x.list_fEnergy[nLUIndex] for x in list_kPointBand ]
    fHOEnergy = max(list_HOBand)
    nHOkPointIndex = list_HOBand.index(fHOEnergy)
    fLUEnergy = min(list_LUBand)
    nLUkPointIndex = list_LUBand.index(fLUEnergy)
    
    print(':BandGap =  %.4f eV' % (fLUEnergy-fHOEnergy))
    if ( nHOkPointIndex <> nLUkPointIndex):
        print("Indirect gap:")
        # also print direct gap at HO and LU k point
        print("\t%.4f eV at VBM k= %s, ik= %d"% (list_kPointBand[nHOkPointIndex].list_fEnergy[nLUIndex]-fHOEnergy,list_kPointBand[nHOkPointIndex].GetkVectorString(),nHOkPointIndex+1))
        print("\t%.4f eV at CBM k= %s, ik= %d"% (0-list_kPointBand[nLUkPointIndex].list_fEnergy[nHOIndex]+fLUEnergy,list_kPointBand[nLUkPointIndex].GetkVectorString(),nLUkPointIndex+1))
        #print(str(list_kPointBand[nHOkPointIndex].list_fEnergy[nLUIndex]-fHOEnergy)+ " eV at VBM k=" + list_kPointBand[nHOkPointIndex].GetkVectorString() +  ' ik=' + str(nHOkPointIndex))
        #print(str(0-list_kPointBand[nLUkPointIndex].list_fEnergy[nHOIndex]+fLUEnergy)+ " eV at CBM k=" + list_kPointBand[nLUkPointIndex].GetkVectorString() +  ' ik=' + str(nLUkPointIndex))
    else:
        print('Direct gap at k=' + list_kPointBand[nLUkPointIndex].GetkVectorString() + ' ik=' + str(nLUkPointIndex))
    
    
    print(":Range of each band (with respect to VBM):")
    print("    n     Bottom       Top         Width")
    for i in range(0,nBandCountMin):
        aBand =  [x.list_fEnergy[i] for x in list_kPointBand]
        fEnergyMax = max(aBand)-fHOEnergy
        fEnergyMin = min(aBand)-fHOEnergy
        print("%5i%12.4f%12.4f%12.4f" % ( i+1, fEnergyMin, fEnergyMax, fEnergyMax-fEnergyMin) )
    
    


def f_GetGap(stFileName,nNumber,fEnergy):
    '''
    Get gap from calculation resuls
    '''
    try:
        fFileEnergy= open(stFileName,"r")
    except IOError:
        print("ERROR: Fail to open the energy file " + stFileName + ".")
        return 0
    
    # skip atom description ( 2 line for each atom)
    for i in range(0,nNumber*2,1):
        fFileEnergy.readline()
        
    
    list_kPointBand = []
        
    while True:
        akPointBand = kPointBand(fFileEnergy)
        if ( akPointBand.bError ):
            break
        else:
            list_kPointBand.append(akPointBand)
        
            
    f_BandAnalyse(list_kPointBand,fEnergy+0.0001)
  


def f_ReadStructType(stFileName):
    '''
    Read struct lattice type from .struct
    '''
    return f_ReadStdout("head -n 2 %s | tail -n 1" % stFileName)[0]

def f_GetFermiEnergyGW(stCasename,stGWDir):
    return float(f_ReadStdout("grep :QP_FERMI_ENERGY %s/%s.outgw | tail -1 | awk '{print $2}'" % (stGWDir,stCasename)))


def f_GetFermiEnergy(stCasename,nFermiFile,nSpin,bGW = False):
    '''
    Read Fermi Energy from file.
    nFermiFile =
    0  -- from *.output2
    1  -- from *.qtl
    2  -- from *.dos
    3  -- from *.scf 
    nSpin =
    0 -- spin-unpolarized
    1 -- spin up
    2 -- spin and down
    '''
    list_stCommand = [
        "grep ':FER' %s.output2%s | cut -c40-48",
        "head -3 %s.qtl%s | tail -1  | cut -c57-66",
        "grep \"#EF\" %s.outputt%s | tail -n 1 | awk '{print $2}'",
        "grep \":FER\" %s.scf | tail -n 1 | awk '{print $NF}'"
    ]
    fFermiEnergy = 0
    
    stSpin = ("" if nSpin == 0 else "up")
    
    if ( nFermiFile == 0):
        print("  Extract Fermi energy from case.output2")
        fFermiEnergy = float(f_ReadStdout(list_stCommand[0] % (stCasename,stSpin)))
    if (nFermiFile == 1):
        print("  Extract Fermi energy from case.qtl")
        fFermiEnergy = float(f_ReadStdout(list_stCommand[1] % (stCasename,stSpin)))
    if (nFermiFile == 2):
        print("  Extract Fermi energy from DOS output")
        fFermiEnergy = float(f_ReadStdout(list_stCommand[2] % (stCasename,stSpin)))
    if (nFermiFile == 3):
        print("  Extract Fermi energy from SCF file")
        fFermiEnergy = float(f_ReadStdout(list_stCommand[3] % (stCasename)))
    
    return fFermiEnergy

def f_Run_GetGap(options):
    '''
    Get band gap from Wien2K files.
    Return 1 if there're any error occured, 0 for OK.
    '''
    stCasename = f_GetCaseName()
    if ( options.casename ):
        stCasename = options.casename
    print("Case name: %s" % stCasename)
    
    # check which file to parse
    if (options.bandfile == 0):
        print("Extract band gap from %s.scf " % stCasename)
        if( not ( options.vbm and options.nb )):
            print('VBM or nb not defined!')
            return 1
        print("Valence band maximum=  %f" % options.vbm)
        print( "Total number of bands= %d" % options.nb)
        fEnergyVBMax = f_ReadStdout("grep 'BAND   :' $file.scf | tail -n $nb |  awk -v vb=$vbm ' $3==vb   {print $5}'")
        fEnergyCBMin = f_ReadStdout("grep 'BAND   :' $file.scf | tail -n $nb |  awk -v vb=$vbm ' $3==vb+1 {print $4}'")
        print(":Egap= %f" % ((fEnergyCBMin-fEnergyVBMax)*Ry2eV) )
        return 0
    
    nNEAtomCount= int(f_GetNonequalAtomCount(stCasename))
    print("  Number of  nonequivalent atoms %d" % nNEAtomCount)
    
    # get fermi energy
    
    #fFermiEnergy = 0
    #if ( options.fermi == 0):
    #    print("  Extract Fermi energy from case.output2")
    #    stCommand = "grep ':FER' %s.output2%s | cut -c40-48" % ( stCasename,  "" if options.spin == 0 else "up"   )
    #    fFermiEnergy = float(f_ReadStdout(stCommand))
    #if (options.fermi == 1):
    #    print("  Extract Fermi energy from case.qtl")
    #    stCommand = "head -3 %s.qtl%s | tail -1  | cut -c57-66" % ( stCasename, "" if options.spin == 0 else "up"   )
    #    fFermiEnergy = float(f_ReadStdout(stCommand))
    #if (options.fermi == 2):
    #    print("  Extract Fermi energy from DOS output")
    #    stCommand = "grep \"#EF\" %s.outputt%s | tail -n 1 | awk '{print $2}'" % ( stCasename, "" if options.spin == 0 else "up"   )
    #    fFermiEnergy = float(f_ReadStdout(stCommand))
    
    fFermiEnergy = f_GetFermiEnergy(stCasename,options.fermi,options.spin)
    print('  Fermi Energy= %f' % fFermiEnergy)
   
    #For compability with pw_band and w2k_band results ( with "-band" suffix), it will automatically detect whether "-band" file exists if energy file is not present.

    if ( not options.so ):
        list_stSpin = [ [""],["up"],["up","dn"] ]
        for stSpin in list_stSpin[options.spin]:
            stEnergyFileName = "%s.energy%s" % ( stCasename,stSpin)
            if ( f_File_IsNullOrEmpty(stEnergyFileName)):
                print("%s not found, try -band files"%stEnergyFileName)
                stEnergyFileName += "-band"
                if ( f_File_IsNullOrEmpty(stEnergyFileName)):
                    raise RuntimeError,"Energy file not found!"
            f_GetGap(stEnergyFileName,nNEAtomCount,fFermiEnergy)
    else:
        stEnergyFileName = "%s.energy%s" % ( stCasename,"so")
        if ( f_File_IsNullOrEmpty(stEnergyFileName)):
            print("%s not found, try -band files"%stEnergyFileName)
            stEnergyFileName += "-band"
            if ( f_File_IsNullOrEmpty(stEnergyFileName)):
                raise RuntimeError,"Energy file not found!"
        f_GetGap(stEnergyFileName,nNEAtomCount,fFermiEnergy)
    
    return 0
    


def Main_w2k_getgap(ArgList):
    usage = "%prog [options]  <nl for each l at each atom>"
    parser = OptionParser(usage)
    parser.add_option("-f",dest="casename",help="Wien2K casename ( needed only when DIR name is not casename ) ")
    parser.add_option("--fermi",type="int",dest="fermi",default=0,help="option for where to extract Fermi energy.\n0  -- from output2\n1  -- from *.qtl\n2  -- from *.dos\n")
    parser.add_option("-o",type="int",dest="bandfile",default=1,help="determine where to get band gap\n0 -- from *.scf ( -vbm and -nb must be present)\n1 -- from *.energy")
    parser.add_option("--vbm",type="float",dest="vbm",help="valence band maximum")
    parser.add_option("--nb",type="int",dest="nb",help="number of bands in *.scf")
    parser.add_option("-s",type="int",dest="spin",default=0,help="option for spin polarized cases.\n0 -- spin-unpolarized\n1 -- spin up\n2 -- spin and down\n")
    parser.add_option("--so",action="store_true",dest="so",default=False,help="for spin-orbit coupling")
    
    (options,args) = parser.parse_args(ArgList)
    if ( len(args) > 1    ):
        parser.print_help()
        parser.error("There are unknown options.")
        
    result = f_Run_GetGap(options)
    
    if ( result <> 0):
        print("ERROR in w2k_getgap")
        
    return result


if __name__ == "__main__":
    Main_w2k_getgap(sys.argv)

