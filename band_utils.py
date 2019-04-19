#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import sys,os,copy,re,commands
import subprocess

import common_caseutil
from py_utils import EmptyFormatter
from common_caseutil import debug_print,Lattice,KPointsT,f_WriteFileFloatValue,f_WriteFileIntValue,f_ReadFileIntValue,f_ReadFileFloatValue,f_Data_ReadMultiCol,f_Data_WriteMultiCol
from chem_utils import f_Element_Z_to_ExtraValenceElectronCount,f_Element_Symbol_to_Z
from py_xmlio import XMLSerilizer
import list_utils as lu
import itertools
import stat
from py_modagr import set_property_auto, split_multiline_to_single

common_caseutil.logLevel = 0
#logLevel = 3

class BandSpinT():
    '''
    Class represent eigenvalues, occupation numbers and orbitals for each band
    '''
    def __init__(self, eig=None, occ=None, orb=None, character=None):
        self.eig = eig #: eigenvalues in [kpt, band index]
        self.occ = occ #: occupation numbers, same above
        self.orb = orb #: list of OrbitalNameT 
        self.character = character #: orbital projection coefficients in [kpt, band index, orbs]

    @property
    def num_band_common(self):
        '''
        Minimum of number of bands in all k-points
        '''
        if (len(self.eig) == 0):
            return 0
        else:
            return min([len(x) for x in self.eig])

    def resolve_crossing(self, list_turnkpt):
        '''
        Try to resolve band crossing by derivatives

        :param list_turnkpt: a list of k-point indicies where treated as turning point
        '''
#Convert occ and character also if exists
        list_aux = []
        dic_aux = {}
        for term in ("occ", "character"):
            if (hasattr(self, term)):
                v = getattr(self, term)
                if (v is not None):
                    list_aux.append(v)
                    dic_aux[len(list_aux)-1] = term

#Resolve it
        list_result, list_aux2 = f_band_resolve_crossing(
                self.eig, 
                func_comp = lambda x,y:x-y,func_eig=lambda x:x,
                list_turnkpt = list_turnkpt,
                list_aux = list_aux)

#Store
        self.eig = list_result
        for key, val in dic_aux.items():
            setattr(self, val, list_aux2[key])

        return



class BandsT():
    '''
    Class represent band structure information
    There are two way to read, BandsT.load_xml and load() (the same as BandsT(dirname)); But the second one is deprecated and not supported
    '''
    align_none = 0
    align_vbm = 1
    align_fermi = 2
    align_auto = 3

    def __init__(self,kpt=None,list_eig=None,vbm=None,fermi=None,num_electron=None,num_spin=1,b_metal=None,gap=None,dirname=None,list_occ=None):
        '''
        Init with possible content

        @todo spin-polarized does not works well with occupation numbers and characters

        :param dirname: if this parameter is specified then try to read from this folder
        :param list_eig: a [kpt,band] 2-D array contains information if num_spin = 1; a list of two arrays if num_spin = 2
        :param list_occ: a [kpt,band] 2-D array contains occupation numbers if num_spin = 1; a list of two arrays if num_spin = 2
        :param num_spin: maybe 1 or 2, default 1
        '''
        self.version = 2
        self.kpt = kpt

#Pay attention that eig, orb, character are saved for each spin
#list_eig is the brute force combination of different spins
#list_orb and list_character should never be used

        if (num_spin == 1):
            self.prop = [BandSpinT(list_eig, list_occ)] #This should be 
            self.list_eig = list_eig
        elif (num_spin == 2): #Spin polarized
            if (list_occ is None):
                list_occ = [None, None]
            self.prop = [BandSpinT(eig, occ) for eig, occ in zip(list_eig, list_occ)]
            self.list_eig = [sorted(up+down) for up,down in zip(*list_eig)]
        else:
            raise ValueError("The number of spin must be 1 or 2 : %s" % (str(num_spin)))

        self.vbm = vbm
        self.fermi = fermi
        self.num_electron = num_electron
        self.num_spin = num_spin
        self.b_metal = b_metal 
        self.gap = gap

        if (dirname is not None):
            self.load(dirname)

    @property
    def num_kpt(self):
        '''
        Number of k-points
        '''
        return len(self.list_eig)


    @property
    def eig_max(self):
        return max([max(x) for x in self.list_eig])

    @property
    def eig_min(self):
        return min([min(x) for x in self.list_eig])

    def save_xml(self,filename_out="band.xml"):
        '''
        Save band structure information into a XML file
        '''
        xp = XMLSerilizer(globals(),"band_info2.xsd")
        xp.SerilizeToFile(self,filename_out)

    @classmethod
    def load_xml(cls,filename_in="band.xml"):
        '''
        Load band structure information from a XML file
        '''
        xp = XMLSerilizer(globals(),"band_info2.xsd")
        obj = xp.Deserilize(filename_in)
#If no version information, use version 1 and convert
        if (not hasattr(obj, "version")):
            xp = XMLSerilizer(globals(),"band_info.xsd")
            obj = xp.Deserilize(filename_in)
#Convert to version 2
            if ( (not hasattr(obj, "num_spin") and not hasattr(obj, "list_eig1"))  or obj.num_spin == 1):
                obj.prop = [BandSpinT(obj.list_eig)]
            elif (obj.num_spin == 2):
                obj.prop = [BandSpinT(obj.list_eig1), BandSpinT(obj.list_eig2)]


        return obj 

    def save(self,dirname_out="band",name_prefix="band",b_plot=True,
            align=align_none, plotinfo=None):
        '''
        Save band related data to specific folder
        :param dirname_out: the directory to save generated files
        :param name_prefix: all data will be saved in prefix.xxx filenames
        :param b_plot: whether to create xmgrace files for plot
        :param b_align: When plotting, whether to align eigenvalues to fermi level (for metal) or VBM (for insulator). The value must be BandsT.align_ auto,fermi,vbm or none
        '''
        #Guess VBM
#       vbm_old = self.vbm
        self.guess_vbm()
#       print("Guess VBM:%f, old VBM:%f" % (self.vbm, vbm_old))

#Try to read the fermi or VBM as reference energy
        if (not self.b_metal):
            print("Treat as an insulator...")
            if (self.vbm is None):
                if (self.fermi is None):
                    print("Cannot guess the VBM...")
                    eig_ref = 0.0
                else:
                    eig_ref = self.fermi
            else:
                eig_ref = self.vbm
        else:
            print("Treat as a metal...")
            if (self.fermi is None):
                if (self.vbm is None):
                    pass
                else:
                    eig_ref = self.vbm
            else:
                eig_ref = self.fermi

#Determine whether we will align the data
        eig_align = None
        if (align == BandsT.align_vbm):
            if (self.vbm is None):
                raise ValueError("Cannot aligned to VBM because of no VBM")
            eig_align = self.vbm
        elif (align == BandsT.align_fermi):
            if (self.fermi is None):
                raise ValueError("Cannot aligned to Fermi energy because of no Fermi energy")
        elif (align == BandsT.align_auto):#Choose whatever possible, VBM for insulator first and Fermi for metal first
            eig_align = eig_ref
        else:
            pass

#       f_Band_SaveData(self.kpt,self.list_eig,eig_ref,dirname_out,self.num_electron,b_plot,eig_ref is None)
        if ( dirname_out == None):
                dirname_out = os.getcwd()

        if ( not os.path.exists(dirname_out) ):
           os.mkdir(dirname_out)
       
        dirname_cur = os.getcwd()
        os.chdir(dirname_out)

        f_Band_WriteToFile(self.list_eig,name_prefix+".dat")
        self.kpt.WriteToFile(name_prefix+".klist",bWriteName=True)
#       #If number of electrons is specified, then calculate VBM and use it
#       if ( nElectron != 0 and nElectron != None):
#           dVBM2 = f_Band_GetVBMFromElectron(listBand,nElectron)
#           if ( dVBM2 != None and dVBM2 != 0):
#               print("Detect VBM: %f" % dVBM2)
#               dVBM = dVBM2

        #Save VBM and electron count
        if (self.b_metal is not None):#If undetermined then just skip
            #Calculate gap for insulator
            if (self.b_metal):
                f_WriteFileFloatValue(self.fermi,name_prefix+".fermi")
            if (not self.b_metal):
                if (self.fermi is not None):#Only store fermi energy if available
                    f_WriteFileFloatValue(self.fermi,name_prefix+".fermi")
                f_WriteFileFloatValue(self.vbm,name_prefix+".vbm")
                self.calc_gap()
                if ( self.gap is not None):
                    f_WriteFileFloatValue(self.gap,name_prefix+".gap")
            if ( self.num_electron != None):
                f_WriteFileIntValue(self.num_electron,name_prefix+".electron")
            pass

    #Save xml
        print("Save XML...")
        self.save_xml(name_prefix+".xml")

        if (b_plot):
            eig_fermi_plot = eig_ref
            if (eig_align is not None):
                print("Plot eigenvalues with %8.4f as the reference..." % eig_ref)
#Where shall we plot the fermi energy?
                eig_fermi_plot = eig_ref - eig_align

            plot_setting = None if plotinfo is None else plotinfo.get("plot_setting")
            f_Band_Plot([self],[eig_align],y_fermi=eig_fermi_plot, plot_setting=plot_setting)

            if (self.prop[0].orb is not None and plotinfo is not None and plotinfo.has_key("orb_map")):
                f_band_plot_character(self, eig_align, y_fermi=eig_fermi_plot,
                        **plotinfo)

        os.chdir(dirname_cur)

    def load(self,dirname,name_prefix="band"):
        '''
        Load all files in specific folder as band structure information saved by save()
        Note, prefix.fermi maybe vbm in read.
        '''
        name1 = os.path.join(dirname,name_prefix)
        self.list_eig = f_Band_ReadFromFile(name1+".dat")
        self.kpt = KPointsT(name1+".klist")
        self.num_electron = f_ReadFileIntValue(name1+".electron")
        self.gap = f_ReadFileFloatValue(name1+".gap")
        self.num_spin = f_ReadFileIntValue(name1+".spin")
        self.fermi = f_ReadFileFloatValue(name1+".fermi")
        self.vbm = f_ReadFileFloatValue(name1+".vbm")
        self.b_metal = None
#Judge whether .fermi file is fermi or vbm energy
        if (self.gap is not None):
            if (self.gap > 0):#Insulator
                self.b_metal = False
                if (self.vbm is None and self.fermi is not None):
                    self.vbm = self.fermi
                    self.fermi = None

    def show_info(self, eig_ref = None):
        '''
        Print band structure basic information
        '''
        f_Band_Analyse(self.list_eig, eig_ref)

    def guess_vbm(self,b_print=False):
        '''
        Try to guess VBM from number of electrons if available. If it is a metal, try to guess the Fermi level.

        :param b_print: Display verbose informations
        :return: Return True if VBM is obtained, False if it is metal, None if it cannot be determined
        '''
        if (self.vbm is not None):
            print("VBM is known so no attempt to guess is made")
            return None
        if (self.num_electron is not None and self.num_electron != 0):
            vbm = f_Band_GetVBMFromElectron(self.list_eig,self.num_electron*self.num_spin,1)
            if (vbm is None):#Metal
                self.b_metal = True
                if (b_print):
                    print("This band structure is metallic!")
                if (self.fermi is not None):
                    print("Fermi is known so no attempt to guess is made")
                    return None
                else:
                    self.fermi = f_Band_GetFermiFromElectron(self.list_eig,[x[-1] for x in self.kpt.listKPt],self.num_electron*self.num_spin,1)
                return False
            else:
                self.vbm = vbm
                self.b_metal = False
                if (b_print):
                    print("Guess VBM: %f" % vbm)
                return True
        else:
            if (b_print):
                print("Cannot guess VBM without number of electrons")
            return None

    def resolve_crossing(self):
        '''
        Resolve the band crossing
        '''
        for band_spin in self.prop:
            band_spin.resolve_crossing(self.kpt.GetTurningPoint())

#Regenerate list_eig
        if (self.num_spin == 1):
            self.list_eig = self.prop[0].eig
        elif (self.num_spin == 2): #Spin polarized
            self.list_eig = [up + down for up,down in zip(self.prop[0].eig, self.prop[1].eig)]

        return
    
    def calc_gap(self):
        '''
        Calculate the gap and put it in self.gap var

        :return: gap if there is a gap, None otherwise
        '''
        if (not self.b_metal):
            if (self.vbm is None):
                self.guess_vbm()
            self.gap = f_Band_GetGap(self.list_eig,self.vbm,bPrint=False)[0]
        else:
            self.gap = None
        return self.gap


def f_Band_ReadFromFile(stFileBand):
    '''
    Read bands from files
    The file format is assumed to be space or tab seperated txt files, with eigenvalues of each k-points in one line.
    '''
    #Read band data
#   f = open(stFileBand,'r')
#   ar_stLine = f.readlines()
#   list_band = [ [ float(x) for x in y.split() ] for y in ar_stLine]
#   f.close()
#   return list_band 
    return f_Data_ReadMultiCol(stFileBand,func=float)

def f_Band_WriteToFile(list_band,stFileName):
    '''
    Write bands to file in format of space seperated format
    '''
    if ( os.path.exists(stFileName)):
        print("Warning: %s exists, overwritten!" % stFileName)
#   f = open(stFileName,'w')
#   for kpt in list_band:
#       for fE in kpt:
#           f.write("%12.7f " % fE)
#       f.write("\n")
#   f.close()
    f_Data_WriteMultiCol(list_band,stFileName,"",lambda x:"%12.7f " %x)

def f_Band_DiffDispersionRMSD(band1,band2,listBandWeight=None,dReverseWeight=1,b_debug=False):
    '''
    Return root-mean-square-deviation of two band structure, based on band dispersions
    If two bands contains different number of bands, then listBandWeight must be set

    :param band1: band data 1
    :param band2: band data 2
    :param listBandWeight:  Weight of each band. If not specified, all band will be treated as 1. 1-D array, start from first band in array.
    :param dReverseWeight: if dispersions are in different sign in band1 and band2, assign a large weight
    :param b_debug: return an extra list of RMSD of each band, sum of all delta**2
    '''
    dErr = 0.0
    if ( listBandWeight == None):
        listBandWeight = [1.0 for x in range(0,len(band1[0]))]

    n_k = len(band1)

    list_rmsd = [0 for x in xrange(len(listBandWeight))]

    for i in xrange(n_k-1):
        for j, dBandWeight in enumerate(listBandWeight):
            disp1 = band1[i][j] - band1[i+1][j] 
            disp2 = band2[i][j] - band2[i+1][j]
#           d1 = (band1[i][j] - band1[i+1][j] - band2[i][j] + band2[i+1][j])**2
            d1 = (disp1 - disp2)**2
            if (disp1 * disp2 < 0):#Assign a specal weight for different signs (generally large)
                d1 *= dReverseWeight
            list_rmsd[j] += d1
            dErr += d1 * dBandWeight

    nBand = 0
    for d1 in listBandWeight:
        if ( d1 != 0.0):
            nBand += 1

#   dErr = (dErr / (n_k-1) / nBand)**0.5

#   for i in xrange(len(listBandWeight)):
#       list_rmsd[i] = (list_rmsd[i] / (n_k-1))
#Note: as dispersion by difference is scaled with 1/k, 
#So to make it comparable with RMSD of eigenvalues, it is not necessary to divided by k again
    dErr = (dErr / nBand)**0.5

    if (b_debug):
        return dErr, list_rmsd
    else:
        return dErr


def f_Band_Diff2RMSD(band1,dVBM1,band2,dVBM2,listKPtWeight=None,listBandWeight=None,nBandIndexStart1=0,nBandIndexStart2=0,b_debug=False, b_debug_print=False):
    '''
    Return root-mean-square-deviation of two band structure.
    The comparison is based on relative value to VBM, which should be specify in input parameter.
    The average is the total summation /  number of k-points with weight != 0 / number of bands with weight != 0
    Warning: two band structures must have same number of k-points and bands,unless weight of k-points and bands are specified. 

    :param band1: band data 1
    :param dVBM1: VBM of band 1
    :param band2: band data 2
    :param dVBM2: VBM of band 2
    :param listKPtWeight: Weight of each k-point. If not specified, all KPt will be treated as 1. 1-D array contains weight only, start from first kpt in array.
    :param listWeight:  Weight of each band. If not specified, all band will be treated as 1. 1-D array, start from first band in array.
    :param b_debug: prints more information and return an extra list for all delta**2
    '''
    if ( len(band1) != len(band2) ):
        raise ValueError,"Two band structure must have same number of k-points to compare"

    dErr = 0.0
    if ( listKPtWeight == None):
        listKPtWeight = [1.0 for x in range(0,len(band1))]
    if ( listBandWeight == None):
        listBandWeight = [1.0 for x in range(0,len(band1[0]))]
    
    if (b_debug):
        list_band_err = [[] for x in listBandWeight]

    for i,dKPtWeight in enumerate(listKPtWeight):
        if (len(band1[i]) < len(listBandWeight) or len(band2[i]) < len(listBandWeight)):
            raise ValueError("Insufficient bands at kpt #%i4: band1 %4i %4i weight %4i" % (i, len(band1[i]),len(band2[i]),len(listBandWeight)))
        for j,dBandWeight in enumerate(listBandWeight):
            dErr += (band1[i][j]-dVBM1-band2[i][j]+dVBM2)**2*dKPtWeight*dBandWeight
            if (b_debug):
                list_band_err[j].append((band1[i][j]-dVBM1-band2[i][j]+dVBM2)**2)
            #print(dErr)
    #Count number of points
    nKPt = 0
    for d1 in listKPtWeight:
        if ( d1 != 0.0):
            nKPt += 1
    nBand = 0
    for d1 in listBandWeight:
        if ( d1 != 0.0):
            nBand += 1
    dErr = (dErr / nKPt / nBand)**0.5

    if (b_debug):
#Print RMSD for each band
        if (b_debug_print):
            for i, l1 in enumerate(list_band_err):
                n_k = len(listKPtWeight)
                print("Band %3i : %10.6f w:%10.6f" % (i+1, (sum(l1)/n_k)**0.5, listBandWeight[i]))
        return dErr, list_band_err
    else:
        return dErr

def f_Band_Diff2RMSD_File(stFileName1,stFileName2,nElectron=0,nBandMin=0,nBandMax=0):
    '''
    Compare two file of band structure file in tmc format with given electron count
    :param nElectron: Number of electron to detect VBM to compare. If equals to 0, then do not align.
    :param nBandMin: Lowest index of band to compare, start from 1
    :param nBandMax: Highest index of band to compare, including itself
    '''
    band1 = f_Band_ReadFromFile(stFileName1)
    band2 = f_Band_ReadFromFile(stFileName2)
    (fVBM1,fVBM2) = (0.0,0.0)
    if ( nElectron != 0):
        fVBM1 = f_Band_GetVBMFromElectron(band1,nElectron)
        fVBM2 = f_Band_GetVBMFromElectron(band2,nElectron)
    if ( nBandMin == 0 ):
        nBand = nElectron / 2
    return f_Band_Diff2RMSD(band1,fVBM1,band2,fVBM2,None,[0.0 for x in range(0,nBandMin-1)] + [1.0 for x in range(nBandMin-1,nBandMax)])


def f_Band_ShiftEnergy(list_band,dAlign):
    '''
    Shift a band structure, set reference energy as 0
    :params dAlign: the reference energy of input band structure
    '''
    for kpt in list_band:
        for i in range(0,len(kpt)):
            kpt[i] -= dAlign


def f_InRange(a,min,max):
    '''
    Return if a <= max and a >= min
    '''
    return ( a <= max and a >= min)


def f_AddListAll(list,func):
    '''
    Return summation of [func(x) for x in list]
    '''
    result = 0
    for i in list:
        result += func(i)
    return result

def f_DiffBandDispersion(listBand,i):
    '''
    return leftDiff, right diff and diff of left/right diff of specific position i
    Note: no overbound is checked
    '''
    nBand = min([len(x) for x in listBand])
    list_leftDiff = [0 for j in range(0,nBand)]
    list_rightDiff = [0 for j in range(0,nBand)]
    list_Diff2 = [0 for j in range(0,nBand)]
    for j in range(0,nBand):
        list_leftDiff[j] = listBand[i][j]-listBand[i-1][j]
        list_rightDiff[j] =  listBand[i+1][j]-listBand[i][j]
        list_Diff2[j] = list_rightDiff[j] - list_leftDiff[j]
    return list_leftDiff,list_rightDiff,list_Diff2

def f_Band_DetectCrossing(listBand,listTurnPt = []):
    '''
    Judge band crossing by left derivative and right derivative. Only energy and k-vector information is used.
    If there are different number of bands in different k-point, only lowest N bands will be considered, where N is the number of bands in k-point that contains fewest bands 
    :param listBand: one element represent a series of energy level at specific k-points. These k-points should form a normal polygon as connected in the input order.
    :param listTurnPt: the k-point index in listBand where is a turning point on the polygon, at where band crossing is ignored. Start or ending point
    '''
    #print("==========Start=============")
    debug_print(3,"Turning points:")
    debug_print(3,listTurnPt)

    listBandNew = copy.deepcopy(listBand)
    nBand = min([len(x) for x in listBand])
    dSlopeLimit = 1.5 #x-range of extend line to detect cross

    #Try to get a reference value to use as parallel determination standard
    listDelta = [ abs(listBand[0][x] - listBand[1][x]) for x in range(0,nBand)] #Delta between first and second kpt
    listDelta.sort()
    #Use the change in the middle 
    dDelta = listDelta[len(listDelta)/2]

    dParallelLimitMin = 0.8 *dDelta /10
    dParallelLimitMax = 1.2 *dDelta /10
    dParallelAbsMax = 0.3 *dDelta /10
    i = 1 # the first column is not moved
    while ( i <len(listBand) -1 ): 
        list_leftDiff,list_rightDiff,list_Diff2 = f_DiffBandDispersion(listBandNew,i)
        dStd = f_AddListAll(list_Diff2,abs) #judge standard: sum of abs(diff of left / right)
        #note the previous k-pt needs included if i > 1, for change to band below may affect right diff of previous k-pt 
        if ( i > 1):
            dStd += f_AddListAll(f_DiffBandDispersion(listBandNew,i-1)[2],abs)
        debug_print(2,"k-point #%d Std: %f" % (i+1,dStd))
        #print(list_Diff2)
        #opposite Diff2 indicate a probably cross
        
        #sort the band at this position
        BandOrder = [(x[1],x[0]) for x in enumerate(listBandNew[i])]
        BandOrder.sort()
        BandOrder = [x[1] for x in BandOrder]
        #print("New Order",BandOrder)        
        for j1_2 in range(0,nBand): #search band-pair from buttom to top, in energy order of energy at that k-point of band
            j1 = BandOrder[j1_2]
            if ( j1_2 == nBand-1): # the last band reached, move to next k-pt
                i += 1
                break

            for j2_2 in range(j1_2+1,nBand):
                j2 = BandOrder[j2_2]
                #print(j1,j2)
                # detect if two bands are near enough
                #standard: extend line of lower band at k-pt is higher than higher band on some point
                #          extend line of upper band at k-pt is lower than lower band on some point
                #There are 4 extend line of each pair
                #if detected, shift the band on the right and calculate again
                nChangeStart = -1 #band cross start point; must be 0 ( this point ) or 1 ( next point )
                dHigher1 = listBandNew[i-1][j2] - ( listBandNew[i][j1] - list_rightDiff[j1]*dSlopeLimit )
                dHigher2 = listBandNew[i+1][j2] - ( listBandNew[i][j1] + list_leftDiff[j1]*dSlopeLimit )
                dLower1 =  -listBandNew[i-1][j1] + ( listBandNew[i][j2] - list_rightDiff[j2]*dSlopeLimit )
                dLower2 =  -listBandNew[i+1][j1] + ( listBandNew[i][j2] + list_leftDiff[j2]*dSlopeLimit )
                if ( dHigher1 < 0 ):
                    if ( dHigher2 < dHigher1):
                        nChangeStart = 1
                    else:
                        nChangeStart = 0
                elif ( dHigher2 < 0): 
                    nChangeStart = 1
                #the other two seems redundant ?
                #whether is not of 
                else:
                    break

                if ( nChangeStart == -1):
                    raise ValueError,"Unexpected error in band crossing calculation"

                debug_print(2,"  %f %f %f / %f %f %f"%(listBandNew[i-1][j1],listBandNew[i][j1],listBandNew[i+1][j1],listBandNew[i-1][j2],listBandNew[i][j2],listBandNew[i+1][j2]))
               
                if  ( (nChangeStart+i) in listTurnPt or  (nChangeStart + i - 1) in listTurnPt  ) : #Ignore at turning point or next value of turning point ( move which will break connecting at turning point )
                    debug_print(2,"  Ignore cross at turning point, band %i / %i" % ( j1+1,j2+1))
                    #i += 1
                    continue

                
                #nearly parallel band may be included but does not require cross, exclude them
#There are two standard of parallel from derivate
#not used now
                #print(list_leftDiff[j1],list_leftDiff[j2],list_rightDiff[j1],list_rightDiff[j2])

                #if (abs(list_leftDiff[j1]-list_leftDiff[j2]) < dParallelAbsMax and abs(list_rightDiff[j1]-list_rightDiff[j2]) < dParallelAbsMax):
                #    debug_print(3,"\tAbsolute Parallel detected")
                #   break                

                #if (f_InRange(list_leftDiff[j1]/list_leftDiff[j2],dParallelLimitMin,dParallelLimitMax) and f_InRange(list_rightDiff[j1]/list_rightDiff[j2],dParallelLimitMin,dParallelLimitMax) ):
                #    debug_print(3,"\tRelative Parallel detected")
                #    break                     
                
                #make the change and see whether it is better, if not, refuse it                
                debug_print(2,"  Seems Cross: #%d k-point, band %d / %d, change %d" % (i+1,j1+1,j2+1,nChangeStart))

                
                # move down one by one and shift the lowest to highest
                # or move up one by one and shift the highest to lowest
                # both are needed to test due to band between cross may change the result
                # and choose the lower std one
                stReverse = "lowest to highest"
                listBandNew2 = copy.deepcopy(listBandNew)
                for i2 in range(i+nChangeStart,len(listBand) ): 
                    for j in range(j1_2,j2_2):
                        listBandNew2[i2][BandOrder[j]] = listBandNew[i2][BandOrder[j+1]]
                    listBandNew2[i2][BandOrder[j2_2]] = listBandNew[i2][BandOrder[j1_2]]
                #calculate standard again
                list_leftDiff_new2,list_rightDiff_new2,list_Diff_new2 = f_DiffBandDispersion(listBandNew2,i)
                dStd2 = f_AddListAll(list_Diff_new2,abs)
                if ( i > 1):
                    dStd2 += f_AddListAll(f_DiffBandDispersion(listBandNew2,i-1)[2],abs)                   
                
                listBandNew3 = copy.deepcopy(listBandNew)
                for i2 in range(i+nChangeStart,len(listBand)):
                    for j in range(j1_2+1,j2_2+1):
                        listBandNew3[i2][BandOrder[j]] = listBandNew[i2][BandOrder[j-1]]
                    listBandNew3[i2][BandOrder[j1_2]] = listBandNew[i2][BandOrder[j2_2]]
                list_leftDiff_new3,list_rightDiff_new3,list_Diff_new3 = f_DiffBandDispersion(listBandNew3,i)
                dStd3 = f_AddListAll(list_Diff_new3,abs)
                if ( i > 1):
                    dStd3 += f_AddListAll(f_DiffBandDispersion(listBandNew3,i-1)[2],abs)                
                
                if ( dStd3 < dStd2):
                    stReverse = "highest to lowest"
                    list_leftDiff_new2,list_rightDiff_new2,list_Diff_new2 = list_leftDiff_new3,list_rightDiff_new3,list_Diff_new3
                    dStd2 = dStd3
                    listBandNew2 = listBandNew3
                
                debug_print(2,"    New Std: %f / Old %f" % (dStd2,dStd))
                #if ( i+1 == 5 and j1+1 == 5 and j2+1 == 3):
                #    pass
                #    return listBandNew2
                
                #Judge whether this change makes result better in distort
                #if not, do not use it
                if ( dStd2 < dStd - dParallelAbsMax / 2):
                    debug_print(1,"      Cross: #%d k-point, band %d / %d, %s" % (i+1,j1+1,j2+1,stReverse))
                    listBandNew = listBandNew2
                    dStd = dStd2
                    list_leftDiff,list_rightDiff,list_Diff2 = list_leftDiff_new2,list_rightDiff_new2,list_Diff_new2
                else:
                    debug_print(1,"      Not a cross")

    #Sort band with lowest enengy in the band

    listBandNew2 = [ [x[i] for x in listBandNew] for i in range(0,nBand)]
    listBandNew2.sort(lambda x,y : cmp(min(x),min(y)))
    listBandNew3 =  [ [x[i] for x in listBandNew2 ] for i in range(0,len(listBandNew))]

    return listBandNew3


def f_band_resolve_crossing(list_band_org,func_comp,func_eig=None,list_turnkpt=None, list_aux=[]):
    '''
    Resovle band crossing by testing derivatives of specific parameters
    We assumed the input band structure is calculated on linear grid except some individual k-points, where the reconnection is ignored.
    :param list_band_org: [k-point, band index] 2D array of list that represents an eigenvalue.
    :param func_comp: a function that accept two list and return the difference
    :param func_eig: return the eigenvalue in unit eV, which is used to restrict search range to speed up
    :param list_turnkpt: list of k-points index where k-grid is a turning points, starts from 0
    :param list_aux: a list of arrays that related to band structure one by one and will be reversed, each array in [k-point, band_index]

    :return: new band structure and aux(optional)
    '''
    eig_diffmax = 1.0 #Only search in range of 1.0eV between near k-points

    def diff_lr_single(x0, x1, x2):
        a = func_comp(x1, x0)
        b = func_comp(x2,x1)
        return (a, b, b-a)

    def diff_lr(list1,i):
        list2 = []
        for x in list1:
            list2.append(diff_lr_single(x[i-1], x[i], x[i+1]))
        return list2

#Note all array we use here shuold be [band index, k-point]
#So we transpose all inputs

#Copy all pointer of each eigenvalue list
    list_band = list(map(list, zip(*list_band_org)))

    for i in range(len(list_aux)):
        list_aux[i] = list(map(list, zip(*(list_aux[i]))))
    
    num_band = len(list_band)
    num_kpt = len(list_band[0])

#Init values 
    if (list_turnkpt is None):
        list_turnkpt = [] 

#Process k-points one by one
    for ix_kpt in xrange(1,num_kpt-2):
#       print("Kpt %i / %i" % (ix_kpt, num_kpt - 2))
        if (ix_kpt in list_turnkpt or ix_kpt+1 in list_turnkpt):
            continue
#if left are fixed, every move of each point affect this and right kpt
#Compute total error of this point and  the right one 
        diff = diff_lr(list_band,ix_kpt)

#Create an order to enumerate pair, which must be unrelated to current condition
#Use a list to track them ( the ix_kpt + 1 column )
        list_current_index = range(num_band)
        
#Enumerate all pair swap between i and i+1 kpt
#i1 and i2 are index of band before any exchange, ix_band1 and ix_band2 are that currently of ix_kpt+1 column
        i1 = 0
        b_restart = False
        last_exchange = None
        while (i1 < num_band-1):
#           print("Loop1 %i / %i" % (i1, num_band))
#           ix_band1 = list_current_index[i1]
            i2 = i1 + 1
            while (i2 < num_band):
#               print("Loop2 %i / %i" % (i2, num_band))
#               print("%i %i %i"%(ix_kpt,i1,i2))
#If func_eig is provided then we test only in range of 
                if (func_eig is not None):
                    if (abs(func_eig(list_band[i1][ix_kpt])-func_eig(list_band[i2][ix_kpt])) > eig_diffmax):
                        i2 += 1
                        continue

#               ix_band2 = list_current_index[i2]

                diff_old1 = func_comp(list_band[i1][ix_kpt+1],list_band[i1][ix_kpt])
                diff_old2 = func_comp(list_band[i2][ix_kpt+1],list_band[i2][ix_kpt])
                diff_old = abs(diff[i1][0]-diff_old1)\
                        +abs(diff[i2][0]-diff_old2)\
                        +abs(list_band[i1][ix_kpt+2]-list_band[i1][ix_kpt+1]-diff_old1)\
                        +abs(list_band[i2][ix_kpt+2]-list_band[i2][ix_kpt+1]-diff_old2)

                diff_new1 = func_comp(list_band[i2][ix_kpt+1],list_band[i1][ix_kpt])
                diff_new2 = func_comp(list_band[i1][ix_kpt+1],list_band[i2][ix_kpt])
                diff_new = abs(diff[i1][0]-diff_new1)\
                        +abs(diff[i2][0]-diff_new2)\
                        +abs(list_band[i1][ix_kpt+2]-list_band[i1][ix_kpt+1]-diff_new2)\
                        +abs(list_band[i2][ix_kpt+2]-list_band[i2][ix_kpt+1]-diff_new1)

                if (diff_new < diff_old - 1e-6):#Smaller, can be reconnected. with a tolerance
#                   print("Exchange: kpt %i: %i %i with %8.4f > %8.4f " % (ix_kpt,i1,i2,diff_old,diff_new))
#                   print("  Old lines1 : (%.3f) %.3f / %.3f / %.3f / %.3f (%.3f)" % tuple(
#                       [diff[i1][0]] + list_band[i1][ix_kpt-1:ix_kpt+3] + [list_band[i1][ix_kpt+2]-list_band[i1][ix_kpt+1]]))
#                   print("  Old lines2 : (%.3f) %.3f / %.3f / %.3f / %.3f (%.3f)" % tuple(
#                       [diff[i2][0]] + list_band[i2][ix_kpt-1:ix_kpt+3] + [list_band[i2][ix_kpt+2]-list_band[i2][ix_kpt+1]]))
                    if (last_exchange is not None and 
                            (last_exchange[0] == ix_kpt and last_exchange[1] == i1 and last_exchange[2] == i2)):
                            print(list_band[ix_band1][ix_kpt], list_band[ix_band2][ix_kpt])
#                           raise ValueError("Cannot resolve the band strcture, dead loop encountered")
                    last_exchange = [ix_kpt, i1, i2]

                    list_tmp = list_band[i1][ix_kpt+1:]
                    list_band[i1][ix_kpt+1:] = list_band[i2][ix_kpt+1:]
                    list_band[i2][ix_kpt+1:] = list_tmp
#Update diff array

#Also change the list_aux
                    for aux1 in list_aux:
                        list_tmp = aux1[i1][ix_kpt+1:]
                        aux1[i1][ix_kpt+1:] = aux1[i2][ix_kpt+1:]
                        aux1[i2][ix_kpt+1:] = list_tmp
                        
#                   tmp = list_current_index[i1]
#                   list_current_index[i1] = list_current_index[i2]
#                   list_current_index[i2] = tmp
#Rescan until converged
                    b_restart = True
                    break
                
                i2 += 1
            if (b_restart):
                b_restart = False
                i1 = 0
                continue
            i1 += 1

#   for eigs in list_band:
#       print(eigs)

#Transpose all back
    list_band = list(zip(*list_band))

    for i in range(len(list_aux)):
        list_aux[i] = list(zip(*(list_aux[i])))

#Return
    return list_band, list_aux



def f_Band_CheckMax(nIndex , nMaxIndex, stName,bPrint=True):
    '''
    Check if nIndex is larger than Max index and show error message. Return True if check passed, vice versa.
    '''
    bPass = True
    if ( nIndex < 0 ):
        if ( bPrint):
            print("ERROR: "+stName+" < 1  !!!")
            print("  --- Check the Fermi energy")
        bPass = False
    if ( nIndex > nMaxIndex ):
        if ( bPrint):
            print("ERROR: "+ stName +" is larger than the number of available bands  !!!")
            print("  --- Check the Fermi energy")
        bPass = False
    
    return bPass


def f_Band_GetVBMFromElectron(listBand,nElectron,spin=1,b_print=False):
    '''
    Get VBM of a semiconductor by electron number
    spin-polarized version is not implemented because it is enough to do this:
    * Combine two band structures of different spins as one
    * pass that band structure with spin=1

    :param listBand: the eigenvalues
    :param nElectron: number of electrons
    :param spin: 1 for unpolarized and 2 for spin-polarized
    :param b_print: display prompts
    '''
    if (spin != 1):
        raise ValueError("Cannot use spin-polarized one here")
    if ( nElectron <= 0) :
        raise ValueError,"Number of electrons must > 0, currently %i" % nElectron
    if ( nElectron / 2 * 2 != nElectron and spin == 1):
        if (b_print):
            print("Number of electrons must be even for an insulator, currently %i" % nElectron)
#       return f_Band_GetFermiFromElectron(listBand,nElectron)
        return None

    #print(nElectron)

    nBand = nElectron / 2
    if ( min([len(x) for x in listBand]) < nBand ):
        raise ValueError,"Band structure has fewer bands (%i) than electron number / 2 (%i)" % ( min([len(x) for x in listBand]) ,nBand)

    vbm = max([ x[nBand-1] for x in listBand])
    cbm = min([ x[nBand] for x in listBand])
    if (vbm <= cbm):
        return vbm
    else:#It is a metal!
        return None

def f_Band_GetFermiFromElectron(list_eig,list_wk,num_electron,num_spin=1):
    '''
    Get Fermi of a metal by electron number. This is only an estimation.
    Weights of k-point is required
    Warning: this is not appropriate because band structure k-points is not suitable for determine Fermi energy 

    spin-polarized version is not implemented, refer to :func:f_Band_GetVBMFromElectron
    '''
    if (num_spin == 2):
        raise NotImplementedError("Spin-polarized cases not supported yet")
#Normalized
    wk1 = [x*1.0/sum(list_wk) for x in list_wk]
#Resort data
    eig1 = list(itertools.chain( *[[ [y,wk1[i]] for y in x] for i,x in enumerate(list_eig)]))
    eig1.sort(key=lambda x:x[0])
    n = 0.0
    for i,eigw in enumerate(eig1):
        n += eigw[1]*2
        if ( n > num_electron):
            break
#Just take the average
    return (eig1[i][0]+eig1[i][0])/2

                               
def f_Band_GetGap(listBandOrg,fFermiEnergy=None,bPrint=True):
    '''
    Get gap from band structure and Fermi energy
    :param fFermiEnergy: the Fermi energy or VBM
    :param bPrint: Print messages, default True
    :return : band gap , detected new VBM ,detected HO band index ( start from 0).    If the system is metal, then return None, None, HO band index
    '''
    resNone = (None,None,None)
    #make a copy, we may shift energy later
    listBand = copy.deepcopy(listBandOrg) 

    if ( fFermiEnergy == None):
        if ( bPrint):
            print("Must specify Fermi energy to determine gap.")
        return resNone

    nBandCountMin = min([len(x) for x in listBand])
    fEnergyMin = min([min(x) for x in listBand])
    fEnergyMax = max([max(x) for x in listBand])
    
    if ( bPrint):
        print("----------------------------------------------------")
        print("         Band Gap       Analysis                    ")
        print("----------------------------------------------------")
        print('Range of bands considered: 1 ' + str(nBandCountMin))
    
    if ( fEnergyMin > fFermiEnergy or fEnergyMax < fFermiEnergy ):
        if ( bPrint):
            print("WARNING from bandanaly: ")
            print(" - Fermi energy outside the energy range of bande!")
            print("  Fermi Energy:   " + str(fFermiEnergy))
            print(" - Maximal energy:    "+str(fEnergyMax))
            print(" - Minimal energy:    "+str(fEnergyMin))
        return resNone
    
    # find HO and LU . Index is the position in array
    nHOIndex = -1
    nLUIndex = 1000

    #Sometimes the Fermi energy is a little lower than the HO, which cause band gap detection error ( the LU index is lowered )
    #So check how many points 
    dErr = 0.01
    nHOIndex2 = -1
    nLUIndex2 = 1000
    
    # note: in Fortran array starts from 0 but in Python from 1
    for listEnergy in listBand:
        #Sk mean single k point
        #nSkHOIndex = -1
        #nSkLUIndex = 1000
        for i in range(0,nBandCountMin):
            #exact standard
            if ( listEnergy[i] <= fFermiEnergy ):
                nHOIndex = i if nHOIndex < i else nHOIndex
            else:
                nLUIndex = i if nLUIndex > i else nLUIndex
            #loose standard
            #Attention: fFermiEnergy may be NaN or Inf, so do not manipulate it !
            if ( listEnergy[i]+dErr <= fFermiEnergy ):
                nHOIndex2 = i if nHOIndex2 < i else nHOIndex2
            if ( listEnergy[i]-dErr >= fFermiEnergy):
                nLUIndex2 = i if nLUIndex2 > i else nLUIndex2
   
    if  not ( f_Band_CheckMax(nHOIndex,nBandCountMin,"HO Band",False) and f_Band_CheckMax(nLUIndex,nBandCountMin,"LU Band",False) ):
        return resNone
    if  not ( f_Band_CheckMax(nHOIndex2,nBandCountMin,"HO Band",False) and f_Band_CheckMax(nLUIndex2,nBandCountMin,"LU Band",False) ):
        return resNone

    #print(nHOIndex,nLUIndex)
    #print(nHOIndex2,nLUIndex2)

#If results from two standard are different, then check
    if ( nHOIndex != nHOIndex2 or nLUIndex != nLUIndex2):
        #If there is a gap between loose standard results, then it should be right; else use previous one.
        if ( min([x[nLUIndex2] for x in listBand]) - max([x[nHOIndex2] for x in listBand]) > dErr):
            if ( bPrint):
                print("WARNING: It seems that specified VBM deviates from the detected VBM, use the detected one")
            nHOIndex = nHOIndex2
            nLUIndex = nLUIndex2
    
    # error checking

    if ( nHOIndex >= nLUIndex):
        if ( bPrint):
            print("Metal system detected")
            print('Highest occupied band: '+str(nHOIndex+1))
            print('Lowest unoccupied band:'+str(nLUIndex+1))
        return 0.0,fFermiEnergy,nHOIndex

#Below part is try to parse band as gap exist, skipped
    if ( nHOIndex > nLUIndex ):
        if ( bPrint):
            print("WARNING: Valence and Conductance bands overlap !!")
            print("Swap HO Band and LU Band...")
        nLUIndex,nHOIndex = nHOIndex,nLUIndex
    else:
# Sometimes when the Fermi energy is not calculated very accurately, judging
# whether it is conductance or valence band only according to its sign is not reliable
# test which dominate in k space and change the other
        if ( nHOIndex == nLUIndex ):
            if ( bPrint):
                print("WARNING: Valence and Conductance bands identical !!")
            nCount = 0
            for aBand in listBand:
                if ( aBand[nHOIndex] > fFermiEnergy ):
                    nCount += 1
            if (nCount >= len(listBand) /2):
                nHOIndex -= 1
            else:
                nLUIndex += 1
    if ( bPrint):
        print('Highest occupied band: '+str(nHOIndex+1))
        print('Lowest unoccupied band:'+str(nLUIndex+1))
    
    # from HO and LU , find energy and gap
    f_Band_ShiftEnergy(listBand,fFermiEnergy)
    
    
    list_HOBand= [ x[nHOIndex] for x in listBand ]
    list_LUBand= [ x[nLUIndex] for x in listBand ]
    fHOEnergy = max(list_HOBand)
    nHOkPointIndex = list_HOBand.index(fHOEnergy)
    fLUEnergy = min(list_LUBand)
    nLUkPointIndex = list_LUBand.index(fLUEnergy)
    
    dGap = fLUEnergy - fHOEnergy

    if ( bPrint):
        print(':BandGap =  %.4f eV' % dGap)

        if ( nHOkPointIndex <> nLUkPointIndex):
            print("Indirect gap:")
            # also print direct gap at HO and LU k point
            # the coordinate of k-vector is not output now
            print("\t%.4f eV at VBM k= %s, ik= %d"% (listBand[nHOkPointIndex][nLUIndex]-fHOEnergy, "skip" ,nHOkPointIndex+1))
            print("\t%.4f eV at CBM k= %s, ik= %d"% (-listBand[nLUkPointIndex][nHOIndex]+fLUEnergy,"skip",nLUkPointIndex+1))
        else:
            print('Direct gap at k=' + "skip" + ' ik=%d' % nLUkPointIndex)
    
    
    #print(":Range of each band (with respect to VBM):")
    #print("    n     Bottom       Top         Width")
    #for i in range(0,nBandCountMin):
    #    aBand =  [x[i] for x in listBand]
    #    fEnergyMax = max(aBand)-fHOEnergy
    #    fEnergyMin = min(aBand)-fHOEnergy
    #    print("%5i%12.4f%12.4f%12.4f" % ( i+1, fEnergyMin, fEnergyMax, fEnergyMax-fEnergyMin) )

    return dGap,fHOEnergy+fFermiEnergy,nHOIndex
           
def f_Band_Analyse(listBand,fRefEnergy=None):
    '''
    Print band structure information summary
    :param fRefEnergy: the reference level which used as 0 ( normally VBM )
    '''

    fHOEnergy = 0
    if ( fRefEnergy is not None):
        fHOEnergy = fRefEnergy
    #print(fHOEnergy)

    nBandCountMin = min([len(x) for x in listBand])
    fEnergyMin = min([min(x) for x in listBand])
    fEnergyMax = max([max(x) for x in listBand])
    
    print("----------------------------------------------------")
    print("         Band Structure Analysis                    ")
    print("----------------------------------------------------")
    print('Range of bands considered: 1 ' + str(nBandCountMin))
    
    print(":Range of each band:")
    print("    n     Bottom       Top         Width")
    for i in range(0,nBandCountMin):
        aBand =  [x[i] for x in listBand]
        fEnergyMax = max(aBand)-fHOEnergy
        fEnergyMin = min(aBand)-fHOEnergy
        print("%5i%12.4f%12.4f%12.4f" % ( i+1, fEnergyMin, fEnergyMax, fEnergyMax-fEnergyMin) )

def f_Band_CombineSOC(listBand):
    '''
    For some system, soc effect is included in calculation, make every orbtial contains only 1 electron in calculation.
    However, some other package only deal with 2-electron-per-orbtial case.
    This function use average value of each 2 band in the order of energy as to create a non-soc band structure
    '''
    listBand2 = []
    for aBand in listBand:
        aBand2 = []
        for i in range(0,len(aBand)/2):
            aBand2.append((aBand[i*2]+aBand[i*2+1])/2)
        listBand2.append(aBand2)

    return listBand2

def f_Band_EstimateNumberOfConductionBand(aCell):
    '''
    Estimate the number of condunction bands one may interest in one system just by atom names and counts
    '''
    nVal = 0
    for Atom in aCell.listAtom:
        nVal += f_Element_Z_to_ExtraValenceElectronCount(f_Element_Symbol_to_Z(Atom[0]))

    return nVal/2

def f_Band_Plot(list_bs,list_shift,ix_kpt=0,y_fermi=None,y_min=None,y_max=None,dirname_out=None,name_prefix="plotband", plot_setting=None):
    '''
    Plot one or more band structure on the same picture to xmgrace.

    :param list_band: A list contains BandsT objects
    :param list_shift: A list contains energy shifts for corresponding band structures. If set to None, no shift will be made. Note this list stored reference energy, which means all new eigenvalues will be (old-shift)
    :param y_min: minimum value of y axis
    :param y_max: maximum value of y axis
    :param name_prefix: all output will be saved as name_prefix.xxx
    :param ix_kpt: for simplicity, only one k-points set will be plotted on the picture, and this parameter specify which one. Default is the first.
    '''
    if (dirname_out is None):
        filename_full = name_prefix
    else:
        if (not os.path.exists(dirname_out)):
            os.mkdir(dirname_out)
        filename_full = os.path.join(dirname_out,name_prefix)

    lines = f_band_plot_xmgrace(list_bs, list_shift, ix_kpt, y_fermi, y_min, y_max, plot_setting)

    with open("%s.agr" % filename_full, 'w') as fXM:
        fXM.writelines(lines)
        
    try:
        subprocess.check_call(["xmgrace", "-hardcopy", "-printfile", "%s_xm.ps" % filename_full, "%s.agr" % filename_full])
        print("Please see result in %s.agr & %s_xm.ps ( xmgrace )" % (name_prefix,name_prefix))
    except OSError:
        print("Cannot invoke xmgrace to generate pictures")
        print("Please see result in %s.agr" % (name_prefix))
    except subprocess.CalledProcessError:
        print("Cannot invoke xmgrace to generate pictures")
        print("Please see result in %s.agr" % (name_prefix))

    return


def f_band_plot_xmgrace(list_bs,list_shift,ix_kpt=0,y_fermi=None,y_min=None,y_max=None, plot_setting=None):
    '''
    Plot one or more band structure on the same picture to xmgrace.

    :param list_band: A list contains BandsT objects
    :param list_shift: A list contains energy shifts for corresponding band structures. If set to None, no shift will be made. Note this list stored reference energy, which means all new eigenvalues will be (old-shift)
    :param y_min: minimum value of y axis
    :param y_max: maximum value of y axis
    :param ix_kpt: for simplicity, only one k-points set will be plotted on the picture, and this parameter specify which one. Default is the first.
    :param plot_setting: stores xmgrace parameters
    '''

#Check band structure
    for band in list_bs:
        if ( len(band.list_eig[0]) <= 1 ):
            raise ValueError,"Band structure contains too few k-points"

#Check shift
    if (list_shift is None):
        list_shift = [0.0 for x in list_bs]
#Set None to 0
    for i in xrange(len(list_shift)):
        if (list_shift[i] is None):
            list_shift[i] = 0.0

#Set x range
    list_dist_all = []
    #convert kpt to pure kpt list and name list
    list_name = [] #specail k-point name
    list_kp = [] #kpoint coordinate

#Split a bandstructure with break k-points as multiple band structures
#Extract eigenvalues from BandsT
#Also copy shift
    list_bs2 = []
    list_shift2 = []
    list_bands_ix_org = [] #Mark the index of unsplitted bandstructure that this bandstructure belongs to

    for u,band in enumerate(list_bs):
        list_dist, list_name1, list_break = band.kpt.get_k_xaxis()
        if (u == ix_kpt):#Only set name for specific set
            list_name = list_name1
        if (len(list_break) > 1):#Split band structure
            list_break.append(len(list_dist))
            for j in range(len(list_break)-1):
                j1 = list_break[j]
                j2 = list_break[j+1]

                list_bs2.append(band.list_eig[j1:j2])
                list_dist_all.append(list_dist[j1:j2])
                list_shift2.append(list_shift[u])
                list_bands_ix_org.append(u)
        else:
            list_bs2.append(band.list_eig)
            list_dist_all.append(list_dist)
            list_shift2.append(list_shift[u])
            list_bands_ix_org.append(u)




#Find y range for the total picture
#   list_emin = []
#   list_emax = []
#   for i,band in enumerate(list_bs):
#       list_emin.append(min([min(x) for x in band.list_eig])-list_shift[i])
#       list_emax.append(max([max(x) for x in band.list_eig])-list_shift[i])
#   Emin = min(list_emin)
#   Emax = max(list_emax)
    Emin = min(map(lambda x,y:x.eig_min-y,list_bs,list_shift))
    Emax = max(map(lambda x,y:x.eig_max-y,list_bs,list_shift))

    xmin = 0
#   xmax = len(list_band)-1
    ymin = Emin-0.5
    ymax = Emax+0.5
    if (y_min != None):
        ymin = y_min
    if (y_max != None):
        ymax = y_max



    #reset x axis
    xmax = max([list_dist[-1] for list_dist in list_dist_all])

    #create xmgrace file
    lines = []
    lines.append('''# Grace project file
#
@version 50000
''')
    
    #Fermi level
    if ( y_fermi != None):
        lines.append('''@with line
@    line on
@    line g0    
@    line linestyle 4
@    line loctype world
@    line %f,%f,%f,%f
@line def
''' % (xmin,y_fermi,xmax,y_fermi) )


#manual axis;prefered method 
    lines.append('''@with g0
@    world %f,%f,%f,%f
@    title "Bandstructure"
@    autoticks
@    xaxis  label char size 1.0
@    xaxis  ticklabel char size 1.0
@    xaxis  tick major grid on
@    xaxis  tick spec type both
@    yaxis  label "Energy(eV)"
@    yaxis  tick off 
@    legend off
@    xaxis  tick spec %d
''' % (xmin,ymin,xmax,ymax,len(list_name)))

    for i,kpt in enumerate(list_name):
        name = kpt[1]
        if ( len(name.replace("'",""))>1):#greek character
            name = "\\x%s" % name[0]
        lines.append('''@    xaxis  tick major %f,%f
@    xaxis  ticklabel  %d, "%s"
''' % (i,kpt[2],i,name) )

#Default style
    list_style = [ \
            [1,1,0,0.0,0],\
            [0,0,1,0.5,2],\
            [0,0,2,0.5,4]\
            ]

    ix = -1

    for i,band in enumerate(list_bs2):
        for j in range(0,len(band[0])):
            ix += 1
            ix_style = list_bands_ix_org[i]
            if (ix_style >= len(list_style)):
                ix_style = ix_style%len(list_style)
            lines.append("@    s%d line type %i\n" % (ix,list_style[ix_style][0]))
            lines.append("@    s%d line color %i\n" % (ix,list_style[ix_style][1]))
            lines.append("@    s%d symbol %i\n" % (ix,list_style[ix_style][2]))
            lines.append("@    s%d symbol size %f\n" % (ix,list_style[ix_style][3]))
            lines.append("@    s%d symbol color %i\n" % (ix,list_style[ix_style][4]))

    ix = -1

    for u,band in enumerate(list_bs2):
        for i in range(0,len(band[0])):
            ix += 1
            lines.append("@target G0.S%d\n" % ix)
            lines.append("@type xy\n")
            for j in range(0,len(band)):
                lines.append("%14.7f %14.7f\n" %  (list_dist_all[u][j],band[j][i]-list_shift2[u]))
            lines.append("&\n")    
   
#lines may contain multi-line string, split them
    lines2 = split_multiline_to_single(lines)

#Adjust the plotting
    if (plot_setting is not None):
        lines2 = set_property_auto(lines2, dic_property=plot_setting)

    return lines2


def f_PlotBand(aKPt,list_band,eig_fermi = 0,stPrefix="plotband",bAlignData=False,list_band2=None,eig_fermi2 = 0,y_min=None,y_max=None):
    '''
    Plot band from k-points list ( name and k-vectors ) and band energy in each k-point.
    @todo gnuplot part does not support second set
    :param aKPT: KPointsT type object represent k-points
    :param list_band: The energy-index list
    :param eig_fermi: The fermi energy of band structure data
    :param bAlignData: if true, move data globally to make E_F = 0
    :param list_band2: if set to any value, then plot second band on the same picture
    :param eig_fermi2: shift the energy of band 2
    :param y_min: minimum value of y axis
    :param y_max: maximum value of y axis
    '''

    #raise error when there is only 0 or 1 kpt
    if ( len(list_band[0]) <= 1 ):
        raise ValueError,"Band structure contains too few k-points"

    if ( bAlignData ): #fermi level for output
        eig_fermi1 = 0
        eig_shift = eig_fermi
        eig_shift2 = eig_fermi2
    else:
        eig_fermi1 = eig_fermi
        eig_shift = 0
        eig_shift2 = 0

    Emin = min([min(x) for x in list_band])
    Emax = max([max(x) for x in list_band])

    if ( bAlignData ):
        Emin = Emin - eig_shift
        Emax = Emax - eig_shift

    xmin = 0
    xmax = len(list_band)-1
    ymin = Emin-0.5
    ymax = Emax+0.5
    if (y_min != None):
        ymin = y_min
    if (y_max != None):
        ymax = y_max

    #list_name: The special k-points list as [index,name]
    #list_kp: The cooridnate of k-points ( the unit is not considered here )

    #convert kpt to pure kpt list and name list
    list_name = [] #specail k-point name
    list_kp = [] #kpoint coordinate
    list_dist = [] #x-axis in plot rely on k-k distance

    for i,kpt in enumerate(aKPt.listKPt):
        list_kp.append(kpt[1:4])
        if ( i == 0):
            list_dist.append(0.0)
        else:
            list_dist.append(list_dist[-1] + lu.f_List_norm( lu.f_List_Op_List(kpt[1:4],"-",aKPt.listKPt[i-1][1:4])))
        if ( kpt[0] != ""):
            name = kpt[0]
            if ( name[0] == '\\' ): #delete first backslash as neither gnuplot nor xmgrace use it. SIESTA use "\" for Latex. 
                name = name[1:]
            list_name.append([i,name,list_dist[-1]])

    #reset x axis
    xmax = list_dist[-1]

    #write to data file
    nXCol = 4
    nDataColStart = 4
    fOut = open(stPrefix + ".dat",'w')
    for i in range(0,len(list_band)):
        fOut.write("%f\t%f\t%f" % tuple(list_kp[i]))
        fOut.write("\t%f" % list_dist[i])
        for band in list_band[i]:
            fOut.write("\t%f" % (band-eig_shift ))
        fOut.write("\n")
    fOut.close()
    #create gnuplot file
    fOut = open(stPrefix + ".gp",'w')
    fOut.write('''#!/usr/bin/env gnuplot
    
set style line 1 lt 1 lw 2 pt 3 ps 0.5 lc rgb "black"

set terminal postscript color enh
set output "%s_gp.ps"
datafile="%s.dat"

xmax = %f
ymin = %f
ymax = %f

set title "Band Structure" font "Times-Roman"
set ylabel "Energy(eV)" font "Times-Roman"
set xrange [0:xmax]
set yrange [ymin:ymax]
set nokey
unset xtics
set ytics nomirror font "Times-Roman"
set border 15

''' %(stPrefix,stPrefix,xmax,ymin,ymax))
    
    #add Fermi level
    if ( eig_fermi1 != None):
        fOut.write("fermi = %f\n" % eig_fermi1)
        fOut.write('set arrow from 0,fermi to xmax,fermi nohead lt 1 lw 2 lc rgb "black"\n')
        fOut.write('set label "E_F" at -xmax/40.0,fermi right offset 0,0 font "Times-Roman" \n')
    
    #add special point
    for kpt in list_name:
        fOut.write('set arrow from %f,ymin to %f,ymax nohead lt 1 lw 2 lc rgb "black"\n' % (kpt[2],kpt[2]))
        name = kpt[1]
        if ( len(name)>1):#greek character
            name = "{/Symbol %s}" % name[0]
        fOut.write('set label "%s" at %f,ymin-(ymax-ymin)/40.0 center font "Times-Roman" \n' % (name,kpt[2]))
    fOut.write('\n')
    
    #plot data
    fOut.write("plot ")
    for i in range(0,len(list_band[0])-1):
        fOut.write("datafile using %d:%d w l ls 1," % (nXCol,i+1+nDataColStart))
    fOut.write("datafile using %d:%d w l ls 1\n\n" % (nXCol,len(list_band[0])+nDataColStart))
        
    fOut.close()
    #plot it 
    os.chmod("%s.gp" % stPrefix, stat.S_IRWXO | stat.S_IRWXG | stat.S_IRWXU )
    commands.getoutput("%s.gp" % stPrefix)
    
    #create xmgrace file
    fXM = open("%s.agr" % stPrefix,'w')
    fXM.write('''# Grace project file
#
@version 50000
''')
    
    #Fermi level
    if ( eig_fermi1 != None):
        fXM.write('''@with line
@    line on
@    line g0    
@    line linestyle 4
@    line loctype world
@    line %f,%f,%f,%f
@line def
''' % (xmin,eig_fermi1,xmax,eig_fermi1) )


#manual axis;prefered method 
    fXM.write('''@with g0
@    world %f,%f,%f,%f
@    autoticks
@    xaxis  label char size 1.5
@    xaxis  ticklabel char size 1.5
@    xaxis  tick major grid on
@    xaxis  tick spec type both
@    yaxis  label "Energy(eV)"
@    yaxis  label char size 1.5
@    yaxis  tick off 
@    yaxis  ticklabel char size 1.5
@    legend off
@    xaxis  tick spec %d
''' % (xmin,ymin,xmax,ymax,len(list_name)))

    for i,kpt in enumerate(list_name):
        name = kpt[1]
        if ( len(name)>1):#greek character
            name = "\\x%s" % name[0]
        fXM.write('''@    xaxis  tick major %f,%f
@    xaxis  ticklabel  %d, "%s"
''' % (i,kpt[2],i,name) )

   
    for i in range(0,len(list_band[0])):
        fXM.write("@    s%d line type 0\n" % i)
        fXM.write("@    s%d symbol 1\n" % i)
        fXM.write("@    s%d symbol size 0.5\n" % i)
        fXM.write("@    s%d symbol color 1\n" % i)

    if (list_band2 != None):
        for i in range(len(list_band[0])+0,len(list_band[0])+len(list_band2[0])):
            fXM.write("@    s%d line color 2\n" % i)
            fXM.write("@    s%d line linewidth 2\n" % i)
   
    for i in range(0,len(list_band[0])):
        fXM.write("@target G0.S%d\n" % i)
        fXM.write("@type xy\n")
        for j in range(0,len(list_band)):
            fXM.write("%14.7f %14.7f\n" %  (list_dist[j],list_band[j][i]-eig_shift ) )
        fXM.write("&\n")    
    if (list_band2 != None):
        for i in range(0,len(list_band2[0])):
            fXM.write("@target G0.S%d\n" % (len(list_band[0])+i))
            fXM.write("@type xy\n")
            for j in range(0,len(list_band2)):
                fXM.write("%14.7f %14.7f\n" %  (list_dist[j],list_band2[j][i]-eig_shift2 ) )
            fXM.write("&\n")    
    fXM.close()
    
    commands.getoutput("xmgrace -hardcopy -printfile %s_xm.ps %s.agr" %(stPrefix,stPrefix))
    
    print("Please see result in %s.gp & %s_gp.ps ( gnuplot ) and %s.agr & %s_xm.ps ( xmgrace )" % (stPrefix ,stPrefix,stPrefix,stPrefix))    

def f_band_plot_diff(band_ref,band_new,dir_out):
    '''
    Plot comparison between two bandstructure to a new folder
    '''
    list_band = [band_ref,band_new]
    for band in list_band:
        if (band.kpt.latt is not None):
            latt = band.kpt.latt
            print("Find usable lattice information")
            break

    if (latt is not None):
        for band in list_band:
            if (band.kpt.latt is None):
                band.kpt.latt = latt
                print("Convert lattice information")
                band.kpt.ConvertUnit("cart")

    if (list_band[0].vbm != None):
        f_Band_Plot(list_band,[x.vbm for x in list_band],y_fermi=0.0,dirname_out=dir_out)
    else:
        f_Band_Plot(list_band,[x.fermi for x in list_band],y_fermi=0.0,dirname_out=dir_out)
    return

def f_band_plot_character(band, 
        eig_shift=0, y_fermi=None,
        name_prefix="plotbandc",
        orb_map=None, plot_setting=None):
    '''
    Plot a band structure with band characters indicated 
    Require matplotlib and output graph files

    How orb_map is organized:
    "cond" refers to a dictionary of conditions like "l":1, and only characters with all conditions satisfied are included.
    If a special condition "none" (value not checked) presents, it is the default without any characters.
    "plot" refers to a dictionary that contains "mode"("line" or symbol"), optional "color", "size", "symbol".  "color" and "size" are for 100% pure state and change with characters.

    Only one line and one symbol of each type will be plotted on the final graph and colors are averaged.


    :param list_eig: eigenvalues in [k, band]
    :param list_orb: list of OrbitalNameT 
    :param list_character: list of characters in [k,band,orb]
    :param eig_shift: shift eigenvalues in graph
    :param y_fermi: Fermi level plotted on the graph
    :param name_prefix:  Name prefix
    :param mode: "line" or "symbol", which to indicate band characters
    :param orb_map: a list of dictionary with two keys: "cond" and "plot".
    '''
    print("Plotting...")
    import numpy as np
    from matplotlib import rc
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.font_manager import FontProperties, FontManager
    from matplotlib.markers import MarkerStyle
    from matplotlib.colors import colorConverter
    from matplotlib.lines import Line2D
    from matplotlib.collections import LineCollection
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    if eig_shift is None:
        eig_shift = 0

#Use linecollection makes each segmented a bit shorter
    b_use_linecollection = True
#Number of points between each two k-points, if 1 then no interpolate
    n_interpolate = 1

    size_font = 12
    fontProperties = {'family':'serif', #'serif' : ['Times'],
            'weight' : 'normal', 'size' : size_font}

#Copy setting from plot_setting.font
    if (plot_setting is not None):
        font1 = plot_setting.get("font")
        if (font1 is not None):
            for key, val in font1.items():
                fontProperties[key] = val
        b0 = plot_setting.get("linesegment")
        if (not b0 is None):
            if (b0 == "short"):
                b_use_linecollection = True
            elif (b0 == "long"):
                b_use_linecollection = False
            else:
                raise ValueError("Value of linesegment must be \"long\" or \"short\"")

        n_interpolate = plot_setting.get("interpolatepoints")
        if (n_interpolate is None):
            n_interpolate = 1

    rc('text', usetex=True)
    rc('font', **fontProperties)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    fig_hidden, ax_hidden = plt.subplots()

#Get RGB values from string
    def get_color(color):
        if (isinstance(color,list)): #4-element RGBA array
            return color
        else:
            return colorConverter.to_rgba(color)

    list_dist, list_name, list_break = band.kpt.get_k_xaxis()

#Convert name to latex if necessary
    for i in xrange(len(list_name)):
        name = list_name[i][1]
        if (len(name) > 1 and not name.startswith("$")):#a Greek character
            list_name[i][1] = "$\\%s$" % name.capitalize()

    def check_orb(orb, orb_cond):
        '''
        This used to check whether an orbital follows the rule to be draw
        '''
        for key, val in orb_cond.iteritems():
            if (not hasattr(orb, key)):
                return False
            val2 = getattr(orb, key)
            if (val2 != val):
                return False
        return True

    for prop in band.prop:
#Shift the eigenvalue once 
        list_eig = np.array(prop.eig) - eig_shift
        list_character = np.array(prop.character)[:, :prop.num_band_common]

#Make a map between orbtials and conditions
        list_link = [ [] for x in orb_map]
        for map1 in orb_map:
            mode = map1["plot"]["mode"]
            if (mode != "line" and mode != "symbol"):
                raise ValueError("Unknown plot mode %s" % mode)

        for ix_orb, orb in enumerate(prop.orb):
            for ix_map, map1 in enumerate(orb_map):
                if (map1["cond"].has_key("none")): #Skip default
                    continue
                if (check_orb(orb, map1["cond"])):
                    list_link[ix_map].append(ix_orb)

#plot with characters


#Plot symbols
        n_point = prop.num_band_common * band.num_kpt
        n_max_mix = len(orb_map) #Maximum size of array, + 4 for none
        dic_plotinfo = {} #Each values is two lists: 
#number of colors/size, default color/size (must be 1-d array), list of other color/size, list of weights (1-d arrray)]
#Values for the line plot are stored in key="line"
        symbol_for_line = "line"
        def add_mix(ix, symbol, mix, percent=None, b_default=False):
            '''
            Add an item into dic_color or dic_size
            '''
            dic = dic_plotinfo
            if (map_plot["mode"] == "line"):
                key = symbol_for_line
            else:
                key = map_plot["symbol"]
            if (not dic.has_key(key)):
                dic[symbol] = [
                        [0, None, [], np.empty((n_max_mix, n_point), dtype=np.float64)], #Plot color info
                        [0, None, [], np.empty((n_max_mix, n_point), dtype=np.float64)], #Plot size info
                        len(dic),   #Plot z-order (first on bottom, z-order larger on top)
                        ]

            l1 = dic[symbol][ix]

            if (b_default):
                l1[1] = mix
            else:
                l1[2].append(mix)
                l1[3][l1[0], :] = percent
                l1[0] += 1

        for map1, link1 in zip(orb_map, list_link):
            map_plot = map1["plot"]
            map_cond = map1["cond"]
            if (map_plot["mode"] == "symbol"):
                symbol = map_plot["symbol"]
            elif (map_plot["mode"] == "line"):
                symbol = symbol_for_line
#Add default
            if (map_cond.has_key("none")):
                add_mix(0, symbol, get_color(map_plot["color"]), b_default=True)
                add_mix(1, symbol, [map_plot["size"]], b_default=True)
            else:
                if (map_plot.has_key("size")): #size dependent on characters
                    ar = list_character[:, :, tuple(link1)].sum(axis=2).flatten()
                    add_mix(1, symbol, [map_plot["size"]], percent=ar)
                if (map_plot.has_key("color")): 
                    ar = list_character[:, :, tuple(link1)].sum(axis=2).flatten()
                    add_mix(0, symbol, get_color(map_plot["color"]), percent=ar)
#Generate x values
        x_scatter = np.repeat(list_dist, prop.num_band_common)
        #Calculate the final result
        ar_color = None
        for symbol, info in dic_plotinfo.iteritems():
            if (symbol != symbol_for_line):
                print("Deal with symbol \"%s\"" % symbol)
            else:
                print("Calculate colors for lines...")
            ar = np.zeros( (n_point, 4), dtype = np.float64)
            list_ar = []
            for i in xrange(2): #Color
                if (info[i][0] == 0): #No character mix, just plot it
                    list_ar.append(np.repeat([info[i][1]], n_point, axis=0))
                else:
                    weight_left = 1 - info[i][3][0:info[i][0],:].sum(axis=0)
                    ar1 = np.dot(weight_left[:,np.newaxis], 
                                np.array(info[i][1])[np.newaxis,:])
                    ar2 = np.dot(info[i][3][0:info[i][0]].T, 
                                np.array(info[i][2])) 
                    list_ar.append(ar1 + ar2)

#Split band with breaks
#4 quantities, x (dist) y(eig) color and size
#And interplote
            list_x = [ ]
            list_y = [ ]
            list_color2 = []
            list_size2 = []
            ar_color2 = np.asarray(list_ar[0]).reshape(len(list_dist), prop.num_band_common, 4)
            ar_size2 = np.asarray(list_ar[1][:,0]).reshape(len(list_dist), prop.num_band_common)
            list_break2 = copy.copy(list_break)
            list_break2.append(band.num_kpt)

            for ix_band in xrange(prop.num_band_common):
                list_x1 = []
                list_y1 = []
                list_color1 = []
                list_size1 = []
                for ix1, i1 in enumerate(list_break):
                    i2 = list_break2[ix1+1]
                    if (n_interpolate == 1):
                        x2 = list_dist[i1:i2]
                        y2 = list_eig[i1:i2, ix_band]
                        ar_color = ar_color2[i1:i2, ix_band, :]
                        ar_size = ar_size2[i1:i2, ix_band]
                    else:
                        x2 = np.zeros(((i2-i1-1)*n_interpolate+1,)) 
                        for i3 in range(i1,i2-1):
                            xleft = list_dist[i3]
                            xright = list_dist[i3+1]
                            for i4 in range(0, n_interpolate):
                                x2[(i3-i1) * n_interpolate+i4] = (xright * i4 + xleft * (n_interpolate-i4)) / n_interpolate
                        x2[-1] = list_dist[i2-1]
                        y2 = np.interp(x2, list_dist[i1:i2], list_eig[i1:i2, ix_band])
                        ar_color_v = []
                        for ix in range(4):
                            ar_color_v.append(np.interp(x2, list_dist[i1:i2], ar_color2[i1:i2, ix_band, ix]))
                        ar_color = np.vstack(tuple(ar_color_v)).T
                        ar_size = np.interp(x2, list_dist[i1:i2], ar_size2[i1:i2, ix_band])

                    list_x1.append(np.asarray(x2))
                    list_y1.append(np.asarray(y2))
                    list_color1.append(np.asarray(ar_color))
                    list_size1.append(np.asarray(ar_size))


                list_x.append(list_x1)
                list_y.append(list_y1)
                list_color2.append(list_color1)
                list_size2.append(list_size1)

            if (symbol != symbol_for_line):
                list_x1 = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(list_x))))
                list_y1 = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(list_y))))
                ar_color = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(list_color2))))
                ar_size = list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(list_size2))))
                ax.scatter(list_x1, list_y1,
                        s=ar_size, #we need 1-D but it is 2D with extra []
                        c=ar_color,
                        marker=symbol,
                        edgecolors="none",
                        zorder=info[2]
                        )
            else: #Store to another array
                print("Plotting lines...")
                for ix_band in xrange(prop.num_band_common):
                    for ix_break in range(len(list_x[0])):
                        list_x1 = list_x[ix_band][ix_break]
                        list_y1 = list_y[ix_band][ix_break]
                        ar_color = list_color2[ix_band][ix_break]
                        ar_size = list_size2[ix_band][ix_break]
                        n_k = len(list_x1)

                        if (b_use_linecollection):
                            points = np.array([list_x1, list_y1]).T.reshape(-1, 1, 2)

                            segments = np.concatenate([points[:-1], points[1:]], axis=1)
                            cmap = LinearSegmentedColormap.from_list("cmmy%i" % ix_band, 
                                    ar_color, N=n_k*4)
                            lc = LineCollection(segments, array=np.linspace(0, n_k, 1),
                                    linewidth = ar_size,
                                    cmap=cmap,
                                    norm=Normalize(0, n_k),
                                    zorder=info[2])
                            ax.add_collection(lc)
                        else:
                            for ix_k in range(n_k-1):
                                x1, x2 = list_x1[ix_k:ix_k+2]
                                y1, y2 = list_y1[ix_k:ix_k+2]
         #Colors of lines are based on average of colors of points
                                lw = (ar_size[ix_k] + ar_size[ix_k+1]) / 2
                                color = (ar_color[ix_k] + ar_color[ix_k+1]) / 2

                                ax.plot([x1, x2], [y1, y2], "-", 
                                        linewidth=lw, 
                                        color=color,
                                        zorder=info[2])

#Add Fermi level
    if (y_fermi is not None):
        ax.axhline(y=y_fermi, xmin=0, xmax=list_dist[-1], 
                linestyle = '--', color='black')

#Adjust y range
    if (plot_setting is not None):
        y_min = plot_setting.get("ymin")
        y_max = plot_setting.get("ymax")
    y_min = list_eig.min() - 0.5 if y_min is None else y_min
    y_max = list_eig.max() + 0.5 if y_max is None else y_max
    ax.set_ylim(y_min, y_max)

#Adjust x range
    ax.set_xlim(0, list_dist[-1])
    ax.xaxis.set_major_locator(
        ticker.FixedLocator([list_dist[x[0]] for x in list_name]))
    ax.xaxis.set_major_formatter(
        ticker.FixedFormatter([x[1] for x in list_name]))

#Plot x grid lines
    ax.grid(True, axis="x", linestyle="-")

#Plot Legend
    list_legend = []
    list_legend_label = []
    for map1 in orb_map:
        if (not "none" in map1["cond"]):
            plot1 = map1["plot"]
            label = plot1.get("label")
            if (label is not None):
                list_legend_label.append(label)
                if (plot1["mode"] == "line"):
                    list_legend.append(Line2D([0,1],[0,1],linewidth=plot1["size"], color=plot1["color"]))
                else:
                    list_legend.append(ax_hidden.scatter([0,1],[0,1],
                        marker=plot1["symbol"], c=plot1["color"], s=plot1["size"], edgecolor="none"))

    if (len(list_legend) > 0):
        ax.legend(list_legend, list_legend_label, scatterpoints=1)

                
    

#Save in two formats
    print("Saving figures...")
#Eps is too cumbersome to deal with
#   fig.savefig(name_prefix + ".eps", frameon=False, bbox_inches="tight")
    fig.savefig(name_prefix + ".png", dpi=300, frameon=False, bbox_inches="tight")

    print("Please see results with characters in %s.eps and %s.png" % (name_prefix, name_prefix))

    return

def Main(ArgList):
    description = '''Solve band crossing by left derivative and right derivative and create a new band. Only energy and k-vector information is used. Unit is not considered in calculation.
'''
    usage = "%prog --band BandFileName --kpt kPointFileName -o OutputBandFileName"
    parser = OptionParser(formatter=EmptyFormatter(),usage=usage,description=description)
    parser.add_option("--band",dest="stBandFileName",help="Specify the energy data filename; This file should be a space or tabular split text file with all states in one k-points in one row")
    parser.add_option("--kpt",dest="stKptFileName",help="Specify the k-point list filename; This file should be a space or tabular split text file with x,y,z coordinate of k-vector in one line. The order must be the same as band file")
    parser.add_option("-o",dest="stOutputBandFileName",help="Output band file name")

    (options,args) = parser.parse_args(ArgList)

    if ( len(args)!= 1):
        parser.error("incorrect number of arguments.")
 
    aKPT = KPointsT(options.stKptFileName)
    list_band = f_Band_ReadFromFile(options.stBandFileName)
    f_Band_WriteToFile(f_Band_DetectCrossing(list_band,aKPT.GetTurningPoint()),options.stOutputBandFileName)
    

if __name__ == "__main__":
    Main(sys.argv)    
