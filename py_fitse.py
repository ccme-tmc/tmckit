#!/usr/bin/env python
#This script is used to fit parameters of semi-empirical calculations to given first-principle calculations
#We can fit multiple parameters with multiple compounds.

#To use this scripts, we must prepare something
#* A function that accepts 1-D array and returns a value which should be minimized for scipy.optimize
#* A library, I/O that contains all the parameters for SE and initial values
#* Perform SE calculations with only parameters given, which means other values should be prepared
#* Perform error checking (like RMSD, gap) with only parameters given

import numpy.random
import scipy
from argparse import ArgumentParser

import sys,subprocess,copy,os
import math
import time
import multiprocessing #Warning : any process called by multiprocessing used pickle to transfer paramters, so we must avoid ( at least be cautious ) to use functions that accept class instance as parameters in Pool.async or Pool.map !
import random

from py_utils import EmptyFormatter
from band_utils import f_Band_ReadFromFile,f_Band_GetVBMFromElectron,f_Band_Diff2RMSD,f_Band_GetGap, f_Band_DiffDispersionRMSD
from optimize import wrap_function
from py_xmlio import XMLSerilizer

from common_caseutil import f_Data_ReadMultiCol
from band_utils import BandsT
from tbgw_utils import TBGWPPLibrary , TBGWInput
from yaeh_utils import yaeh_input, yaeh_ReadBand

##@var nProcess
# global parameter to determine how many process can be used
nProcess = 0


class FitSingleInput(object):
    '''
    Input of fitting with single compound 
    This class does not have initial value, so it must be deserilized from xml to use
    '''
    def __init__(self):
#Below part will not be serilized or de-
        ## @var paraBase 
        # Link to database of parameters, which support method "get_para" to return a list as variable parameters
        self.para_base = None
        #:var st_program The name of the program used, inherited from MultiInput
        self.st_program = ""
## @var calcIn
# Describe the input file(s) for the semi-empirical program
# Must support 2 method: write, load_para. We never read when do the fitting so only write is necessary
        self.calcIn = None
        self.bandRef = None
        self.dVBMRef = None
        self.dGapRef = None
        self.listKPtWeight = None
        self.listBandWeight = None
#Below part is used to store result
        self.dGap = None
        self.dRMSD = None
        self.dDispersionRMSD = None
        pass

    def load(self,st_program):
        '''
        Load external content with given parameters
        Note : the program-related input object should be initilized manually
        '''
        self.st_program = st_program
#Customize
        if (st_program == "tbgwpw"):
           self.calcIn = TBGWInput(self.stInputFile)
           self.calcIn.b_warn_write = False
        elif (st_program == "yaeh" or st_program == "yaehg"):
            self.calcIn = yaeh_input(self.stInputFile)

#End Customize

        self.bandRef = f_Band_ReadFromFile(self.stBandRefFile)
        self.listBandWeight = [0.001 for x in range(0,self.nBandStartIndex)] + [ 1.0 for x in range(self.nBandStartIndex,self.nBandStartIndex+self.nBandFitCount)]
        #Additional weight for VBM and CBM if range we considered include VBM and CBM 
        if ( len(self.listBandWeight) >= self.nElectron /2 + 1):
            print("Unoccupied states included in calculation, use different weight for HOMO and LUMO")
            if ( self.dBandHighWeight == None):
                self.dBandHighWeight = 10
            for i in xrange(0,-self.nBandHighVB):
                self.listBandWeight[self.nElectron/2-1-i] = self.dBandHighWeight
            for i in xrange(0, self.nBandHighCB):
                self.listBandWeight[self.nElectron/2] = self.dBandHighWeight - self.dCBMLessWeight 

        if ( self.nBandTotal == None): #Not specified
            #Additional very small weight to make result not fluctuate too much
            self.listBandWeight.append(0.001)
        else:
            #Add all additional band to make result stabilize
            self.listBandWeight = self.listBandWeight + [0.001 for x in range(0,self.nBandTotal-len(self.listBandWeight))]

        self.dVBMRef = f_Band_GetVBMFromElectron(self.bandRef,self.nElectron)
        self.dGapRef,dt,nt = f_Band_GetGap(self.bandRef,self.dVBMRef,False)

#Use weights 
        if (self.bUseKPtWeight):
            self.band = BandsT.load_xml(self.stDataRefFile)
            self.listKPtWeight = [kpt[4] for kpt in self.band.kpt.listKPt]
#Normalize to each k, each band = 1
            d_total = sum(self.listKPtWeight)
            n_len = len(self.listKPtWeight)
            self.listKPtWeight = [ x / d_total * n_len for x in self.listKPtWeight]

    def calc_error(self):
        '''
        Calculate the error based on weights
        '''
        dRMSD = 0 if self.dRMSD is None else self.dRMSD
        dDispRMSD = 0 if self.dDispersionRMSD is None else self.dDispersionRMSD

        if (self.dGapRef is None): #Gap never used
            delta_gap = 0
        elif (self.dGap is not None): #Two gaps
            delta_gap = self.dGap-self.dGapRef
        else: #We need gap but what calculated is metal
            delta_gap = -self.dGapRef
#       print("Error: %s %6.2f %6.2f" % (self.stCompound,self.dRMSD,delta_gap))

        err = (dRMSD**2 + delta_gap**2*self.dGapWeight + dDispRMSD**2*self.dDispersionWeight)**0.5
#       print(dRMSD,  delta_gap, self.dGapWeight, dDispRMSD, self.dDispersionWeight)

        return err

class FitSingleGlobal:
    '''
    Global options for FitSingleInput
    '''
    def __init__(self):
        pass


class FitMultiInput:
    '''
    Input of multi-compound parameter fit
    '''
    def __init__(self):
        self.CalcProgram = None
        self.ParaFileInput=""
        self.ParaFileOutput=""
        self.listCompoundConfig=[]
#Non-serilize part
        self.params_database = None #A class deal with parameters I/O

    def apply_global(self):
        '''
        Put all global parameters into FitSingle
        Everything in FitSingle will be overwritten
        '''
#Get all attributes not None to apply
        print("Applying global options for compounds...")
        for key in dir(self.GlobalOptions):
            if ("__" in key):
                continue
            val = getattr(self.GlobalOptions, key)
            if (val is None or callable(val)):
                continue
            print("... %s ..." % key)
#Apply
            for sIn in self.listCompoundConfig:
                setattr(sIn, key, val)

        return


    def load(self,params_database):
        '''
        Load all single-case data and cache them after sin
        '''
        stCwd = os.getcwd()
#Load Parameters
        self.params_database = params_database
       
#Deal with single calculation in each folder
        for sIn in self.listCompoundConfig:
#Load indivdual file
            stFolder = sIn.stCompound
            print("Set up %s ..."% stFolder)
            os.chdir(stFolder)

            sIn.load(self.CalcProgram)
            sIn.para_base = self.params_database #Link every single-compound calculation to this multi-c
#Clear output
            f = open("single.out","w")
            f.close()

            os.chdir(stCwd)

def f_query_tmpfile(filename):
    '''
    Create a temporary file in /tmp/py_fitse/*, and link filename to i
    Existsed filename will be deleted
    guarantee to not overwrite anyother files
    '''
    tmpdir = "/tmp/fitse"
    tmpname = None

    if (not os.path.exists(tmpdir)):
        os.mkdir(tmpdir)


    while (tmpname is None or os.path.exists(tmpname)):
        tmpname = "/tmp/fitse/%08i" % random.randint(10000000,99999999)
    f = open(tmpname, 'w')
    f.close()

#Note symlink is not a file!
    if (os.path.exists(filename) or os.path.islink(filename)):
        os.remove(filename)
    os.symlink(tmpname, filename)

    return tmpname


def f_CalcSingle_Prepare(sIn,stFolder="."):
    '''
    Prepare input files for yaeh calculation
    :param sIn: the MiniSingle object to indicate calculation
    :param stFolder: the folder where the calculation contains
    '''
    stCurrent = os.getcwd()
    os.chdir(stFolder)

    sIn.calcIn.load_para(sIn.para_base)#To update calculation, update parameter database
    #print("Start calculation under %s" % os.getcwd())
    sIn.calcIn.write(sIn.stOutputFile)
    os.chdir(stCurrent)

    return 0

def f_CalcSingle_Run(stCompound,stInputFile,nElectron,bandRef,dVBMRef,\
        listBandWeight,listKPtWeight,dGapWeight,dDispersionReverseWeight=1,\
        stFolder=".",st_program="unknown", n_process=1, b_plot_band=False):
    '''
    Run calculations
    The parameters are that in a FitSingleInput; to avoid multiprocess problem, every parameter must be a basis types, transfered one-by-one instead of the class
    :param stFolder: the folder where the calculation contains
    :return : a 3-element tuple with Compound Name ( which is used to mark which case is calculated ), Gap and RMSD in this calculation
    '''
    #print(stCompound,stInputFile,nElectron,stFolder)
    stCurrent = os.getcwd()
    os.chdir(stFolder)
#Customize
#   subprocess.check_call("/export/home/wf/programs/TBGW/PW/pw %s &> /dev/null" % stInputFile,shell=True)
    b_error = False
    f = open(os.devnull,'w')
    if (st_program == "tbgwpw"):
        subprocess.check_call(["/export/home/wf/programs/TBGW/PW/pw",stInputFile],stdout=f,stderr=f)
        bandNew = f_Data_ReadMultiCol("pwband.dat",func=float)
    elif (st_program == "yaeh" or st_program == "yaehg"):
        dic_path = {"yaeh" : "/export/home/wf/programs/yaehmop/tightbind/ehbind",
                "yaehg" : "/export/home/wf/programs/yaehmop-org/tightbind/ehbind-mpi"}
#Create symbol link
        list_tmp1 = [stInputFile + x for x in [".out", ".status", ".band"]]
        if (n_process >= 1):
            for i in xrange(1,n_process):
                list_tmp1.append("%s.status.%i" % (stInputFile, i))
                list_tmp1.append("%s.out.%i" % (stInputFile, i))
        list_tmp2 = [f_query_tmpfile(x) for x in list_tmp1]
        try:
            if (n_process <= 1):#Do not use mpirun if specified something < 1
                subprocess.check_call([dic_path[st_program],stInputFile],stdout=f,stderr=f)
            else:
                subprocess.check_call(["mpirun","-np",str(n_process),dic_path[st_program],stInputFile],stdout=f,stderr=f)
            band2 = yaeh_ReadBand("%s.band" % stInputFile)
            bandNew = band2.list_eig
        except subprocess.CalledProcessError:
            print("Cannot run the program with given parameters %s" % stFolder)
            b_error = True
#Delete link and files
        for filename in (list_tmp1 + list_tmp2):
            os.remove(filename)
    else:
        raise ValueError("Unknown Program %s" % st_program)
    f.close()
#End Customize
#
#Return immediately if wrong
    if (b_error):
        dDiff = 1000.0
        dGap = 1000.0
        dDisp = 1000
    else:
        dVBMNew = f_Band_GetVBMFromElectron(bandNew,nElectron)
        if (dVBMNew is None or dVBMRef is None):#Metal
            dVBMNew = bandNew[0][0]
            dVBMRef = bandRef[0][0]
        dDiff = f_Band_Diff2RMSD(bandNew,dVBMNew,bandRef,dVBMRef,listKPtWeight,listBandWeight)
        if ( math.isnan(dDiff) ): #Unphysical result; give it to a fixed value. NaN will stop minimize routine immediately
            print("NaN eigen value @ %s" % stCompound )
            dDiff = 1000
            dDisp = 1000
        else:
            dDisp = f_Band_DiffDispersionRMSD(bandNew, bandRef, listBandWeight, dDispersionReverseWeight)

        #Calculate gap
        f = open("single.out","a") #Single-c case output
        dGap =0.0
        if ( dGapWeight != 0.0):
            dFermi = f_Band_GetVBMFromElectron(bandNew,nElectron)
            dGap,dt,n1 = f_Band_GetGap(bandNew,dFermi,False)
            if ( dGap == None): #Error in get gap, set to -10
                dGap = -10
            f.write("Current Gap: %f\n" % dGap)

        f.write("Current RMSD: %f\n" % dDiff)

        f.close()

    os.chdir(stCurrent)

    return (stCompound, dDiff, dGap, dDisp)

def f_CalcSingle(sIn,stFolder="."):
    f_CalcSingle_Prepare(sIn,stFolder)
    print("Calc %s ..." % sIn.stCompound)
    (stCompound, dDiff, dGap, dDisp) = f_CalcSingle_Run(sIn.stCompound,sIn.stOutputFile,sIn.nElectron,sIn.bandRef,sIn.dVBMRef,sIn.listBandWeight,sIn.listKPtWeight,sIn.dGapWeight,sIn.dDispersionReverseWeight,stFolder,sIn.st_program,1)
    #store result to input object 
    sIn.dGap = dGap
    sIn.dRMSD = dDiff
    sIn.dDispersionRMSD = dDisp
    return sIn

def f_CalcMulti(paraNew,mIn,bSkipCompareOnly=True):
    '''
    Calculate average RMSD of multiple compound
    If new parameter is None, then use parameter currently stored in paraNew
    :param bSkipViewOnly: do not calculate single case with tag "CompareOnly". Should be true for fit and false for compare result display
    :param nProcess: the number of process available in the calculation
    '''
    global nProcess

#Deal with new parameters
    if ( paraNew is not None):
       mIn.params_database.set_para(paraNew)
#   print(paraNew)

#Multi process
#Multi process here is actually not necessary as we invoke yaeh also in extra process, which creates two times of number of yaeh processs.
#However the listening process is rather cheap so do not mind it.
    #print("Processes: %i" % nProcess)
    if ( nProcess != 0):
        pool = multiprocessing.Pool(processes=nProcess)

#Deal with single calculation in each folder
    dDiff = 0.0
    stCwd = os.getcwd()
    listDiff = []

    for sIn in mIn.listCompoundConfig:
        if ( bSkipCompareOnly and sIn.bCompareOnly):
            continue

        stFolder = sIn.stCompound

        #os.chdir(stFolder)
        if (nProcess == 0): #Seqential calculation
            listDiff.append(f_CalcSingle(sIn,stFolder))
        else:
#Note if may difference
            #print("New process...")
            f_CalcSingle_Prepare(sIn,stFolder)
            listDiff.append(pool.apply_async(f_CalcSingle_Run,(sIn.stCompound,sIn.stOutputFile,sIn.nElectron,sIn.bandRef,sIn.dVBMRef,sIn.listBandWeight,sIn.listKPtWeight,sIn.dGapWeight,sIn.dDispersionReverseWeight,stFolder,sIn.st_program,sIn.nProcess)))

        #os.chdir(stCwd)

#Wait until complete
    if ( nProcess != 0):
        pool.close()
        pool.join()
        listDiff = [ x.get() for x in listDiff]
#Update SingleInput in MultiInput, with correct order as when multi-process it is not updated
        listDiff2 = []
        for result in listDiff:
            for sIn in mIn.listCompoundConfig:
                if ( result[0] == sIn.stCompound):
                    sIn.dRMSD = result[1]
                    sIn.dGap = result[2]
                    sIn.dDispersionRMSD = result[3]
                    #Reconstruct listDiff to make it the same as that one in serial mode
                    listDiff2.append(sIn)
                    break
        listDiff = listDiff2

#Calculate average error
    listDiff = [sIn.calc_error() for sIn in listDiff] #SingleInput to error

    dDiff = sum([x**2 for x in listDiff]) #Sum error

    dDiff = (dDiff/len(mIn.listCompoundConfig))**0.5 #Average error
    print("Current multi-compounds RMSD: %f" % dDiff)

    return dDiff

def CallBackFunc(para):
    print("Current parameter: ")
    print(para)

def CallBackFuncGM(para, f, accept):
    '''
    Callback for global minimization in each step
    '''
    print("Current minimum accept %s @ %f" % (accept, f))
    print(para)

def f_PrintCalcMulti(mIn):
    '''
    Print detailed difference of each single compounds
    Note: there are two kinds of weights of bands, one is that in mini.in and another is p-electron only, which is obsolete.
    '''
    b_band_p_only = False
    list_bak = []
    if (b_band_p_only):
        for sIn in mIn.listCompoundConfig:
#set weight of only top valence band and lowest condunction band.
#3/4 of VB ( exclude s from sp only band ) will be set to weight=1 and CBM set to weight=1, otherwise 0
            list_bak.append(sIn.listBandWeight)
            sIn.listBandWeight = [ 0 for x in range(0,sIn.nElectron/8)] + [ 1 for x in range(0,sIn.nElectron/2 - sIn.nElectron/8) ] #VB
            sIn.listBandWeight.append(1) #CBM
            print(sIn.listBandWeight)

#Calculation
    f_CalcMulti(None,mIn,False)

#restore and print result
    print("Error summary:")
    print("%-30s%7s%7s%7s%7s%7s"%("Compound","Gap","Ref","Error","RMSD","Total"))
    i = 0
    for sIn in mIn.listCompoundConfig:
        if (b_band_p_only):
            sIn.listBandWeight = list_bak[i]
        i += 1
        if (sIn.dGapRef is None):
            print("%-30s%7s%7s%7s%7.2f" % (sIn.stCompound,"","","",sIn.dRMSD) )
        else:
            print("%-30s%7.2f%7.2f%7.2f%7.2f%7.2f" % (sIn.stCompound,sIn.dGap,sIn.dGapRef,sIn.dGap-sIn.dGapRef,sIn.dRMSD,sIn.calc_error()) )

    print("Total Error: %9.4f" % (sum([sIn.calc_error()**2 for sIn in mIn.listCompoundConfig])/len(mIn.listCompoundConfig))**0.5)

    return

def Main(ArgList):
    description = '''This program can be used to find optimized orbital parameters for yaeh program to make it approach to standard band structure.
    Several fitting scheme based on band RMSD can be used.
    In general, weight of each band included in fitting is 1, and unused band is 0.001 to avoid fluctuation
    1. Use lowest N band, and set HOMO and LUMO weight to 10 if conduction band is included.
    2. Use N band start from lowest i-th band. Weights of HOMO and LUMO are still 1.
    3. Beside band RMSD, included gap with weight 1.

    Notes on parallization: writing a jobs scheduling system is not easy, so we use a very simple way that N compounds can be computed simultaneously, and after anyone exited another one will start. Also, N_i processes will be used for i-th compound.
    So the processes used at the same time is variable, and may exceed total number of processes assigned, which may reduce the efficiency. To avoid this, set up the number of processes for i-th compound carefully.
    The best way is assign enough process in PBS and run all.
    '''
    global nProcess

    parser = ArgumentParser(description=description)
    parser.add_argument("-p",dest="ProcessCount",type=int,default=0,help="The number of tasks can be executed in parallel. Note p=0 will disable MPI for programs and do not use Pool, but p=1 use Pool")
    parser.add_argument("-m",dest='method_min',type=str,default='simplex',help="The minimization algorithm, possible values: simplex, cg, powell, bfgs, basinhopping")
    parser.add_argument("--tol",dest='tol_min',type=float,default=1e-5,help="The tolerance for minimization, only used in basin-hopping minimization steps now")
    parser.add_argument("--iter",dest="n_iter",type=int,default=100,help="Number of basin-hopping iterations, default 100")
    parser.add_argument("filename",nargs=1,help="The xml input file")

    options = parser.parse_args()

    tStart = time.time()

    nProcess = options.ProcessCount
    if ( nProcess < 0 ):
        nProcess = 0
    if (nProcess != 0):
        print("Parallel tasks: %i" % nProcess)
        if (nProcess == 1):
            print("Warning: only 1 process is used in parallel calculation, note the implentation is different from serial one!")

#Prepare calculation
    xpMulti = XMLSerilizer(globals(),"fit_input.xsd")
    mIn = xpMulti.Deserilize(options.filename[0])
    mIn.apply_global()

#This is customized by programs
#Set parameters I/O
    if (mIn.CalcProgram == "tbgwpw"):
        params_base = TBGWPPLibrary(filename=mIn.ParaFileInput)
        mIn.load(params_base)
        x0 = params_base.get_para()
    elif (mIn.CalcProgram == "yaeh" or mIn.CalcProgram == "yaehg"):
        params_base = yaeh_input(mIn.ParaFileInput,b_write_fix=True)
        mIn.load(params_base)
        x0 = params_base.get_para()
        print(x0)

    MinFunc = wrap_function(f_CalcMulti,tuple([mIn]))

#Check parameters
    nProcessUse = sum([x.nProcess for x in mIn.listCompoundConfig])
    if (nProcess > nProcessUse):
        print("Warning: %i processes are specified when maximally possible parallel process is %i, you can lower to it without performance loss" % (nProcess, nProcessUse)) 

    dic_func_min = {"cg": scipy.optimize.fmin_cg,
            "bfgs": scipy.optimize.fmin_bfgs,
            "simplex" : scipy.optimize.fmin,
            "powell" : scipy.optimize.fmin_powell}

    if ( not mIn.bCompareOnly):
#Global minimization
#We always use powell as the local minimization here
        if (options.method_min == "basinhopping"):
            def take_step(x):
                '''
                Random place a parameters guess, linear random
                a steps are taken as 20% ( if  abs(x) < 1 ) 
                or 20% and at least 2 if (abs(x) > 1)
                Abnormal values (abs(x) > 100) will be shrinked to 1.0
                This is 
                '''
                x2 = []
                for v0 in x:
                    if (abs(v0) < 1 and v0 > 0): #Guess as exponetial part
                        v2 = (numpy.random.random() * 0.8 + 0.6) * v0
                    elif (abs(v0) < 100):
#Change at least this value
                        vd1 = max(abs(v0), 5.0)
                        vd = (numpy.random.random() * 0.8 - 0.4) * vd1
                        v2 = v0 + vd
                    else:
                        v2 = 1.0
                    x2.append(v2)
                print("Jump")
                print(x2)
                return x2


            print("Method: Basin-hopping with tolerance=%f" % options.tol_min)
            res = scipy.optimize.basinhopping(MinFunc, x0=x0, niter=options.n_iter, T=0.5, 
                    minimizer_kwargs = {"method":"Powell","tol":options.tol_min, "callback":CallBackFunc,
                        "options":{'xtol':options.tol_min,"ftol":options.tol_min}},
                    take_step = take_step,
                    callback = CallBackFuncGM)
            res = res.x
        else:
            #Local minimization only
            func_min = dic_func_min[options.method_min]
            if (options.method_min == "cg" or options.method_min == "bfgs"): #1st derivatives (finite difference)
                res = func_min(MinFunc,x0=x0,fprime=None,epsilon=0.001,gtol=1e-02,maxiter=1000,callback=CallBackFunc)
            else:#Value only
                res = func_min(MinFunc,x0=x0,xtol=1e-3,ftol=1e-05,maxiter=10000,callback=CallBackFunc)

#Store result
        mIn.params_database.set_para(res)
        mIn.params_database.write(mIn.ParaFileOutput)

#Show result
    f_PrintCalcMulti(mIn)

    print("Time used: %.2f s"% (time.time()-tStart) )

if __name__ == "__main__":
    Main(sys.argv);
