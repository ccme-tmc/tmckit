#!/usr/bin/env python
#The module to deal with dos
import os
import subprocess 
from math import *

from constants import Ry2eV,Ha2eV,Pi
from common_caseutil import f_WriteFileFloatValue
from py_modagr import set_property_auto, split_multiline_to_single

def f_DOS_FromEigens(listBand,listWeight=None,dMin=None,dMax=None,dStep=0.01,dGauss=0.1,stUnit="eV"):
    '''
    Calculate DOS from eigenvalues of each k-point and weight
    :param listBand: the eigenvalues, one array for each k-point
    :param listWeight: the weight of k-points 
    :param dMin: the minimum value of DOS range, default the minimum-dGauss*6 in eigenvalues.
    :param dMax: the maximum value of DOS range, default the maximum+dGauss*6 in eigenvalues.
    :param dStep: the distance between two data points
    :param dGauss: the Gaussing smearing
    :param stUnit: the unit of listBand, default eV, possible values include "Ry","Rydberg","Ha" and "Hartree". All other parameters are always assumed to be eV. 
    :return: list of E-DOS data in 2 lists E and DOS. If unit is specified, then output unit is eV.
    '''
    if ( listWeight is not None):
        if ( len(listBand) != len(listWeight)):
            raise ValueError,"The number of k-points is different between eigenvalues and weights"

#Convert unit
    if ( stUnit == "Ha"  or stUnit == "Hartree"):
        listBand2 = [ [ x*Ha2eV for x in y] for y in listBand]
    elif ( stUnit == "Ry" or stUnit == "Rydberg"):
        listBand2 = [ [ x*Ry2eV for x in y] for y in listBand]
    else:
        listBand2 = listBand

#Set range
    if ( dMin == None):
        dMin = min([min(x) for x in listBand2])
        dMin = int(dMin-1-6*dGauss)
    if ( dMax == None): 
        dMax = max([max(x) for x in listBand2])
        dMax = int(dMax+1+6*dGauss)


    dCoef = 1.0/2/dGauss/dGauss
    dCoef2 = 1/dGauss/sqrt(2*Pi)

    dEnergy = dMin
    listResult = []
    while dEnergy <= dMax:
        dTot = 0.0
        for i,listEnergy in enumerate(listBand2):
            dTotK = 0.0
#Summation  all eigenvalues contribution in one k
            for d1 in listEnergy:
                dTotK += dCoef2 * exp(-dCoef*(d1-dEnergy)**2)
#Multiply weight
            if ( listWeight is not None):
                dTotK *= listWeight[i]
#Sum k*weight
            dTot += dTotK
        listResult.append([dEnergy,dTot])
#Next energy point
        dEnergy += dStep
    return [x[0] for x in listResult],[x[1] for x in listResult]

def f_DOS_SaveData(listName,listEnergy,listDOS,dVBM=None,stOutDir=None,bPlot=True,bAlign=True, plotinfo=None):
    '''
    Save DOS related data to specific folder

    :param stOutDir: the directory to put generated files. If set to None, will use current directory
    :param plotinfo: an object represent how to control output
    '''
    if ( stOutDir == None):
        stOutDir = os.getcwd()

    if ( not os.path.exists(stOutDir) ):
       os.mkdir(stOutDir)
   
    stDirCur = os.getcwd()
    os.chdir(stOutDir)

    f = open("dos.dat","w")
#Write name
    f.write("#Energy")
    for st in listName:
        f.write("\t%s" % st)
    f.write("\n")
#Write energy-DOS
    for i,dEnergy in enumerate(listEnergy):
        f.write("%10.6f" % dEnergy)
        for j in range(0,len(listName)):
            f.write("\t%10.6f" % listDOS[j][i] )
        f.write("\n")

    if ( dVBM is not None):
        f_WriteFileFloatValue(dVBM,"dos.fermi")

    if ( bPlot ):
        plot_setting = None if plotinfo is None else plotinfo.get("plot_setting")
        f_DOS_Plot([x for x in listName],listEnergy,listDOS, dVBM, bAlignData=bAlign and (dVBM is not None), plot_setting=plot_setting)
   
    f.close()

    os.chdir(stDirCur)


class OrbitalNameT:
    '''
    The class to represent the notation of an orbital. possible cases included z,m,l,n,i(atom index),s(atom species).
    Also the name can be specified explicitly with a string.
    If any of l, m, z not specied, then it is treated as summation of all lowers
    '''
    qn_none = -999 #Quantum number which represents not specified
    dic_spin = {-1:"down",0:"",1:"up"}
    dic_l = {-1:"",
            0:"s",1:"p",2:"d",3:"f",4:"g",-999:"X"}
    dic_l_rev = {"s":0,"p":1,"d":2,"f":3,"g":4}
#Deprecated
    comparer_default = lambda x,y : not x.is_custom and not y.is_custom and x.spin == y.spin
    comparer_custom = lambda x,y : not x.is_custom and not y.is_custom 
    comparer_common =  lambda x,y : not x.is_custom and not y.is_custom and x.total == y.total and x.inter == y.inter
    #Function compare whether two objects have same quantum number
    list_comparer_unit = {
            "z" : lambda x,y: False,
            "m" : lambda x,y: x.m_name == y.m_name if x.m_name is not None else (x.m == y.m) ,
            "l" : lambda x,y: x.l == y.l,
            "n" : lambda x,y: x.n == y.n,
            "i" : lambda x,y: x.i == y.i,
            "s" : lambda x,y: x.species == y.species,
            "t" : lambda x,y: x.total == y.total,
            "spin" : lambda x,y: x.spin == y.spin
    }
#Deprecated
    list_comparer = {
            "z" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and False,
            "m" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and (x.m_name == y.m_name if x.m_name is not None else x.m == y.m) and x.l == y.l and x.n == y.n and x.i == y.i and x.species == y.species,
            "l" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and x.l == y.l and x.n == y.n and x.i == y.i and x.species == y.species,
            "n" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and x.n == y.n and x.i == y.i and x.species == y.species,
            "i" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and x.i == y.i and x.species == y.species,
            "s" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and x.species == y.species,
            "t" : lambda x,y:  OrbitalNameT.comparer_default(x,y) and x.total == y.total
            }
    dic_level = { 
            "z" : 2,
            "m" : 3,
            "l" : 4,
            "n" : 5,
            "i" : 6,
            "s" : 7,
            "t" : 8,
            } #The output level


    def __init__(self,content=[],name=None):
        '''
        :param content: The content of indexes, must be an array in the order of s,i,n,l,m,z
        '''
        self.name = "" #custom name. If it is not empty then always use it in later processing
        self.spin = 0 #spin
        self.z = OrbitalNameT.qn_none #zeta
        self.m = OrbitalNameT.qn_none #m
        self.l = OrbitalNameT.qn_none #l
        self.n = OrbitalNameT.qn_none #n
        self.i = OrbitalNameT.qn_none #atom index
        self.species = "" #atom species(This should be a string)
        self.total = False #Is it the total of whole system?
        self.m_name = None # Custom string for m. This function is designed as m may be arbitrally choosed by a program
        self.inter = False #Is it belong to regions outside of RMT, interstitial region, in WIEN2k ?
        if ( len(content) == 7):
            (self.species,self.i,self.n,self.l,self.m,self.z,self.spin) = content
        if ( name is not None):
            self.name= name

    @property
    def is_custom(self):
        '''
        Whether this name is customized, which means the qunatum numbers are unknown.
        '''
        return (self.name != "" and self.name is not None)

    def ToStringFormat(self,format):
        '''
        Return the string to represent the name
        :param format: Format string, %x is used to represent information of orbitals. Possbile choices are %z, %m, %l, %n, %i, %s, %t, %spin(spin)
        '''
        if ( self.is_custom ):
            return self.name
#Special for total DOS
        if (self.total):
            if ("%spin" in format):
                st = "tot%s"%(OrbitalNameT.dic_spin[self.spin])
            else:
                st = "tot"
            return st
        st = format
        st = st.replace("%spin",OrbitalNameT.dic_spin[self.spin])
        st = st.replace("%z",str(self.z))
        st = st.replace("%m",self.m_name if self.m_name is not None else str(self.m))
        st = st.replace("%l",OrbitalNameT.dic_l[self.l])
        st = st.replace("%n",str(self.n))
        st = st.replace("%i",str(self.i))
        st = st.replace("%s",self.species)
        st = st.replace("%t","sum")
        return st

    def __str__(self):
        return self.ToString("z")

    def __repr__(self):
        return self.__str__()

    def ToString(self,stLevel):
        '''
        Return the string to represent the name
        :param stLevel: indicate to which quantum number should be printed. Possible level included z,m,l,n,i,s (spin is controlled by its number)
        '''
        if ( self.is_custom ):
            return self.name

        if ( self.total): #Total: ignore level
            stName = "tot%s"%(OrbitalNameT.dic_spin[self.spin])
            return stName

#If a value is -1 it is allowed to be skipped
        st_z = str(self.z) if self.z != OrbitalNameT.qn_none else ""
        st_n = str(self.n) if self.n != OrbitalNameT.qn_none else ""
        st_i = str(self.i) if self.i != OrbitalNameT.qn_none else ""
#m maybe negative
        if ( self.m_name is not None):
            st_m = self.m_name
        else:
            st_m = str(self.m) 

        #String in bracket
        #if it is empty, then bracket part is neglected
        st2 = ""

        if ( stLevel == "z"):
            st2 = "%s%s%s"%(st_n,OrbitalNameT.dic_l[self.l],st_m)
        elif ( stLevel == "m"):
            st2 = "%s%s%s"%(st_n,OrbitalNameT.dic_l[self.l],st_m)
        elif( stLevel == "l"):
            st2 = "%s%s"%(st_n,OrbitalNameT.dic_l[self.l])
        elif( stLevel == "n"):
            st2 = "%s"%(st_n)

        if ( st2 != ""):
            st2 = "(%s)" % st2

        if ( stLevel == "z"):
            stName = "%s%s%s%s"%(self.species,st_i,st2,st_z)
        elif ( stLevel == "m"):
            stName = "%s%s%s"%(self.species,st_i,st2)
        elif( stLevel == "l"):
            stName = "%s%s%s"%(self.species,st_i,st2)
        elif( stLevel == "n"):
            stName = "%s%s%s"%(self.species,st_i,st2)
        elif( stLevel == "i"): 
            stName = "%s%s"%(self.species,st_i)
        elif( stLevel == "s"):
            stName = "%s%s"%(self.species,st_i)
        elif( stLevel == "t"):
            stName = "sum"
        else:
            raise ValueError,"Unknown orbital level %s" % stLevel

#Add spin
        stName += OrbitalNameT.dic_spin[self.spin]

        return stName

def f_DOS_CombinePDOS(listName,listDOS,level=None,func_same=None):
    '''
    Combine multi-zeta dos to one DOS
    :param level: to which level the DOS should be combined. "m" means only qunatum number m is used to distinguish (of course i,n,l are considered) but zeta is not, and "l" means only l, "n" means only n, "i" means only atom,"s" means atom species,"t" means only total. 
    :param func_same: Custom function to determine whether two orbitals are same and should be combined. 
    '''
    if (func_same is None):
        func_same = OrbitalNameT.list_comparer[level]

    listDup = []
#Detect whether two states belongs to one
    for i,stName in enumerate(listName):
        if ( len(listDup) == 0):
            listDup.append([i])
            continue
        #stName2 = stName[:-1]#Get name without zeta information
        bFound = False
        for dup in listDup:
            for j in dup:
                #if ( stName2 == listName[j][:-1]):#Same 
                if ( func_same(stName,listName[j])):
#                   print("Same: %s / %s" % (stName, listName[j]))
                    dup.append(i)
                    bFound = True
                    break
            if ( bFound):
                break
        if ( not bFound ):
            listDup.append([i])
#Combine DOS
    listName2 = []
    listDOS2 = []
    #print(listDup)
    for dup in listDup:
        listName2.append(listName[dup[0]])
        listDOS2.append( [sum(y) for y in (zip(*[ listDOS[x] for x in dup]))])

    return listName2,listDOS2

def f_DOS_CompareName(name1,name2):
    '''
    Compare the order of two OrbitalNameT object for sorting
    '''
    if ( name1.is_custom or name2.is_custom ):#Do not sort if anyone is customized
        return 0
    if ( name1.total == True):
        return -1
    elif ( name2.total == True):
        return 1

    return cmp([name1.i,name1.n,name1.l,name1.m,name1.z,name1.spin],
            [name2.i,name2.n,name2.l,name2.m,name2.z,name2.spin],
            )

def f_DOS_SortPDOS(listName,listPDOS):
    '''
    Sort PDOS data by name
    '''
    list = zip(listName,listPDOS)
    list.sort(cmp=lambda x,y:f_DOS_CompareName(x[0],y[0]))
    return zip(*list)


def f_DOS_Plot(listName,listEnergy,listDOS,dFermi=None, bAlignData=False, dirname_out=None, name_prefix="plotdos", plot_setting=None):
    '''
    Plot density of states data from states names and raw data.

    :param listName: the name of states ( mostly orbital names like Ta(1)px
    :param listEnergy: the list of energy ( x-axis)
    :param listDOS: the list of DOS data, first dimension as listName, and second dimension as listEnergy
    :param dFermi: the Fermi energy
    :param name_prefix: the file name prefix
    :param bAlignData: whether the energy should be shift to make Fermi = 0 
    :param plot_setting: An object to control how to display plots
    '''
    if (dirname_out is None):
        filename_full = name_prefix
    else:
        if (not os.path.exists(dirname_out)):
            os.mkdir(dirname_out)
        filename_full = os.path.join(dirname_out,name_prefix)

    lines = f_dos_plot_xmgrace(listName, listEnergy, listDOS, dFermi, bAlignData, plot_setting)

    with open("%s.agr" % filename_full, 'w') as fXM:
        fXM.writelines(lines)
        
    try:
        subprocess.check_call(["xmgrace", "-hardcopy", "-printfile", "%s_xm.ps" % filename_full, "%s.agr" % filename_full])
        print("Please see result in %s.agr & %s_xm.ps ( xmgrace )" % (name_prefix,name_prefix))
    except subprocess.CalledProcessError:
        print("Cannot invoke xmgrace to generate pictures")
        print("Please see result in %s.agr" % (name_prefix))

    return



def f_dos_plot_xmgrace(listName, listEnergy, listDOS, dFermi=None, bAlignData=False, plot_setting=None):
    '''
    Plot density of states data from states names and raw data.

    :param listName: the name of states ( mostly orbital names like Ta(1)px
    :param listEnergy: the list of energy ( x-axis)
    :param listDOS: the list of DOS data, first dimension as listName, and second dimension as listEnergy
    :param dFermi: the Fermi energy
    :param bAlignData: whether the energy should be shift to make Fermi = 0 
    :param plot_setting: An object to control how to display plots
    '''
    if ( bAlignData ): #fermi level for output
        if ( dFermi == None):
            raise ValueError,"Fermi level must be specified to align energy values"
        dFermi1 = 0
        dShift = dFermi
    else:
        dFermi1 = dFermi
        dShift = 0

    Emin = min(listEnergy)-1
    Emax = max(listEnergy)+1

    if ( bAlignData ):
        Emin = Emin - dShift
        Emax = Emax - dShift

    xmin = Emin-0.5
    xmax = Emax+0.5
    ymin = min([min(x) for x in listDOS])
    if ( ymin > 0.0 ):
        ymin = 0.0
    ymax = max([max(x) for x in listDOS])

#Create xmgrace input file
    lines = []
    lines.append('''# Grace project file
#
@version 50000
''')
    
    #Fermi level
    if ( dFermi1 is not None):
        lines.append('''@with line
@    line on
@    line g0    
@    line loctype world
@    line %f,%f,%f,%f
@line def
''' % (dFermi1,ymin,dFermi1,ymax) )
        lines.append('''@with string
@    string on
@    string g0
@    string loctype world
@    string %f,%f
@    string def "E\sF"
''' % (dFermi1,ymin-(ymax-ymin)/20 ))

#Graph format
    lines.append('''@with g0
@    world %f,%f,%f,%f
@    title "DOS"
@    autoticks
@    xaxis  label char size 1.0
@    xaxis  ticklabel char size 1.0
@    xaxis  label "Energy(eV)"
@    yaxis  label "Density of States (arb. unit)"
@    yaxis  tick off 
@    legend on
''' % (xmin,ymin,xmax,ymax))


    #Line format 
    for i in range(0,len(listName)):
        lines.append("@    s%d line color %d\n" % (i,i%15+1 ))
        lines.append("@    s%d legend \"%s\"\n" % (i,listName[i]))
    
    #Data
    for i in range(0,len(listName)):
        lines.append("@target G0.S%d\n" % i)
        lines.append("@type xy\n")
        for j in range(0,len(listEnergy)):
            lines.append("%14.7f %14.7f\n" %  (listEnergy[j] - dShift, listDOS[i][j]))
        lines.append("&\n")    

#lines may contain multi-line string, split them
    lines2 = split_multiline_to_single(lines)

#Adjust the plotting
    if (plot_setting is not None):
        lines2 = set_property_auto(lines2, dic_property=plot_setting)

    return lines2

