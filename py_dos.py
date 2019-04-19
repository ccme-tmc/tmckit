#!/usr/bin/env python

import os,sys,copy
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import json

from py_xmlio import XMLSerilizer

from package_utils import  uni_Package_Detect
from common_caseutil import f_WriteFileFloatValue
from dos_utils import OrbitalNameT,f_DOS_CombinePDOS,f_DOS_SortPDOS,f_DOS_SaveData
from vasp_utils import vasp_ReadPDOS
from qesp_utils import QESP_ReadPDOS
from w2k_caseutil import w2k_ReadPDOS



class pdos():
    def __init__(self):
        pass

class energy_values_type():
    def __init__(self):
        pass

class orbital_type():
    def __init__(self):
        pass

class Coord():
    def __init__(self):
        pass


def siesta_read_pdos(stFileName):
    '''
    Read PDOS from siesta pdos.xml
    '''
    #Read XML
    xs = XMLSerilizer(globals(),"siesta_pdos.xsd")
    m = xs.Deserilize(stFileName)
#As SIESTA seperate values in several lines, we combine them

#Read energy first
    listEnergy = []
    for l1 in m.energy_values.__base__:
        listEnergy.extend(l1)
    nCount = len(listEnergy)
#Combine DOS
    listDOS =[ reduce(lambda x,y: x+y,z.data) for z in m.orbital]
    #Get Name
    #listName = [ f_GetOrbitalName("%s%i"%(z.species,z.atom_index),z.n,z.l,z.m) +str(z.z)for z in m.orbital]
    listName = [ OrbitalNameT([z.species,z.atom_index,z.n,z.l,z.m,z.z,0]) for z in m.orbital ]
    dVBM = None



#Deal with spin-polarized cases
    if ( m.nspin == 2):
        listName2 = []
        listDOS2 = []
        for j in range(0,len(listName)):
            #listName2.append(listName[j].replace(")",")(up)"))
            #listName2.append(listName[j].replace(")",")(down)"))
            nameup = copy.copy(listName[j])
            nameup.nspin = 1
            namedown = copy.copy(listName[j])
            namedown.nspin = -1
            listName2.append(nameup)
            listName2.append(namedown)
            DOS1 = []
            DOS2 = []
            for i in range(0,len(listDOS[0])/2):
                DOS1.append(listDOS[j][i*2])
                DOS2.append(listDOS[j][i*2+1])
            listDOS2.append(DOS1)
            listDOS2.append(DOS2)
        listName = listName2
        listDOS = listDOS2

    return listName,listEnergy,listDOS,dVBM

def uni_ReadPDOS(stProgram="auto",para=[]):
    '''
    Read DOS from specific program result
    :param stProgram: program name, "auto" means automatically detected
    :return:  names of states, energy points, DOS,VBM
    '''
    stProgram = uni_Package_Detect(stProgram,para)

#    dicFunc = {"auto":[0,None],"w2k":[0,w2k_ReadBand],"yaeh":[1,yaeh_ReadBand],"siesta":[3,siesta_ReadBand],"vasp":[0,vasp_ReadBand],"qe":[1,QESP_ReadBand]}
    dicFunc = {"auto":[0,None],"siesta":[1,siesta_read_pdos],"vasp":[0,vasp_ReadPDOS],"qe":[1,QESP_ReadPDOS],"w2k":[0,w2k_ReadPDOS]}

    if ( not dicFunc.has_key(stProgram) ):
        raise ValueError,"Package %s is not supported!" % stProgram

    if ( len(para)  != dicFunc[stProgram][0]):
        raise ValueError,"incorrect number of arguments."
    
    #(aKPt,listBand,dVBM) = dicFunc[stProgram][1](*para)
    return  dicFunc[stProgram][1](*para)

def f_pdos_filter(x,st_filter):
    '''
    Test whether a orbital name follows specific condition
    :param x: OrbitalName
    :param st_filter: the string, contains "x" as orbital, will be directly executed !!!
    :return:  boolean
    '''
    result = eval(st_filter)
    return result

def f_pdos_combine(list_level):
    '''
    Wrap multiple combine function
    :param list_level: list of strings that represent combination function level
    '''
    def tmp(x,y):
        '''
        :param x: OrbitalNameT
        :param y: OrbitalNameT
        '''
        if (not OrbitalNameT.comparer_common(x,y)):
            return False
        for level in list_level:
            if (not OrbitalNameT.list_comparer_unit[level](x,y)):
                return False
        return True
    return tmp

def Main(ArgList):
    description = '''Plot DOS of specific package in specific format, also create one-file plain tabular seperated text for DOS.
    File required:
        SIESTA :  pdos.xml or $CASE$.PDOS
        VASP : none, require to run in the case folder
        WIEN2K :  none , require to run in the case folder
        QE :  prefix of output filenames, require to run in the folder contains all output of projwfc.x
    By default, output folder is current folder, so it is sugguested to run this with -d option to redirect output files to another folder as it will create plenty of files.
    Example:
    py_dos.py -d dos --format "%s%i(%l)(%spin)" --filter "x.species=='O'"
    '''
    usage = "%prog -i PROGRAM [program related parameters like $CASE.fdf $CASE.bands] [-g] [-s] [-d OutputFolder] -c [Combine] --filter [Filter]"
    
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("extra_args", type=str, nargs="*", help="Arguments passed to package-specific procedures")
    parser.add_argument("-i",dest="Program",default="auto",help="The package name")
    parser.add_argument("-g",dest="Align",default=False,action="store_true",help="Shift band position to make VBM = 0 in graph")
    parser.add_argument("-s",dest="FormatOnly",default=False,action="store_true",help="Do not plot band structure, just write band structure data and k-point list in space-seperated format")
    parser.add_argument("-d",dest="OutputDir",default="dos",help="The name of folder used to store output files")
    parser.add_argument("-f",dest="ForceOverwrite",default=False,action="store_true",help="Control whether to overwrite if the output directory already exists, default is not ")
    parser.add_argument("--format",dest="Format",default="z",help="Combine multiple DOS into a single one in specific level. Levels include %%spin,%%m,%%l,%%n,%%i(atom index),%%s(species),%%t(total)." + \
    " This string also indicates display names, format strings like %%x will be replaced while others will be kept")
    parser.add_argument("--filter",dest="Filter",default=None,help="The command to judge whether a specific orbital should be included in final results, like 'x.l ==1 and x.i < 70'. Warning: This string will be EXECCUTED DIRECTLY so be cautious!")
    parser.add_argument("--plotinfo", dest="FilePlotInfo", default=None, help="The input file for plotting parameters")

    options = parser.parse_args()

    stProgram = uni_Package_Detect(options.Program, options.extra_args)
    print("Reading PDOS...")
    listName,listEnergy,listDOS,dVBM = uni_ReadPDOS(stProgram, options.extra_args)
    print("Read PDOS OK...")
    print("Parsing PDOS ...")
    listName,listDOS=f_DOS_SortPDOS(listName,listDOS)


#Read plotinfo
    info_plot = None
    if (options.FilePlotInfo is not None):
        with open(options.FilePlotInfo) as f:
            info_plot = json.loads(f.read())

    
#Filter
    if (options.Filter != None):
        list = zip(listName,listDOS)
        list = [x for x in list if f_pdos_filter(x[0],options.Filter)]
        if (len(list) ==0):
            print("All orbitals are filtered! Nothing will be printed.")
            sys.exit(0)
        (listName,listDOS) = zip(*list)
#Combine
#Old-style
#    if ( options.CombineLevel != "z"):
#        listName,listDOS = f_DOS_CombinePDOS(listName,listDOS,options.CombineLevel)
#New-style
#Extract indicators from Format
    b_split_spin = "%spin" in options.Format
    st = options.Format
    list_mark = []
    if (b_split_spin):
        list_mark.append("spin")
        st = st.replace("%spin","")
    i2 = 0
    while (True):
        i2 = st.find("%",i2)
        if (i2 == -1):
            break
        list_mark.append(st[i2+1])
        i2 += 1
    
    func_same = f_pdos_combine(list_mark)
    listName,listDOS = f_DOS_CombinePDOS(listName,listDOS,func_same=func_same)
#Debug line
#   for name in listName:
#       print(name.ToString("z"))

#Reverse spin-down data to negative
    if (b_split_spin):
        print("Note: spin-down data are stored as negative values")
        list = zip(listName,listDOS)
        list = [[x[0],[-y for y in x[1]] if x[0].spin == -1 else x[1]]  for x in list]
        (listName,listDOS) = zip(*list)

        #Convert name
#Create name string
#   listName = [ OrbitalNameT.ToString(x,options.CombineLevel) for x in listName]
    listName = [ OrbitalNameT.ToStringFormat(x,options.Format) for x in listName]

    #Save data
    if (os.path.exists(options.OutputDir) and not options.ForceOverwrite):
        print("Warning:Directory %s already exists, band structure information will not be saved!" % options.OutputDir)
        return

    f_DOS_SaveData(listName,listEnergy,listDOS,dVBM,options.OutputDir,not options.FormatOnly,options.Align, info_plot)
   
if __name__ == "__main__":
    Main(sys.argv);

