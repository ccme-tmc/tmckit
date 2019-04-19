#!/usr/bin/env python

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from common_caseutil import *
from package_utils import uni_Package_Detect
from band_utils import BandsT
import json
import copy

from w2k_caseutil import w2k_ReadBand, w2k_read_band_character, \
        w2k_rescale_character
from yaeh_utils import yaeh_ReadBand
from vasp_utils import vasp_ReadBand, vasp_read_band_character
from qesp_utils import qesp_read_band
from abinit_utils import abi_ReadBand
from wann_utils import wann_read_band
from bdf_utils import bdf_read_band
from yambo_utils import yambo_read_band

#Globally control all plot 
#bAlign = False

def siesta_read_electron(stOutput):
    '''
    Read number of electrons from siesta output
    '''
    stLine = f_ReadStdout("grep -i \"Total number of electrons\" %s"  % stOutput)
    if ( stLine != None):
        if ( stLine.strip() != ""):
            return int(float(stLine.strip().split()[-1]))

    raise ValueError,"%s does not contain total number of electron information (Not a SIESTA output?)" % stOutput


def siesta_read_band(stFDFFileName,stOutputfile,stBandFileName):
    '''
    Read band structure of siesta
    Number of electrons is deduced from fermi energy, which may give wrong value when system is a metal
    @todo when spin=2 , two bands does not seperate
    '''

    dicScale = {"pi/a":"piba","ReciprocalLatticeVectors":"crystal"}

    aKPt = KPointsT()
    
    f=open(stFDFFileName)
    list_stLine = f.readlines()
    nInfoLine = 0 #The information line in .bands file is 3 if BandLines is used, 2 if BandPoints is used
    f.close()
    for i in range(0,len(list_stLine)):
        stLine = list_stLine[i].strip()
        if ( len(stLine) == 0 ): #Skip empty line
            continue
        if ( stLine[0] == '#'): #Skip comments
            continue
        if ( "BandLinesScale " in stLine):
            aKPt.stMode = dicScale[stLine.split()[1]]
        if ( "%block BandLines" in stLine):
            j = 1
            listPt = []
            listCount = []
            arLine = list_stLine[i+j].split()
            while(not "%endblock" in arLine[0]):
                listPt.append(arLine[1:])
                listCount.append(int(arLine[0]))
                j += 1
                arLine = list_stLine[i+j].split()
            aKPt.CreateFromPointLine(listPt,listCount[1:])
            nInfoLine = 3
            break
        elif ( "%block BandPoints" in stLine):
            j = 1
            listPt = []
            arLine = list_stLine[i+j].split()
            while ( not "%endblock" in arLine[0]):
                listPt.append( [float(x) for x in arLine[:3]])
                j += 1
                arLine = list_stLine[i+j].split()
            aKPt.ReadFromList(listPt)
            nInfoLine = 2
            break

    #Convert aKPt to tpiba mode data if it is piba
    if ( aKPt.stMode == "piba"):
        aKPt.stMode = "tpiba"
        for ar in aKPt.listKPt:
            ar[1:4] = [ x * 0.5 for x in ar[1:4]]

    nElectron = siesta_read_electron(stOutputfile)

    #read band structure
    f=open(stBandFileName)
    list_stLine = f.readlines()
    f.close()
    

    fFermi = float(list_stLine[0])
    (nBand,nSpin,nKPt) = [int(x) for x in list_stLine[nInfoLine].split()]

    if ( nKPt != len(aKPt.listKPt)):
        raise ValueError,"k-points in %s and %s are not same!" % ( stFDFFileName,stBandFileName)

    listBand = []
    (dVBM,fCBM) = [float(x) for x in list_stLine[nInfoLine-1].split()]

    #read energy until the count of energy reaches nBand
    arLine = []
    i = nInfoLine + 1
    nSkip = 1 if nInfoLine == 3 else 3 # Skip first term( band line) if BandLines is used, Skip first three term (k-pt coord) if BandPoints is used
    delta = 0.0002
    #nHOIndex = -1
    while( len(listBand) < nKPt ):
        arLineNow = [float(x) for x in list_stLine[i].split()] #the first 10 is preserved for k-point coordinate
        if ( len(arLine) != 0): #a start one k-point
            arLine += arLineNow
        else:
            arLine = arLineNow[nSkip:]
#Cut
        if ( len(arLine) == nBand*nSpin): #Read one k-point end
            listBand.append(arLine)
            dVBM = max([x for x in arLine if x < fFermi+delta]+[dVBM])
            fCBM = min([x for x in arLine if x > fFermi-delta]+[fCBM])
            arLine = []
            #if ( dVBM in listBand[-1]):
                #nHOIndex = listBand[-1].index(dVBM)
                #print("VBM:%i %f" % ( nHOIndex,dVBM))
                
        i += 1
    
    #use VBM as zero point; if either VBM or CBM is very close to fFermi, then use fFermi
    if ( abs(dVBM-fFermi) < delta/2 or abs(fCBM-fFermi) < delta/2):
        dVBM = fFermi

    print("Detect VBM from SIESTA Fermi: %f" % dVBM)

#    return aKPt,listBand,dVBM,nElectron,nSpin
    return BandsT(aKPt,listBand,dVBM,None,nElectron,nSpin)

def read_eigs(filename):
    '''
    Read eigenvalues only as a band structure
    '''
    list_eig = f_Data_ReadMultiCol(filename,func=float)
    kpt = KPointsT()
    kpt.listKPt = [["",0,0,x,0] for x in xrange(len(list_eig))]
    band = BandsT(kpt=kpt, list_eig=list_eig)
    return band

def uni_ReadBand(stProgram="auto",para=[]):
    '''
    Read band structure from specific program result
    :param stProgram: program name, "auto" means automatically detected
    :param para: arguments for specific programs
    :return: BandsT object
    '''
#   :param nElectron: electron number, default 0. If specified, then VBM will be calculated for plot
#   :param bPlot: Whether plot bandstructure
#   :param bAlign: Whether shift band energy to make VBM = 0 in plot
#    :return:  k-points,eigenvalues,VBM,Electron count,Spin count
#    '''

    stProgram = uni_Package_Detect(stProgram,para)

    dicFunc = {\
            "auto":[(0,),None],\
            "w2k":[(0,1,),w2k_ReadBand],\
            "yaeh":[(1,),yaeh_ReadBand],\
            "siesta":[(3,),siesta_read_band],\
            "vasp":[(0,),vasp_ReadBand],\
            "vasp_out":[(1,),vasp_ReadBand],\
            "qe":[(1,), qesp_read_band],\
            "abinit":[(2,),abi_ReadBand],\
            "gap":[(1,),w2k_ReadBand],\
            "w90":[(1,),wann_read_band],\
            "bdf":[(1,),bdf_read_band],\
            "yambo":[(1,2,),yambo_read_band],\
            "raw":[(1,),lambda x:BandsT(dirname=x)],\
            "band":[(1,),read_eigs]\
            }

    if (not dicFunc.has_key(stProgram) ):
        raise ValueError,"Package %s is not supported!" % stProgram

    if (not len(para) in dicFunc[stProgram][0]):
        raise ValueError,"incorrect number of arguments : %i" % len(para)
    
    #(aKPt,listBand,dVBM) = dicFunc[stProgram][1](*para)
    return  dicFunc[stProgram][1](*para)

def uni_read_band_character(stProgram="auto",para=[]):
    '''
    Read band character informations

    :return: two lists, first is list of components in OrbitalT, second is [kpt,band] nested array of all components
    '''
    stProgram = uni_Package_Detect(stProgram,para)

    dicFunc = {
            "auto":[(0,),None],
            "w2k":[(0,),w2k_read_band_character],
            "vasp":[(0,),vasp_read_band_character]
            }

    if ( not dicFunc.has_key(stProgram) ):
        raise ValueError,"Package %s is not supported!" % stProgram

    if ( not len(para) in dicFunc[stProgram][0]):
        raise ValueError,"incorrect number of arguments."

    return dicFunc[stProgram][1](*para)

def f_band_remove_k(band1, n):
    '''
    Remove first n k-points from bands
    '''
    if (n == 0):
        return band1

    band = copy.deepcopy(band1)

    band.list_eig = band.list_eig[n:]

    for prop in band.prop:
        prop.eig = prop.eig[n:]
        if (prop.character is not None):
            prop.character = prop.character[n:]
        if (prop.occ is not None):
            prop.occ = prop.occ[n:]

    band.kpt.listKPt = band.kpt.listKPt[n:]

    return band



def Main(ArgList):
    description = '''Plot band strcture of specific package in specific format, also create plain tabular seperated text of band structure data and k-points list.
    File required:
        SIESTA : .fdf , stdout, .bands
        VASP : none, require to run in the case folder
        WIEN2K :  none , require to run in the case folder
        FHI-GAP : energy file name suffix (_gw-band or _gw0-band) , require to run in the WIEN2k case folder
        QE : data-file.xml in the save folder, or pw.x output
        ABINIT : _EIG output or _GW output
        Wannier90: _w90_band.dat file ( other files are necessary but will be selected according to this file name)
        Yambo : o-*.qp file
        RAW: folder name contains all files generated by %prog 
    By default, output folder is current folder, so it is sugguested to run this in an empty folder to avoid overwritten as it will create plenty of files.
    Band characteres can also be plotted. To do this add "--plotinfo File". The input file contains extra format informations.
    '''
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    parser.add_argument("extra_args", type=str, nargs="*", help="Arguments passed to package-specific procedures")
    parser.add_argument("-i",dest="Program",default="auto",help="The package name")
    parser.add_argument("-g",dest="Align",default=False,action="store_true",help="Shift band position to make VBM = 0 in graph")
    parser.add_argument("--electron",dest="Electron",default="0",help="The number of electrons in the system of the band structure. This is used to determine VBM/Fermi energy. If not specified, this program will try to find from input files")
    parser.add_argument("-s",dest="FormatOnly",default=False,action="store_true",help="Do not plot band structure, just write band structure data and k-point list in space-seperated format")
    parser.add_argument("-d",dest="OutputDir",default="band",help="The name of folder used to store output files")
    parser.add_argument("-k",dest="KPointsUnit",default="cart",help="The unit of k-points in the output, also determine x-axis distance. Possible options include 'cart','crystal' and 'default', where default means use what the program read.")
    parser.add_argument("-f",dest="ForceOverwrite",default=False,action="store_true",help="Control whether to overwrite if the output directory already exists, default is not ")
    parser.add_argument("--cross",dest="ResolveCorssing",default=False,action="store_true",help="Contrl whether to resolving band-crossing by derivatives")
    parser.add_argument("--plotinfo", dest="FilePlotInfo", default=None, help="The input file for plotting parameters, include band structure with band characters indicated")
    parser.add_argument("--copyk", dest="CopyKPoints", default=None, help="Copy k-points information from given band.xml, instead of use what is read. This to make two plots have the same x-coord")
    parser.add_argument("--removek", dest="RemoveKPoints", type=int, default=0, help="Remove the first N k-points in the bandstructure, especially useful when plotting HSE / metaGGA bands that must be calculated with full M-P k-points self-consistently")
    parser.add_argument("--kname", dest="NameK", default=None, help="Provide a list of names of k-points that will be plotted, in format \" 1 Gamma 30 X 50 L\", indicies start from 1, this works before --removek")
    parser.add_argument("--alignenergy", dest="AlignEnergy", default=None, help="Manually specify the energy to align instead of using VBM/Fermi energy automatically read. The unit must be the same as what the package outputs.") 

    options = parser.parse_args()

    stProgram = uni_Package_Detect(options.Program,options.extra_args)
    band = uni_ReadBand(stProgram,options.extra_args)
#Read k-points
    if (options.CopyKPoints is not None):
        band_copyk = BandsT.load_xml(options.CopyKPoints)
        band.kpt = band_copyk.kpt

    if (options.NameK is not None):
        list_name = options.NameK.split()
        for i in xrange(len(list_name) / 2):
            ix1 = int(list_name[i*2]) - 1
            name1 = list_name[i*2+1]
            band.kpt.listKPt[ix1][0] = name1

#Read character informations
    info_plot = None
    if (options.FilePlotInfo is not None):
        with open(options.FilePlotInfo) as f:
            info_plot = json.loads(f.read())
        if (info_plot.has_key("orb_map")):
            list_orb, list_character = uni_read_band_character(stProgram)
            if (stProgram == "w2k"): #Wien2K qtl is not good for plotting
                list_character = w2k_rescale_character(band.list_orb, band.list_character)
#Add to band object
            for prop, orb, character in zip(band.prop, list_orb, list_character):
                prop.orb = orb
                prop.character = character

#Remove extra k-points
    if (options.RemoveKPoints > 0):
        band = f_band_remove_k(band, options.RemoveKPoints)

#Resolving crossing
    if (options.ResolveCorssing):
        band.resolve_crossing()

#Unit conversion
    unit_k = options.KPointsUnit.lower()
    if (unit_k != "default" and band.kpt.stMode != unit_k):
        print("Convert k-point unit from %s to %s" % (band.kpt.stMode,unit_k))
        #band.kpt.ConvertUnit(unit_k,band.kpt.latt.PrimitiveCellVector)
        band.kpt.ConvertUnit(unit_k)

    nElectronInput = int(options.Electron)
    if ( nElectronInput != 0):
        band.num_electron = nElectronInput

#Align energy
    if (options.AlignEnergy is not None):
        e1 = float(options.AlignEnergy)
        band.fermi = e1
        band.vbm = e1

#   f_Band_Analyse(listBand)
    band.show_info( (band.fermi if band.b_metal else band.vbm) if options.Align else None)

#   f_Band_SaveData(aKPt,listBand,dVBM,options.OutputDir,nElectron,not options.FormatOnly,options.Align)
    if (os.path.exists(options.OutputDir) and not options.ForceOverwrite):
        print("Warning:Directory %s already exists, band structure information will not be saved!" % options.OutputDir)
        return
    
    band.save(options.OutputDir,b_plot=not options.FormatOnly,align=BandsT.align_auto if options.Align else BandsT.align_none, plotinfo=info_plot)

if __name__ == "__main__":
    Main(sys.argv);
