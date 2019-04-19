#!/usr/bin/env python
#This module is used to plot band structure and DOS from raw data

from __future__ import print_function
from optparse import OptionParser
import sys,os,copy,re,commands
import stat

from common_caseutil import KPointsT,f_ReadFileFloatValue
from py_utils import EmptyFormatter
from band_utils import f_Band_ReadFromFile,f_PlotBand
import list_utils as lu



def Main(ArgList):
    description = '''Plot band structure or DOS from data including kpoint-energy data, kpoint list, kpoint name list ( for band )  or energy-DOS data, PDOS name ( for DOS). Fermi energy can be used to adjust position in the process.
    %prog --fermi 0.5
'''
    usage = "%prog "
    parser = OptionParser(formatter=EmptyFormatter(),usage=usage,description=description)
    parser.add_option("--band",dest="BandFileName",help="Specify the energy data filename; This file should be a space or tabular split text file with all states in one k-points in one row")
    parser.add_option("--kpt",dest="KptFileName",help="Specify the k-point list filename; This file should be a space or tabular split text file with x,y,z coordinate of k-vector in one line. The order must be the same as band file")
    parser.add_option("--fermi",dest="FermiFileName",help="Specify the fermi energy filename.")
    parser.add_option("-g",dest="Align",action="store_true",default=False,help="If set to true, data will be modified to make E_F=0, otherwise a fermi level will be pointed out at its actual position")
    parser.add_option("--band2",dest="Band2FileName",default=None,help="Specify the second energy data filename; the second energy data file will be plotted in dot mode")
    parser.add_option("--fermi2",dest="Fermi2FileName",help="Specify the second fermi energy filename.")
    parser.add_option("--ymax",dest="Ymax",default=None,help="Maximum of Y axis")
    parser.add_option("--ymin",dest="Ymin",default=None,help="Minimum of Y axis")


    (options,args) = parser.parse_args(ArgList)

    if ( len(args)!= 1):
        parser.error("incorrect number of arguments.")

    #aKPT,list_band = f_ReadBandFile(options.stBandFileName,options.stKptFileName)
    list_band = f_Band_ReadFromFile(options.BandFileName)
    aKPT = KPointsT(options.KptFileName)

    if (options.Band2FileName != None):
        list_band2 = f_Band_ReadFromFile(options.Band2FileName)
    else:
        list_band2 = None

    eig_fermi = 0.0
    if ( options.FermiFileName != None):
        eig_fermi = f_ReadFileFloatValue(options.FermiFileName) 
        if ( eig_fermi == None ):
            print("Warning: fermi file not found, use 0.0 instead")
            eig_fermi = 0.0
    eig_fermi2 = 0.0
    if ( options.Fermi2FileName != None):
        eig_fermi2 = f_ReadFileFloatValue(options.Fermi2FileName) 
        if ( eig_fermi == None ):
            print("Warning: second fermi file not found, use 0.0 instead")
            eig_fermi2 = 0.0
    
    eig_max = None if options.Ymax == None else float(options.Ymax)
    eig_min = None if options.Ymin == None else float(options.Ymin)

    f_PlotBand(aKPT,list_band,eig_fermi,bAlignData=options.Align,list_band2=list_band2,eig_fermi2=eig_fermi2,y_min=eig_min,y_max=eig_max)

if __name__ == "__main__":
    Main(sys.argv)    

