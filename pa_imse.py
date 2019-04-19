#!/usr/bin/env python
from __future__ import print_function
import os,sys
from argparse import ArgumentParser
from constants import Ha2eV
from py_pade import pade_multi,ac_split
from common_caseutil import f_Data_WriteMultiCol,f_Split_RangeString

parser = ArgumentParser(description="Read the data on the imaginary axis in ABINIT AC calculations and extrapolate to the real axis. Only the first set it used.")
parser.add_argument('filename',type=str,help="The ABINIT stdout")
parser.add_argument('dirname',type=str,help="The output directory")
parser.add_argument('-r',dest="energyrange",type=str,default="-20~20~0.1",help="The energy window on the real axis to extrapolate in the format of min~max~step, unit eV")

args = parser.parse_args()

f = open(args.filename)
lines = f.readlines()
f.close()

if (not os.path.exists(args.dirname)):
    os.mkdir(args.dirname)

freq_real = f_Split_RangeString(args.energyrange)

i = 0
while i<len(lines):
    line = lines[i]
    if ("Sigmac for pade" in line):
        nk = int(line[21:25])
        nb = int(line[31:35])
        list_se = []
        i += 1
        while (lines[i][0] != 'U'):
            list_se.append([float(x)*Ha2eV for x in lines[i].split()])
            i += 1
#       se_real,z_real = pade_multi( [x[0]+1j*x[1] for x in list_se],[x[2]+1j*x[3] for x in list_se], freq_real)
        se_real,z_real = ac_split([x[0]+1j*x[1] for x in list_se],[x[2]+1j*x[3] for x in list_se], freq_real,pade_multi)
        filename1 = os.path.join(args.dirname,"%i-%i-imag.dat" % (nk,nb))
        filename2 = os.path.join(args.dirname,"%i-%i-real.dat" % (nk,nb))
        if (os.path.exists(filename1)):
            print("Reached next dataset, stop.")
            break
        f_Data_WriteMultiCol( [x[1:] for x in list_se],filename1)
        data2 = list(zip(freq_real,[x.real for x in se_real], [x.imag for x in se_real]))
        f_Data_WriteMultiCol(data2,filename2)
    i += 1

