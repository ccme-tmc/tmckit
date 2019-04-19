#!/usr/bin/env python
'''@package bdf_utils
Contains classes and functions for BDF
'''
from band_utils import BandsT
from common_caseutil import KPointsT, f_GetReciprocalLattice

def bdf_read_band(filename):
    '''
    Read a band structure from bdf output
    '''
    f = open(filename)
#Read reciporocal vectors
    lines = f.readlines()
    f.close()
    mK = [ [float(y) for y in x.split()[1:4]] for x in lines[2:5] ]
    eig_fermi = float(lines[6].split()[-1])
    print(eig_fermi)
    list_kpt = []
    list_eig = []
    for line in lines[7:]:
        ar = line.split()
        #"A" is marked as non-special k-points
        if (ar[1] == "A"):
            ar[1] = ""
        list_kpt.append(ar[1:2] + [float(x) for x in ar[2:5]]+[1.0])
        list_eig.append([float(x) for x in ar[5:]])
    kpt = KPointsT()
    kpt.stMode = "crystal"
    kpt.ReadFromList(list_kpt)
    kpt.ConvertUnit("cart",mLatt=f_GetReciprocalLattice(mK))
    band = BandsT(kpt,list_eig,eig_fermi,eig_fermi)
    band.b_metal = True #Just set to plot fermi
    return band

