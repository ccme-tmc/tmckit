#!/usr/bin/env python
#This module read yambo output
from common_caseutil import Lattice, KPointsT
from band_utils import BandsT
from qesp_utils import qesp_read_band

def yambo_read_band(filename, filename_pw_xml=None):
    '''
    Read band structure from yambo o-*.qp file
    Additional data-file.xml can be used to indicate information
    Lattice is a fake and not used
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

    ar = [[float(x) for x in line.split()] for line in lines if line[0] != "#"]
    for a1 in ar:
        a1[0] = int(a1[0])
        a1[1] = int(a1[1])
    ix_k = [x[0] for x in ar]
    ix_b = [x[1] for x in ar]
    ix_k_min = min(ix_k)
    ix_k_max = max(ix_k)
    ix_b_min = min(ix_b)
    ix_b_max = max(ix_b)

#If spin=2, there is 6th column for spin
#If not, only 5 columns
    if (len(ar[0]) == 5):
        num_spin = 1
    elif (len(ar[0]) == 6):
        num_spin = 2
    else:
        raise ValueError("Unknown number of columns %i" % len(ar[0]))

    if (filename_pw_xml is not None):
        band = qesp_read_band(filename_pw_xml)
        if (band.num_spin != num_spin):
            raise ValueError("Inconsistent spin for QE file and Yambo file")
        #Modify bands
#Set spin
        if (num_spin == 1):
            for a1 in ar:
                a1.append(0)
        else:
            for a1 in ar:
                a1[-1] = int((1-a1[-1])/2)
#Set GW eigenvalues
#Note yambo make a shift to input DFT bands to make KS VBM=0, 
#We shuold get it back
#Method 1 : calc shift and apply it to Yambo data
#       a = ar[0]
#       shift = a[2] - band.prop[a[-1]].eig[a[0]-1][a[1]-1]

#       for a in ar:
#           band.prop[a[-1]].eig[a[0]-1][a[1]-1] = a[2] + a[3] - shift
        
#Method 2 : add the correction to DFT instead of reading yambo DFT
        for a in ar:
            band.prop[a[-1]].eig[a[0]-1][a[1]-1] += a[3]

#Calculate properties
        if (band.vbm is not None):
#Clear and recalculate
            band.vbm = None
            band.guess_vbm()

    else:
        if (num_spin == 1):
            list_eig = [[None] * (ix_b_max - ix_b_min + 1) for x in range(ix_k_max-ix_k_min+1)]
            for x in ar: #Eo + (E-Eo)
                list_eig[x[0] - ix_k_min][x[1] - ix_b_min] = x[2] + x[3]
        elif (num_spin == 2):
            num_spin = 2
            list_eig = [
                    [[None] * (ix_b_max - ix_b_min + 1) for x in range(ix_k_max-ix_k_min+1)],
                    [[None] * (ix_b_max - ix_b_min + 1) for x in range(ix_k_max-ix_k_min+1)]]
            for x in ar: #Eo + (E-Eo)
                list_eig[int((1-x[5])/2)][x[0] - ix_k_min][x[1] - ix_b_min] = x[2] + x[3]

#Dummy kpt
        latt = Lattice()
        latt.ReadFromRawCellVector([[1,0,0],[0,1,0],[0,0,1]])
        kpt = KPointsT()
        kpt.ReadFromList([[x, 0, 0] for x in range(len(list_eig))])
        kpt.unit = "crystal"
        kpt.latt = latt

        band = BandsT(kpt, list_eig, num_spin = num_spin)
        
    return band


