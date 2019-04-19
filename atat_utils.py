#!/usr/bin/env python
#Tools for ATAT
from constants import Ang2Bohr
from common_caseutil import Lattice
from list_utils import f_Matrix_dot, f_Matrix_transpose

def atat_read_struct(filename):
    '''
    Read structure from ATAT input file
    Must be 3-line conventional cell + 3-line primitive cell
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

#Convert to Bohr
    mat_conv = [[ float(x) * Ang2Bohr for x in line.split()] for line in lines[:3]]
#Convert mat_prim from conv unit to cart unit
    mat_prim = [[ float(x) for x in line.split()] for line in lines[3:6]]
    mat_prim = f_Matrix_transpose(f_Matrix_dot(f_Matrix_transpose(mat_conv), f_Matrix_transpose(mat_prim)))

    latt = Lattice()
    latt.ReadFromRawCellVector(mat_conv, mat_prim)
    for line in lines[6:]:
        ar = line.split()
        if ("," in ar[-1]):
            ar[-1] = ar[-1].split(",")[0]
        latt.AddAtom([ [ar[-1]] + [float(x) for x in ar[:3]] ], unit="crystal", latt="conv")

    return latt

