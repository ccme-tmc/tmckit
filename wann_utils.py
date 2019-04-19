#!/usr/bin/env python
#post-processing for wannier 90
from constants import Ang2Bohr
from common_caseutil import KPointsT,f_Data_ReadTwoCol
from band_utils import BandsT
from io_utils import io_grep_lines

def wann_read_length_unit(filename):
    '''
    Read the length unit used in wannier90 from stdout
    '''
    return io_grep_lines(filename,"Length Unit",-1,5)

def wann_read_band(filename_band,filename_kpt=None,filename_gnu=None,filename_out=None):
    '''
    Read band structure information from wannier90 output
    Indeed wannier90 output can be directly plotted but sometimes we need to combine two picture!
    :param filename_band: the name of eigenvalues ( normally seedname_band.dat )
    :param filename_kpt: the name of eigenvalues ( normally seedname_band.kpt;if not provided treated same as filename_band )
    :param filename_gnu: the name of gnuplot script ( if not provided treated same as band filename)
    :param filename_out: the name of stdout ( if not provided deduced from band filename)
    '''
    seedname = filename_band[:-9]
    if (filename_kpt == None):
        filename_kpt = seedname + "_band.kpt"
    if (filename_gnu == None):
        filename_gnu = seedname + "_band.gnu"
    if (filename_out == None):
        filename_out = seedname + ".wout"

    unit = wann_read_length_unit(filename_out)

#K-points
    kpt1 = KPointsT()
    f = open(filename_kpt)
    n = int(f.readline())
    list_k = []
    for i in xrange(n):
#Weight is ignored
        list_k.append([float(x) for x in f.readline().split()][0:3])
    f.close()
    kpt1.ReadFromList(list_k)
#Bands
    list_band = f_Data_ReadTwoCol(filename_band)
    list_kcoord = [float(x[0]) for x in list_band]
    list_band = [ [float(y) for y in x[1:]] for x in list_band]

#Special k-points from gnuplot
    f = open(filename_gnu)
    for line in f:
        if ("xtics" in line):
            line2 = line[line.index("(")+1:line.index(")")]
            ar = line2.split(",")
            for spec in ar:
                kcoord_last = 0.0
                kcoord_now = 0.0
                spec_name = spec[spec.index('"')+1:spec.rindex('"')].strip()
                spec_value = float(spec[spec.rindex('"')+1:])
#               print(spec_name,spec_value)
#Search nearest point
                for i,value in enumerate(list_kcoord):
#                   print(value)
                    if (spec_value <= value): #
                        if (abs(spec_value-value) <= abs(kcoord_last-spec_value)):
                            kpt1.listKPt[i][0] = spec_name
                        else:
                            kpt1.listKPt[i-1][0] = spec_name
                        break
                    kcoord_last = value
                    #Detect the last one if not used
                    if (i == len(list_kcoord)-1):
                        kpt1.listKPt[-1][0] = spec_name
#               print("Set kpt%i to %s" % ( i,spec_name))
            break
#Convert unit to Bohr
#However list_kcoord is not used later, maybe useful later
    if (unit == "Ang"):
        print("Convert unit from Angstrom to Bohr")
        list_kcoord = [x*Ang2Bohr for x in list_kcoord]

    return BandsT(kpt1,list_band,None,None,None,1)
#    return kpt1,list_band,None,None,1


def wann_read_hr(filename):
    '''
    Read the Hamiltonian from hr file
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

    num_wann = int(lines[1])
    num_ws = int(lines[2])
    list_ws = []
    n_line_ws = (num_ws - 1) // 15 + 1
    for i in xrange(0, n_line_ws):
        list_ws += [int(x) for x in lines[3+i].split()]

    list_R = []
    list_h = []
    ix_line = n_line_ws + 3
    for ix in xrange(0, num_ws):
        list_h_ws = []
        for i in xrange(num_wann):
            list_h_ws_i = []
            for j in xrange(num_wann):
                ar = lines[ix_line].split()
                if (i == 0 and j == 0):
                    list_R.append(ar[:3])
                ar = ar[-2:]
                ix_line += 1
                list_h_ws_i.append(float(ar[0]) + float(ar[1]) * 1j)
            list_h_ws.append(list_h_ws_i)
        list_h.append(list_h_ws)

    return num_wann, list_ws, list_R, list_h

def wann_read_centers(filename):
    '''
    Read centers.dat

    :return two list, one for Wannier function centres and another for atom positions
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()

    num_center = int(lines[0])
    l1 = [line.split() for line in lines[2:2+num_center]]
    l1 = [[x[0], float(x[1]), float(x[2]), float(x[3])] for x in l1]

    l2 = [line.split() for line in lines[2+num_center:]]
    l2 = [[x[0], float(x[1]), float(x[2]), float(x[3])] for x in l2]

    if (len(l2) == 0):#Atoms are included in l1, I don't know when will this occur
#Finished with atom name X
        for i, line in enumerate(l1):
            if (line[0] != "X"):
                break
        l2 = l1[i:]
        l1 = l1[:i]


    return l1, l2

def wann_read_latt(filename):
    '''
    Read unit_cell_cart from wannier90 input file
    Unit is converted to Bohr
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if ("begin unit_cell_cart" in line.lower()):
            ix1 = i+1
            break

#Unit line
    st_unit = lines[ix1].strip().lower()
    if (not st_unit.isalpha()):
        st_unit  = "ang"
        ix1 -= 1

    coef = Ang2Bohr if st_unit == "ang" else 1

    mat = []
    for ix in xrange(ix1+1, ix1+4):
        mat.append([float(x) * coef for x in lines[ix].split()])
    
    return mat

def wann_read_atom_pos(filename):
    '''
    Read atoms_cart from wannier90 input file
    Unit converted to Bohr
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if ("begin atoms_cart" in line.lower()):
            ix1 = i+1
            break

#Unit line
    st_unit = lines[ix1].strip().lower()
    if (not st_unit.isalpha()):
        st_unit  = "ang"
        ix1 -= 1

    coef = Ang2Bohr if st_unit == "ang" else 1

    list_atom = []
    ix = ix1 + 1
    while ("end" not in lines[ix]):
        list_atom.append([float(x) * coef for x in lines[ix].split()[1:]])
        ix += 1

    return list_atom

