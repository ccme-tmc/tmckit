#!/usr/bin/env python
from math import *
elements = ['H' , 'He', 'Li', 'Be', 'B',  'C',  'N',  'O' , 'F' , 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P' , 'S',  'Cl', 'Ar', 'K' , 'Ca',\
            'Sc', 'Ti', 'V' , 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',\
            'Rb', 'Sr', 'Y' , 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I' , 'Xe',\
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',\
                  'Hf', 'Ta', 'W' , 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',\
            'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U' , 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',\
                  'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mg', 'Ds', 'Rg']
elements_info = [ \
     ('H ',   1.00),\
     ('He',   4.00),\
     ('Li',   6.94),\
     ('Be',   9.01),\
     ('B ',  10.81),\
     ('C ',  12.01),\
     ('N ',  14.01),\
     ('O ',  16.00),\
     ('F ',  19.00),\
     ('Ne',  20.18),\
     ('Na',  22.99),\
     ('Mg',  24.31),\
     ('Al',  26.98),\
     ('Si',  28.09),\
     ('P ',  30.97),\
     ('S ',  32.07),\
     ('Cl',  35.45),\
     ('Ar',  39.95),\
     ('K ',  39.09),\
     ('Ca',  40.08),\
     ('Sc',  44.96),\
     ('Ti',  47.87),\
     ('V ',  50.94),\
     ('Cr',  52.00),\
     ('Mn',  54.94),\
     ('Fe',  55.85),\
     ('Co',  58.93),\
     ('Ni',  58.69),\
     ('Cu',  63.55),\
     ('Zn',  65.41),\
     ('Ga',  69.72),\
     ('Ge',  72.64),\
     ('As',  74.92),\
     ('Se',  78.96),\
     ('Br',  79.90),\
     ('Kr',  83.80),\
     ('Rb',  85.47),\
     ('Sr',  87.62),\
     ('Y ',  88.91),\
     ('Zr',  91.22),\
     ('Nb',  92.91),\
     ('Mo',  95.94),\
     ('Tc',  97.91),\
     ('Ru', 101.07),\
     ('Rh', 102.91),\
     ('Pd', 106.42),\
     ('Ag', 107.87),\
     ('Cd', 112.41),\
     ('In', 114.82),\
     ('Sn', 118.71),\
     ('Sb', 121.76),\
     ('Te', 127.60),\
     ('I ', 126.90),\
     ('Xe', 131.29),\
     ('Cs', 132.91),\
     ('Ba', 137.33),\
     ('La', 138.91),\
     ('Ce', 140.12),\
     ('Pr', 140.91),\
     ('Nd', 144.24),\
     ('Pm', 144.91),\
     ('Sm', 150.36),\
     ('Eu', 151.96),\
     ('Gd', 157.25),\
     ('Tb', 158.93),\
     ('Dy', 162.50),\
     ('Ho', 164.93),\
     ('Er', 167.26),\
     ('Tm', 168.93),\
     ('Yb', 173.04),\
     ('Lu', 174.97),\
     ('Hf', 178.49),\
     ('Ta', 180.95),\
     ('W ', 183.84),\
     ('Re', 186.21),\
     ('Os', 190.23),\
     ('Ir', 192.22),\
     ('Pt', 195.08),\
     ('Au', 196.97),\
     ('Hg', 200.59),\
     ('Tl', 204.38),\
     ('Pb', 207.20),\
     ('Bi', 208.98),\
     ('Po', 208.98),\
     ('At', 209.99),\
     ('Rn', 222.02),\
     ('Fr', 223.02),\
     ('Ra', 226.03),\
     ('Ac', 227.03),\
     ('Th', 232.04),\
     ('Pa', 231.04),\
     ('U ', 238.03),\
     ('Np', 237.05),\
     ('Pu', 244.06),\
     ('Am', 243.06),\
     ('Cm', 247.07),\
     ('Bk', 247.07),\
     ('Cf', 251.08),\
     ('Es', 252.08),\
     ('Fm', 257.10),\
     ('Md', 258.10),\
     ('No', 259.10),\
     ('Lr', 260.11),\
     ('Rf', 261.11),\
     ('Db', 262.11),\
     ('Sg', 263.12),\
     ('Bh', 264.12),\
     ('Hs', 265.13),\
     ('Mt', 266.13),\
     ('Ds', 269.00),
     ('Rg', 272.00)]  


def f_Element_Symbol_to_Z(symbl):
  try: 
    z = elements.index(symbl.strip()) + 1
  except:
    print "WARNING in f_Element_Symbol_to_Z: %s is not an element symbol" %(symbl)
    z = 0 
  return z

def f_Element_Symbol_to_Mass(symbl):
  z = f_Element_Symbol_to_Z(symbl) 
  return elements_info[z-1][1]

def f_Element_Z_to_Symbol(znucl):
  return elements[znucl-1]


element_configuration=[
    [],
    ['',['1s',1]],
    ['',['1s',2]],
    ['He',['2s',1]],
    ['He',['2s',2]],
    ['He',['2s',2],['2p',1]],
    ['He',['2s',2],['2p',2]],
    ['He',['2s',2],['2p',3]],
    ['He',['2s',2],['2p',4]],
    ['He',['2s',2],['2p',5]],
    ['He',['2s',2],['2p',6]],
    ['Ne',['3s',1]],
    ['Ne',['3s',2]],
    ['Ne',['3s',2],['3p',1]],
    ['Ne',['3s',2],['3p',2]],
    ['Ne',['3s',2],['3p',3]],
    ['Ne',['3s',2],['3p',4]],
    ['Ne',['3s',2],['3p',5]],
    ['Ne',['3s',2],['3p',6]],
    ['Ar',['4s',1]],
    ['Ar',['4s',2]],
    ['Ar',['3d',1],['4s',2]],
    ['Ar',['3d',2],['4s',2]],
    ['Ar',['3d',3],['4s',2]],
    ['Ar',['3d',5],['4s',1]],
    ['Ar',['3d',5],['4s',2]],
    ['Ar',['3d',6],['4s',2]],
    ['Ar',['3d',7],['4s',2]],
    ['Ar',['3d',8],['4s',2]],
    ['Ar',['3d',10],['4s',1]],
    ['Ar',['3d',10],['4s',2]],
    ['Ar',['3d',10],['4s',2],['4p',1]],
    ['Ar',['3d',10],['4s',2],['4p',2]],
    ['Ar',['3d',10],['4s',2],['4p',3]],
    ['Ar',['3d',10],['4s',2],['4p',4]],
    ['Ar',['3d',10],['4s',2],['4p',5]],
    ['Ar',['3d',10],['4s',2],['4p',6]],
    ['Kr',['5s',1]],
    ['Kr',['5s',2]],
    ['Kr',['4d',1],['5s',2]],
    ['Kr',['4d',2],['5s',2]],
    ['Kr',['4d',4],['5s',1]],
    ['Kr',['4d',5],['5s',1]],
    ['Kr',['4d',5],['5s',2]],
    ['Kr',['4d',7],['5s',1]],
    ['Kr',['4d',8],['5s',1]],
    ['Kr',['4d',10],['5s',0]], # 5s0 is for valence auto-detect
    ['Kr',['4d',10],['5s',1]],
    ['Kr',['4d',10],['5s',2]],
    ['Kr',['4d',10],['5s',2],['5p',1]],
    ['Kr',['4d',10],['5s',2],['5p',2]],
    ['Kr',['4d',10],['5s',2],['5p',3]],
    ['Kr',['4d',10],['5s',2],['5p',4]],
    ['Kr',['4d',10],['5s',2],['5p',5]],
    ['Kr',['4d',10],['5s',2],['5p',6]],
    ['Xe',['6s',1]],
    ['Xe',['6s',2]],
    ['Xe',['5d',1],['6s',2]],
    ['Xe',['4f',1],['5d',1],['6s',2]],
    ['Xe',['4f',3],['6s',2]],
    ['Xe',['4f',4],['6s',2]],
    ['Xe',['4f',5],['6s',2]],
    ['Xe',['4f',6],['6s',2]],
    ['Xe',['4f',7],['6s',2]],
    ['Xe',['4f',7],['5d',1],['6s',2]],
    ['Xe',['4f',9],['6s',2]],
    ['Xe',['4f',10],['6s',2]],
    ['Xe',['4f',11],['6s',2]],
    ['Xe',['4f',12],['6s',2]],
    ['Xe',['4f',13],['6s',2]],
    ['Xe',['4f',14],['6s',2]],
    ['Xe',['4f',14],['5d',1],['6s',2]],
    ['Xe',['4f',14],['5d',2],['6s',2]],
    ['Xe',['4f',14],['5d',3],['6s',2]],
    ['Xe',['4f',14],['5d',4],['6s',2]],
    ['Xe',['4f',14],['5d',5],['6s',2]],
    ['Xe',['4f',14],['5d',6],['6s',2]],
    ['Xe',['4f',14],['5d',7],['6s',2]],
    ['Xe',['4f',14],['5d',9],['6s',1]],
    ['Xe',['4f',14],['5d',10],['6s',1]],
    ['Xe',['4f',14],['5d',10],['6s',2]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',1]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',2]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',3]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',4]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',5]],
    ['Xe',['4f',14],['5d',10],['6s',2],['6p',6]],
    ['Rn',['7s',1]],
    ['Rn',['7s',2]],
    ['Rn',['6d',1],['7s',2]],
    ['Rn',['6d',2],['7s',2]],
    ['Rn',['5f',2],['6d',1],['7s',2]],
    ['Rn',['5f',3],['6d',1],['7s',2]],
    ['Rn',['5f',4],['6d',1],['7s',2]],
    ['Rn',['5f',6],['7s',2]],
    ['Rn',['5f',7],['7s',2]],
    ['Rn',['5f',7],['6d',1],['7s',2]],
    ['Rn',['5f',9],['7s',2]],
    ['Rn',['5f',10],['7s',2]],
    ['Rn',['5f',11],['7s',2]],
    ['Rn',['5f',12],['7s',2]],
    ['Rn',['5f',13],['7s',2]],
    ['Rn',['5f',14],['7s',2]],
]

def f_Element_Z_to_ValenceElectronCount(nZ):
    '''
    Return count of valence electron of specified element Z
    Definition of valence electron : s,p are always valence, d is valence if semi-occupied OR no p exists, nf is valence f is semi-occupied OR (n+1)d2 exists
    '''
    confVal = element_configuration[nZ][1:] 
    nVal = 0
    b_d = True
    b_f = True
    for  conf in confVal: #First cycle, add s/p electron, detect d/f electron
        #print(conf)
        if ( conf[0][1] == 's' ):
            nVal += conf[1]
        elif ( conf[0][1] == 'p' ):
            if ( conf[1] >= 1):
                b_d = False
            nVal += conf[1]
        elif ( conf[0][1] == 'd' ):
            if ( conf[1] >= 2):
                b_f = False
    for conf in confVal: #Second cycle, add d/f electron
        if ( conf[0][1] == 'd'):
            if ( b_d or conf[1] < 10 ):
                nVal += conf[1]
        elif ( conf[0][1] =='f'):
            if ( b_f or conf[1] < 14):
                nVal += conf[1]

    return nVal

def f_Element_Z_to_ExtraValenceElectronCount(nZ):
    '''
    Return number of unoccupied "electrons" most chemical or physical properties related to. 
    For example, O has 2 "unoccupied" and F has 1. Ti has 8 ( d10 + s2 is what we concern )
    '''
    confVal = element_configuration[nZ][1:]
    confVal2 = element_configuration[nZ+1][1:]
    nVal = 0
    b_p = False
    b_d = False
    b_f = False
    b_p_exist = False
    b_d_exist = False
#first turn to determine valence
    for conf in confVal:
        if ( conf[0][1] == 's' ):
            nVal += 2 - conf[1]
            if ( int(conf[0][0]) >= 2):
                b_p_exist = True
            if ( int(conf[0][0]) >= 4):
                b_d_exist = True
        elif ( conf[0][1] == 'p' ):
            nVal += 6-conf[1]
            b_p = conf[1] > 0
        elif ( conf[0][1] == 'd' ):
            nVal += 10-conf[1]
            b_d =  conf[1] > 0
        elif ( conf[0][1] == 'f' ):
            nVal += 14 -conf[1]
            b_f = conf[1] > 0

#second turn to determine whether we need extra p / d orbital
# IIA/IB/IIB needs p, IIA needs d ( n>3 ), XIIIA needs no more
    if ( b_f and not b_d  ): # f-region
        nVal += 10
    elif ( nVal <= 1): #I or II
        if ( not b_p ): #Skip XIIIA
            if ( b_d ): #IB/IIB
                nVal +=  6
            else: #IIA or XIIIA
                if ( nVal == 0 and b_d_exist ): #IA does not require d
                    nVal += 10
                if ( not b_f and b_p_exist ):#Skip f-region
                    nVal += 6



    return nVal


   

        

