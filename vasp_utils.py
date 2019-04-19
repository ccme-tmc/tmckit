#!/usr/bin/env python
#VASP package utility

import copy
import commands,sys,string,math,re, subprocess,shutil 

from math import pi
from constants import *
from common_caseutil import deprecated,KPointsT,Lattice,f_file_ensuredir,RawFileIO
from band_utils import f_Band_GetVBMFromElectron,f_Band_GetGap,BandsT
from dos_utils import OrbitalNameT
from io_utils import *
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
from py_xmlio import *

class VaspInput():
    '''
    Represent input files of VASP, including INCAR, POSCAR and KPOINTS only currently
    '''
    def __init__(self,dirname=None):
        '''
        Initialize from specific directory
        '''
        self.incar =  None
        self.poscar = None

        self.kpt_comment = None #First line
        self.kpt_num = None #Second line
        self.kpt_mode = None #Third line for kpt list and fourth line for line-mode
        self.kpt_mp = None #M-P parameters, including number of points in x,y,z 
        self.kpt_shift = None #M-P shift of x,y,z
        self.kpt_mesh_vec = None # a 3x3 array as grid vector for mesh generating
        self.kpt_list = None #list of k-points in list mode
        self.kpt_line = None #list of line in line mode
        self.kpt_cut = None #Cut of k-points in auto mode (Gamma-centered MP)
        self.kpt_extra = None #other part which just saved

        if (dirname is not None):
            self.read(dirname)

    def __trim__(self,line):
        '''
        Remove comment part started with # or !
        '''
        i = line.find("#")
        if (i != -1):
            line = line[:i]
        i = line.find("!")
        if (i != -1):
            line = line[:i]
        return line
    
    
    def read(self,dirname=None):
        '''
        Read from a folder
        '''
        dirname = f_file_ensuredir(dirname)
        self.incar = RawFileIO(os.path.join(dirname,"INCAR"))
        self.poscar = RawFileIO(os.path.join(dirname,"POSCAR"))
        self.read_kpt(os.path.join(dirname,"KPOINTS"))

#ASE IO
#       self.atoms = self.ase.io.read(os.path.join(dirname,"POSCAR"),"vasp")


    def write(self,dirname=None):
        '''
        Write into a folder
        '''
        dirname = f_file_ensuredir(dirname)
        self.incar.write(os.path.join(dirname,"INCAR"))
        self.poscar.write(os.path.join(dirname,"POSCAR"))
        self.write_kpt(os.path.join(dirname,"KPOINTS"))

    def read_kpt(self,filename=None):
        '''
        Read the k-points
        '''
        if (filename is None):
            filename = "KPOINTS"
        lines = RawFileIO(filename).lines
        self.kpt_comment = lines[0].strip()
        self.kpt_num = int(lines[1].split()[0])
        self.kpt_mode = lines[2].split()[0]
        m = self.kpt_mode[0].lower()
        if (self.kpt_num != 0):
            if (m != 'l'):
                self.kpt_list = []
                i = 3
                while (len(self.kpt_list) < self.kpt_num):
                    kpt = [float(x) for x in self.__trim__(lines[i]).split()]
                    if (len(kpt) == 4):#Only add non-empty line
                        self.kpt_list.append(kpt)
                    i += 1
                self.kpt_extra = lines[i:]
            else:
                self.kpt_mode = lines[3].split()[0]
                self.kpt_line = []
                i = 4
                while (i < len(lines)):
                    kpt = [float(x) for x in self.__trim__(lines[i]).split()]
                    if (len(kpt) == 3):#Only add non-empty line
                        self.kpt_line.append(kpt)
                    i += 1
        elif (m =='a'):#Auto
            self.kpt_cut = int(self.__trim__(lines[3]).split()[0])
        elif (m == 'g' or m == 'm'): #MP
            self.kpt_mp = [int(x) for x in self.__trim__(lines[3]).split()]
            if (len(lines) > 4 and len(self.__trim__(lines[4])) > 0):
                self.kpt_shift = [float(x) for x in self.__trim__(lines[4]).split()]
            else:
                self.kpt_shift = None
        else:#Grid vector
            self.kpt_mesh_vec = [[float(x) for x in self.__trim__(y)] for y in lines[3:6]]
            self.kpt_shift = [float(x) for x in self.__trim__(lines[6]).split()]
        return
    
    def write_kpt(self,filename=None):
        '''
        Write KPOINTS file
        '''
        if (filename is None):
            filename = "KPOINTS"
        f = open(filename,'w')
        f.write("%s\n" % self.kpt_comment)
        f.write("%i\n" % self.kpt_num)
#if kpt_line is not None then use line mode
        if (self.kpt_line is not None):
            f.write("Line-mode\n")
        f.write("%s\n" % self.kpt_mode)
        m = self.kpt_mode[0].lower()

        if (self.kpt_line is not None):
            for i in xrange(len(self.kpt_line)/2):
                f.write("%f %f %f\n" % tuple(self.kpt_line[i*2]))
                f.write("%f %f %f\n\n" % tuple(self.kpt_line[i*2+1]))
        elif (self.kpt_num != 0): #List of k mode
            for kpt in self.kpt_list:
                f.write("%f %f %f %f \n" % tuple(kpt))
        else:
            if (m == 'a'):
                f.write("%i\n" % self.kpt_cut)
            elif (m == 'g' or m == 'm'):
                f.write("%i %i %i\n" % tuple(self.kpt_mp))
                if (self.kpt_shift is not None):
                    f.write("%f %f %f\n" % tuple(self.kpt_shift))
            else:
                for i in xrange(3):
                    f.write("%f %f %f\n" % tuple(self.kpt_mesh_vec[i]))
                f.write("%f %f %f\n" % tuple(self.kpt_shift))
        f.close()

        

def vasp_ReadLattVec(stFileName="OUTCAR"):
    '''
    Read real-vector and k-vector used in vasp
    '''
    f = open(stFileName)
    ar_stLine = f.readlines()
    f.close()

    for i in range(len(ar_stLine)-1,0,-1):
#Find the last appearance of rec lattice vectors
        if ( "reciprocal lattice vectors" in ar_stLine[i]):
            break

    matR = []
    matK = []
    for j in range(i+1,i+4):
        #Note: unit is angstorm in VASP, convert to Bohr
        matR.append([float(x)*Ang2Bohr for x in ar_stLine[j].split()[0:3]])
        matK.append([2*pi*float(x)/Ang2Bohr for x in ar_stLine[j].split()[3:6]])

    aCell = Lattice()
    aCell.ReadFromPrimitiveCellVector(matR)

    return matR,matK
    

def vasp_ReadBand(stFileName="EIGENVAL",filetype=None):
    '''
    Read k-point and band from vasp EIGENVAL file or OUTCAR file
    WARNING: EIGENVAL does not handle non-integer electrons! Use OUTCAR or vasprun.xml if it is that case
    Note the structure is always read from vasprun.xml

    Another problem is the occuptation numbers are 1 in vasprun.xml  and 2 in OUTCAR, if spin unpolarized

    :param stFileName: the file to read bandstructure
    :param filetype: the file type to read, can be OUTCAR or EIGENVAL. Automatically inferred from stFileName if set to None
    '''
    if (filetype is None):
        if ("EIGENVAL" in stFileName):
            filetype = "EIGENVAL"
        elif ("vasprun.xml" in stFileName):
            filetype = "vasprun.xml"
        else:
            filetype = "OUTCAR"

    f = open(stFileName)
    ar_stLine = f.readlines()
    f.close()

    list_kp = []
    list_band = []
    list_band2 = []
    list_occ = None
    bElectronFraction = False

    if (filetype == "EIGENVAL"):
#Note we cannot determine LNONCOLLINEAR in EIGENVAL, assume false
        b_noncoll = False

        num_spin = int(ar_stLine[0].split()[3])

        (nElectron,nKPt,nbnd) = [int(x) for x in ar_stLine[5].split()]

        for i in range(0,nKPt):
            i2 = i * (nbnd+2) + 7
            list_kp.append([float(x) for x in ar_stLine[i2].split()])
            ar2 =[x.split() for x in ar_stLine[i2+1:i2+nbnd+1]]
            list_band.append([float(x[1]) for x in ar2])
            if (num_spin == 2):
                list_band2.append([float(x[2]) for x in ar2])
        if (num_spin == 2):
            list_band = [list_band,list_band2]
    elif (filetype == "OUTCAR"):

        def read_tag(tag, i0):
            val = None
            for i in xrange(i0, len(ar_stLine)):
                line = ar_stLine[i]
                if (tag in line):
                    val = line.split()[2]
                    break
            if (val is None):
                raise ValueError("Cannot find tag %s" % tag)

            return val, i


        list_occ = []
        list_occ2 = []

        i = 0
        num_spin, i = read_tag("ISPIN", i)
        num_spin = int(num_spin)
        b_noncoll, i = read_tag("LNONCOLLINEAR", i)
        b_noncoll = b_noncoll == "T"
        nElectron, i = read_tag("NELECT", i)
        dElectron = float(nElectron)
        nElectron = int(dElectron)
        bElectronFraction = abs(dElectron != nElectron) > 1e-7
#       if (bElectronFraction):
#           print("The number of electrons is not an integer.")

        efermi, i = read_tag("E-fermi", i)
        efermi = float(efermi)
        i += 2
        if (num_spin != 1):
            i += 2
        list_band_now = list_band
        list_occ_now = list_occ
        while (len(ar_stLine[i+1])>4):
            i += 1
            if (list_band_now == list_band):
                list_kp.append([float(x) for x in ar_stLine[i].split()[-3:]])
            i += 2
            ar2 = []
            ar_occ = []
            line = ar_stLine[i]
            while (len(line) > 4):
                ar_line = line.split()
                ar2.append(float(ar_line[1]))
                ar_occ.append(float(ar_line[2]))
                i += 1
                line = ar_stLine[i]
            list_band_now.append(ar2)
            list_occ_now.append(ar_occ)
            if (num_spin != 1 and "spin" in ar_stLine[i+1]):
                list_band_now = list_band2
                list_occ_now = list_occ2
                i += 2

        nKPt = len(list_band)
        nbnd = len(list_band[0])
#       print(nElectron, nKPt, nbnd)
        if (num_spin == 2):
            list_band = [list_band, list_band2]
            list_occ = [list_occ, list_occ2]
    elif (filetype == "vasprun.xml"):
#Read from vasprun.xml
        for event, elem in ET.iterparse(stFileName,events=('end',)):
            if (elem.tag == "parameters"):
                node_parameters = elem
            elif (elem.tag == "eigenvalues"):
                node_eigenvalues = elem
            elif (elem.tag == "kpoints"):
                node_kpoints = elem
            elif (elem.tag == "dos"):
                node_dos = elem 
                break

#Read LNONCOLLINEAR
        b_noncoll = node_parameters.find('./separator[@name="electronic"]/separator[@name="electronic spin"]/i[@name="LNONCOLLINEAR"]').text
        b_noncoll = b_noncoll.strip() == "T"

#Read NELECT
        
        dElectron = float(node_parameters.find('./separator[@name="electronic"]/i[@name="NELECT"]').text)
        nElectron = int(dElectron)
        bElectronFraction = abs(dElectron != nElectron) > 1e-7

#Read band
        node_band = node_eigenvalues
        list_band = []
        list_occ = []
        list_occ2 = []
#3-level nested , hard to read so split them
        for node1 in node_band.findall("./array/set/set"): #Spin
            band_spin = []
            occ_spin = []
            for node2 in node1.findall("./set"): #Kpt
                data_kpt = [node3.text.split() for node3 in node2.findall("./r")]
                band_spin.append([float(x[0]) for x in data_kpt])
                occ_spin.append([float(x[1]) for x in data_kpt])
            list_band.append(band_spin)
            list_occ.append(occ_spin)
#       list_band = reduce(lambda x,y: x+y, list_band)
#       list_occ = reduce(lambda x,y: x+y, list_occ)
#Read kpoints coordinates
        list_kp = [[float(x) for x in nodek.text.split()] for nodek in node_kpoints.findall('./varray[@name="kpointlist"]/v')]

#Read other properties
        num_spin = int(node_parameters.find('./separator[@name="electronic"]/separator[@name="electronic spin"]/i[@name="ISPIN"]').text)
        efermi = float(node_dos.find('./i[@name="efermi"]').text)
#Do not use array for each spin if only 1 spin
        if (num_spin == 1):
            list_band = list_band[0]
            list_occ = list_occ[0]

    aKPT = KPointsT()
    aKPT.ReadFromList(list_kp)
    aKPT.stMode = "crystal"
    
    #OUTCAR is not necasseary, it is used to convert k-point unit
    try:
        latt = vasp_read_latt("initialpos","vasprun.xml")
    except ET.ParseError:
        matR,matK = vasp_ReadLattVec("OUTCAR")
        latt = Lattice()
        latt.ReadFromRawCellVector(matR)

    aKPT.latt = latt

#   aKPT.ConvertUnit("cart",matR)
    aKPT.ConvertUnit("cart")

#Each band is assumed to be double occupied in BandsT
#So we use 2 times electron in NONCOLLINEAR where only single occupied
    if (b_noncoll):
        print("Noncollinear calculation, band.electron = 2*actual electrons")
        nElectron *= 2 

    band = BandsT(aKPT,list_band,None,None,nElectron,num_spin, list_occ=list_occ)

    #Calculate VBM
#Only do if the number of electrons is an integer
    if (bElectronFraction):
        band.fermi = efermi
    else:
        if (not band.guess_vbm()):
            band.fermi = vasp_getout("efer")
    

#   return aKPT,list_band,dVBM,nElectron,1
    return band

def vasp_read_xml_atom_species(node):
    '''
    Read atom species from "atominfo" node in vasprun.xml
    If more than two species come with the same name, then they will be appended with suffix 1,2,3...

    :return: a list of atom species names
    '''
    atominfo = node.findall("./array[@name='atoms']/set")[0]
    list_atom = [int(x.text.strip()) for x in atominfo.findall("./rc/c[2]")]
#Test if there are repeated atoms
    list_atom_type = []
    dic_atom_type = {}
    pp = node.findall("./array[@name='atomtypes']/set/rc/c[2]")
    for node1 in pp:
        name = node1.text.strip()
        if (name in dic_atom_type):
            dic_atom_type[name] += 1
        else:
            dic_atom_type[name] = 1
        list_atom_type.append((name, dic_atom_type[name]))
#Convert name/index to name+suffix; if index are always 1, no suffix
    list_atom_type = map(lambda x:x[0] if dic_atom_type[x[0]] == 1 else ("%s%i" % x), list_atom_type)
    list_atom = [list_atom_type[x-1] for x in list_atom]
    return list_atom
    


def vasp_ReadPDOS(stFileName="vasprun.xml"):
    '''
    Read projected dos from vasprun.xml
    Please ensure the projected dos is calculated.
    '''
    dos2 = None
    for event, elem in ET.iterparse(stFileName,events=('end',)):
        if (elem.tag == "parameters"):
            node_parameters = elem
        elif (elem.tag == "atominfo"):
            node_atominfo = elem
        elif (elem.tag == "dos"):
            dos2 = elem 
            break

    if (dos2 is None):
        raise ValueError("This VASP result does not contains PDOS")

    dFermi = float(dos2.findall("./i[@name='efermi']")[0].text)

#Read Total DOS

    list_total = []
    list_total_comment = []
#    for set in total.getElementsByTagName("set"):
    for set1 in dos2.findall("./total/array/set/set"):
        list_one = []
        #for r in set.getElementsByTagName("r"):
        for r in set1.findall("./r"):
            list_one.append([float(x) for x in r.text.split()])
#Transpose
        list_total.append(zip(*list_one))
        list_total_comment.append(set1.attrib["comment"])

#Read partial DOS

    list_part = []

    #Read PDOS orbital name    
    list_part_name = []
    for part in dos2.findall("./partial/array/field"):
        list_part_name.append(part.text.strip())
    #Read PDOS Data
    for set in dos2.findall("./partial/array/set/set"):
        list_set2 = set.findall("./set") #Data belongs to one atom
        for set2 in list_set2:
            list_one = []
            for r in set2.findall("./r"):
                list_one.append( [float(x) for x in r.text.split()])
            list_part.append([str(set.attrib["comment"]),str(set2.attrib["comment"]),zip(*(list_one))])

#Read atom name
    list_atom = vasp_read_xml_atom_species(node_atominfo)

#Read spin
    nSpin = int(node_parameters.findall("./separator/separator/i[@name='ISPIN']")[0].text)

#Format 
    listName = []
#Energy from Total DOS
    listEnergy = list_total[0][0]  

    listLM = [ [0,0],[1,-1],[1,0],[1,1],[2,-2],[2,-1],[2,0],[2,1],[2,2]]


    def get_dosutil_spin_from_vasp(comment):
        '''
        Read vasp spin 1/2 to dos_util spin 1/-1
        '''
        spin = int(comment.split()[-1])
     #Spin in vasp is marked as 1 and 2, convert to 1 and -1
        spin = -2 * spin + 3
        return spin

    listDOS = []
#Add total dos
    for i, (total_comment, energy_dos) in enumerate(zip(list_total_comment, list_total)):
        name = OrbitalNameT()
        name.total = True
#One or two spin?
        if ( nSpin == 2):#Only set spin for spin-polarized calculation
            name.spin = get_dosutil_spin_from_vasp(total_comment)

        listName.append(name)
        listDOS.append(energy_dos[1])
    for set1 in list_part:
        #Create name
        name = OrbitalNameT()
        name.i = int(set1[0].split()[-1])
        if ( nSpin == 2):#Only set spin for spin-polarized calculation
            name.spin = get_dosutil_spin_from_vasp(set1[1])
        name.n = ""# no n
        name.z = ""# no z
        name.species = list_atom[name.i-1] #index start from 1
        #We use custom m here
#        print list_part_name
        for i,st_name in enumerate(list_part_name):
            name2 = copy.copy(name)
            if ( st_name == "energy"): #Skip energy
                continue
#            print st_name[0]
            name2.l = OrbitalNameT.dic_l_rev[st_name[0]]
            name2.m_name = st_name[1:]
            listName.append(name2)

            listDOS.append(set1[2][i])

#       for i in xrange(1,9+1): #Read s -> d
#           name2 = copy.copy(name)
#           name2.l = listLM[i-1][0]
#           name2.m = listLM[i-1][1]
#           listName.append(name2)

#           listDOS.append(set1[2][i])

    return listName,listEnergy,listDOS,dFermi

def vasp_read_band_character(filename="vasprun.xml"):
    '''
    Read orbitals and components from vasprun.xml
    '''
    proj = None
    for event, elem in ET.iterparse(filename,events=('end',)):
        if (elem.tag == "parameters"):
            node_parameters = elem
        elif (elem.tag == "atominfo"):
            node_atominfo = elem
        elif (elem.tag == "projected"):
            proj = elem 
            break

#Read atom name
#Read atom name
    list_atom = vasp_read_xml_atom_species(node_atominfo)

#Read spin
    nSpin = int(node_parameters.findall("./separator/separator/i[@name='ISPIN']")[0].text)

#Read orbital name (atom information not included)
    list_orb = []
    for af in proj.findall("./array/field"):
        orb_name_vasp = af.text.strip()
        #Create name
        name = OrbitalNameT()
        name.i = 0 #Set atom index later
        #name.spin set later
        name.n = ""# no n
        name.z = ""# no z
        #name.species set later
        #We use custom m here
        name.l = OrbitalNameT.dic_l_rev[orb_name_vasp[0]]
        name.m_name = orb_name_vasp
        list_orb.append(name)
#Add atoms
    list_orb_atom = []
    if (nSpin == 2):
        ar_spin = [1, -1]
    else:
        ar_spin = [0]

    for spin in ar_spin:
        ar = []
        for ix, atom in enumerate(list_atom):
            for orb in list_orb:
                orb2 = copy.deepcopy(orb)
                orb2.spin = spin
                orb2.i = ix
                orb2.species = atom
                ar.append(orb2)
        list_orb_atom.append(ar)

    n_orb_one_spin = len(list_orb_atom[0])

        
#Read components
    list_character = []
    for set_spin in proj.findall("./array/set/set"):
        spin = int(set_spin.get("comment")[-1])
        list_kpt = []
        for set_kpt in set_spin.findall("./set"):
            list_band = []
            for set_band in set_kpt.findall("./set"):
                ar = []
                for r in set_band.findall("./r"):
                    ar += [float(x) for x in r.text.split()]
                list_band.append(ar)
            list_kpt.append(list_band)
        list_character.append(list_kpt)
                
    return list_orb_atom, list_character

def vasp_read_atom_species(filename="POTCAR"):
    '''
    Read atom types from POTCAR
    '''
    f = open(filename)
    list_species = []
    for line in f:
        if ("VRHFIN" in line):
            st = line.split(":")[-1].split(":")[0]
            list_species.append(st)
    f.close()

def vasp_read_latt2(name_struct="initialpos",filename="vasprun.xml"):
    '''
    Read lattice information from vasprum.xml 
    Note the unit is angstrom, convert to Bohr here
    This reads the whole xml and is very slow if it is large. Use `func:vasrp_read_latt` instead.

    :param name_struct: the structure name, default "initialpos" for POSCAR, or "finalpos" for CONTCAR
    :param filename: the xml filename, default "vasprun.xml"
    '''
    latt = Lattice()
    f = open(filename)
    line_all = f.read()
    f.close()
    line_all = line_all.replace("\0","")
    doc2 = ET.fromstring(line_all)
    doc_struct = doc2.find('./structure[@name="%s"]' % name_struct)
    node_vec = doc_struct.find('./crystal/varray[@name="basis"]')
    m_vec = []
    for v1 in node_vec.findall("./v"):
        m_vec.append([float(x)*Ang2Bohr for x in v1.text.split()])

    latt.ReadFromRawCellVector(m_vec)
#Atoms
    list_species = doc2.findall('./atominfo/array[@name="atoms"]/set/rc') 
    list_pos = doc_struct.findall('./varray[@name="positions"]/v')

    list_atoms = [ [species.findall('./c')[0].text.strip()] + 
                [float(x) for x in pos.text.split()]
            for species,pos in zip(list_species,list_pos)]
    latt.AddAtom(list_atoms,unit="crystal",latt="prim")

    return latt

def vasp_read_latt(name_struct="initialpos",filename="vasprun.xml"):
    '''
    Read lattice information from vasprum.xml 
    Note the unit is angstrom, convert to Bohr here
    This reads the whole xml and is very slow if it is large

    :param name_struct: the structure name, default "initialpos" for POSCAR, or "finalpos" for CONTCAR
    :param filename: the xml filename, default "vasprun.xml"
    '''
#Just read these two, structure must be after atominfo
    for event, elem in ET.iterparse(filename,events=('end',)):
        if (elem.tag == "atominfo"):
            doc_atom = elem
        if (elem.tag == "structure"):
            if (elem.get("name") == name_struct):
                doc_struct = elem
                break

    latt = Lattice()
    node_vec = doc_struct.find('./crystal/varray[@name="basis"]')
    m_vec = []
    for v1 in node_vec.findall("./v"):
        m_vec.append([float(x)*Ang2Bohr for x in v1.text.split()])

    latt.ReadFromRawCellVector(m_vec)
#Atoms
    list_species = doc_atom.findall('./array[@name="atoms"]/set/rc') 
    list_pos = doc_struct.findall('./varray[@name="positions"]/v')

    list_atoms = [ [species.findall('./c')[0].text.strip()] + 
                [float(x) for x in pos.text.split()]
            for species,pos in zip(list_species,list_pos)]
    latt.AddAtom(list_atoms,unit="crystal",latt="prim")

    return latt

def vasp_create_poscar(latt,filename,unit="primitive",lattregion="conv", alat=1.0):
    '''
    Create a POSCAR from Lattice object

    :param latt: Lattice object
    :param unit: maybe "crystal" or "ang"
    :param alat: the length unit on the second line
    '''
    import list_utils as lu

    f = open(filename,'w')
    f.write(" Created by tmckit\n")
    f.write(" %14.7f\n" % 1.0)

#VASP require the volume of product of a1 \dot (a2 \times a3) is positive
#So we check first
    pcv = latt.PrimitiveCellVector
    vol = lu.f_List_dot3(pcv[0], lu.f_List_cross3(pcv[1], pcv[2]))
    if (vol < 0): #Reverse the third one to make it positive
        pcv[2] = lu.f_List_Op_Scalar(pcv[2], "*", -1)

#Convert unit from bohr to ang
    for vec in pcv:
        f.write("%14.7f %14.7f %14.7f\n" % tuple([x * Bohr2Ang for x in vec]))

#Sort atoms with species name
    list2 = latt.GetAtomList(unit,lattregion)
    list2.sort(key=lambda x:x[0])
    name = None
    list3 = []
    list_group = []
    for atom in list2:
        if (name is None):
            name = atom[0]
            list3.append(atom)
            continue
        if (name != atom[0]):
            list_group.append(list3)
            name = atom[0]
            list3 = [atom]
        else:
            list3.append(atom)
    list_group.append(list3)

    for group in list_group:
        f.write(" %s" % group[0][0])
    f.write("\n")
    for group in list_group:
        f.write(" %4i" % len(group))
    f.write("\n")

    if (unit == "prim" or unit == "primitive"):
        unit_vasp = "Direct"
    elif ( unit == "ang" or unit == "angstrom"):
        unit_vasp = "Cartesian"
    else:
        raise ValueError("Unsupported length unit %s" % unit)
    f.write("%s\n" % unit_vasp)

    for group in list_group:
        for atom in group:
            f.write("%14.7f %14.7f %14.7f\n" % tuple(atom[1:4]))

    return list_group

list_vaspfile_input = ["INCAR", "KPOINTS", "POSCAR", "POTCAR"]
list_vaspfile_text = list_vaspfile_input + ["OUTCAR", "DOSCAR", "EIGENVAL", "CONTCAR", "IBZKPT", "vasprun.xml"]
list_vaspfile_scf = list_vaspfile_text + ["CHGCAR"]
list_vaspfile_all = list_vaspfile_text + ["WAVECAR" + "WAVEDER"]
dic_vaspfile = {"input":list_vaspfile_input, "text":list_vaspfile_text, "scf":list_vaspfile_scf, "all":list_vaspfile_all}

def vasp_copy(dirname_src=None,suffix_src=None,dirname_dst=None,suffix_dst=None,mode="text"):
    '''
    Copy related VASP files from one folder to another with given suffix
    "input" mode: INCAR,KPOINTS,POSCAR,POTCAR
    "text" mode: "input" + OUTCAR,DOSCAR,EIGENVAL,CONTCAR,IBZKPT,vasprun.xml
    "scf" mode: above + CHGCAR
    "all" mode : above + WAVECAR
    :param dirname_src: The folder to copy files from
    :param suffix_src: a suffix append as -suffix after filenames in the source folder
    :param dirname_dst: The folder to copy files to
    :param suffix_dst: a suffix append as -suffix after filenames in the destination folder
    '''
    if (dirname_src == dirname_dst and suffix_src is suffix_dst):
        raise ValueError("Cannot copy to the same folder with the same names")

    dirname_src = f_file_ensuredir(dirname_src)
    dirname_dst = f_file_ensuredir(dirname_dst)

    if (suffix_src is None):
        suffix_src = ""
    else:
        suffix_src = "_" + suffix_src

    if (suffix_dst is None):
        suffix_dst = ""
    else:
        suffix_dst = "_" + suffix_dst

    for filename in dic_vaspfile[mode]:
        filename_src = os.path.join(dirname_src,filename+suffix_src)
        filename_dst = os.path.join(dirname_dst,filename+suffix_dst)
#       print(filename_src)
        if (os.path.exists(filename_src)):
#           print(filename_src)
            if (os.path.exists(filename_dst)):
                os.remove(filename_dst) #Avoid read-only 
            shutil.copy(filename_src,filename_dst)
#           print("Copy file: %s" % filename_dst)

    return

def vasp_save(dirname=None,suffix=None,mode="text"):
    '''
    Save vasp calculation results to specific folder.
    "input" mode: INCAR,KPOINTS,POSCAR,POTCAR
    "text" mode: "input" + OUTCAR,DOSCAR,EIGENVAL,CONTCAR,IBZKPT,vasprun.xml
    "scf" mode: above + CHGCAR
    "all" mode : above + WAVECAR
    :param dirname: The folder to copy files to
    :param suffix: a suffix append as -suffix after filenames
    :param mode: "text","scf" or "all", indicating which files to save
    '''
    if (dirname is None and suffix is None):
        raise ValueError("Cannot save to the same folder with the same names")

    vasp_copy(dirname_dst=dirname, suffix_dst=suffix, mode=mode)

    return

def vasp_load(dirname=None,suffix=None,mode="text"):
    '''
    Load vasp calculation results from specific folder.
    "input" mode: INCAR,KPOINTS,POSCAR,POTCAR
    "text" mode: "input" + OUTCAR,DOSCAR,EIGENVAL,CONTCAR,IBZKPT,vasprun.xml
    "scf" mode: above + CHGCAR
    "all" mode : above + WAVECAR
    :param dirname: The folder to copy files from
    :param suffix: a suffix append as -suffix after filenames
    :param mode: "text","scf" or "all", indicating which files to load
    '''
    if (dirname is None and suffix is None):
        raise ValueError("Cannot load from the same folder with the same names")

    vasp_copy(dirname_src=dirname, suffix_src=suffix, mode=mode)

    return

def vasp_check_finish(filename="OUTCAR"):
    '''
    Check the OUTCAR to see whether the calculations is finished

    :param filename: The OUTCAR file to read
    '''
    with open(filename, 'r') as f:
        for line in f:
            if ("General timing and accounting informations for this job" in line):
                return True
    return False

def vasp_readband_out(outfile='OUTCAR',debug=False):
  """
  Read band energies from OUTCAR
  """
  # first get some dimension parameters
  nk = vasp_getout('nk',outfile)
  nb = vasp_getout('nb',outfile)
  efer = vasp_getout('efer',outfile)


  lines = io_read_lines_tagged(outfile,'E-fermi :',-1,'-------------',mode=-1)
  del lines[0:2]
#  if debug: print lines 

  kvecs=[]
  ebands=[]
  for ik in range(nk):
    i = ik * (nb+3)
    kvecs.append([float(x) for x in lines[i].split()[3:6]])
    ebands.append([float(x.split()[1]) for x in lines[i+2:i+nb+2]])

  return kvecs,ebands,efer

def vasp_write_kpoints(nks,mode='G',sh=None):
  
  if mode == 'G':
    ofile = open('KPOINTS','w')
    ofile.write("K-Points\n")
    ofile.write("0 \n")
    ofile.write("Gamma \n")
    ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
    if sh is None:
      ofile.write("0 0 0\n")
    else:
      ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))

  elif mode == 'M':
    ofile = open('KPOINTS','w')
    ofile.write("K-Points\n")
    ofile.write("0 \n")
    ofile.write("Monkhorst-Pack \n")
    ofile.write("%d %d %d\n"%(nks[0],nks[1],nks[2]))
    if sh is None:
      ofile.write("0 0 0\n")
    else:
      ofile.write("%8.4f %8.4f %8.4f\n"%(sh[0],sh[1],sh[2]))

    ofile.close()


def vasp_run(np=1,out=None,err=None):
  """
  Run vasp 
  """

  ierr = 0
  if np == 1:
    vasp_cmd = "vasp"
  else:
    vasp_cmd = "mpirun -np %d vasp"%(np)

  failure,output = commands.getstatusoutput( vasp_cmd )
  
  if failure:
    print "ERROR when running " + vasp_cmd
    ierr = 1

  print output
  return ierr 

def vasp_getout(var0,outfile='OUTCAR',debug=False):
  """
  Extract information from vasp OUTCAR file
  """
  var = var0.lower()
  if var == 'nb':
    val = int(io_grep_lines(outfile,'NBANDS=',1,-1,debug=debug))
    if debug: print var+"= ",val
  elif var == 'nat':
    val = int(io_grep_lines(outfile,'NIONS =',1,-1,debug=debug)) 
  elif var == 'nk':
    val = int(  io_grep_lines(outfile,'NKPTS =',1,4,debug=debug))
  elif var == 'nsp':
    val = int(io_grep_lines(outfile,'ISPIN  =',1,3,debug=debug))
  elif var == 'cput':
    val = float(io_grep_lines(outfile,'Total CPU time used',-1,-1,debug=debug))
  elif var == 'mag':
    val = float(io_grep_lines(outfile,'number of electron',-1,-1,debug=debug))

  elif var == 'efer':
    val = float(io_grep_lines(outfile,'E-fermi',-1,3,debug=debug))
  elif var == 'epsm':
    lines = io_read_lines_tagged(outfile,"REAL DIELECTRIC FUNCTION",3)[0]
    tmp = lines[-1].split()
    val = []
    for i in range(3):
      val.append(float(tmp[i+1]))
  elif var == 'epsm_pol': # read the epsm from the polarization approach
    lines = io_read_lines_tagged(outfile,"MACROSCOPIC STATIC DIELECTRIC TENSOR",4,mode=1)
    val = []
    for i in range(3):
      tmp = lines[i+1].split()
      val.append(float(tmp[i]))

  elif var == 'etot':
    val = float(io_grep_lines(outfile,'energy(sigma->0)',-1,-1,debug=debug))
  elif var == 'vol':
    val = float(io_grep_lines(outfile,'volume of cell',-1,-1,debug=debug))
  
  elif var == 'lat':
    lines = io_read_lines_tagged(outfile,"direct lattice vectors",iop=3,mode=-1)
    val = []
    for i in range(3):
      xyz=[]
      for j in range(3): 
        xyz.append(float(lines[i].split()[j]))
      val.append(xyz)
  elif var == 'gap':
    kvecs,ebands,efer = vasp_readband_out(outfile,debug=debug,mode=1)
    val = band_gap_analysis([ebands],efer,kvecs)
  
  return val

