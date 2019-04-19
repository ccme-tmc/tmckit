from ase import Atoms
from ase.constraints import FixAtoms
from ase.build import fcc111, hcp0001

mol_C = Atoms('C', positions=[(0., 0., 0.)])
mol_H = Atoms('H', positions=[(0., 0., 0.)])
mol_O = Atoms('O', positions=[(0., 0., 0.)])
mol_CO = Atoms('CO', positions=[(0., 0., 0.), (0., 0., 1.3)])
mol_OH = Atoms('OH', positions=[(0., 0., 0.), (0., 0., 1.1)])
mol_CH = Atoms('CH', positions=[(0., 0., 0.), (0., 0., 1.1)])
mol_CH2 = Atoms('CH2', positions=[(0., 0., 0.),(0.521, 0.322, 0.925),(-0.935,-0.552,0.422)])
mol_CH2_b = Atoms('CH2', positions=[(0., 0., 0.),(-0.92, 0.0, 0.60),(0.92,0,0.60)])

mol_CH3 = Atoms('CH3', positions=[
            (0., 0., 0.),
            (   0.895,   -0.517,    0.434),
            (  -0.895,   -0.517,    0.435),
            (   0.000,    1.034,    0.435)
            ])
mol_CHO    =  Atoms('CHO', positions=[(0., 0., 0.),(-0.373,-0.61, 0.866),(0.627,1.067,0.279)])
mol_COH    =  Atoms('COH', positions=[(0., 0., 0.),(0.627,1.067,0.279),(-0.373,-0.61, 0.866)])
mol_CO2    = Atoms('CO2', positions=[(0., 0., 0.),(-1.18, 0.0,0.0),(1.18, 0.0, 0.0)])
mol_CO2_ch = Atoms('CO2', positions=[(0., 0., 0.),(-0.92,0.0,0.78), (1.28,0.0,0.0)])

mol_CH4 = Atoms('CH4', positions=[
        (   0.000,    0.000,    0.000),
        (   0.519,   -0.741,    0.632),
        (  -0.519,    0.741,    0.632),
        (   0.738,    0.516,   -0.640),
        (  -0.738,   -0.516,   -0.640)
        ])

mol_set = {
        'C':    [mol_C,  1.24],
        'O':    [mol_O,  1.52],
        'H':    [mol_H,  1.13],
        'OH':   [mol_OH, 1.60],
        'CH':   [mol_CH, 1.35],
        'CH2':  [mol_CH2,1.52],
        'CH2-b':[mol_CH2_b,1.60],
        'CH3': [mol_CH3,1.82],
        'CH4': [mol_CH4,3.50],
        'CO':  [mol_CO, 1.53],
        'CHO': [mol_CHO,1.71],
        'COH': [mol_COH,1.71],
        'CO2':    [mol_CO2,3.50],
        'CO2-ch': [mol_CO2_ch,2.03],
        }

def create_bimetal_surf(composition,surf_type,size,latt='',fix_layers=2,vacuum=12.0):
    """
    Create bimetallic surfaces 
       composition:  "M1" - single-metal 
                     "M1,M2,config" - bimetallic surface with the configuration 
       latt: 

    """

    # extract lattice constants if available
    if latt == '':
        if composition=='Co' and surf_type == 'fcc111':
            a_latt = 3.544 
        else:
            a_latt, c_latt = (None, None)
    else:
        tmp = latt.split(',')
        a_latt = float(tmp[0])
        if len(tmp) > 1:
            c_latt = float(tmp[1])


    # set the surface type name 
    if surf_type == 'fcc111':
        surf_type_name = '111'
    elif surf_type == 'hcp0001':
        surf_type_name = '0001'

    # set the metallic surface slab
    tmp = composition.split(',')
    M1_name = tmp[0]
    if len(tmp) > 1:
        M2_name = tmp[1]
        config = tmp[2]
    else:
        M2_name = ''
        config = ''

    if not "bulk" in config:
        if surf_type == 'fcc111':
            surf = fcc111(M1_name, a=a_latt, size=size, vacuum=vacuum)
        elif surf_type =='hcp0001':
            surf = hcp0001(M1_name, a=a_latt, c=c_latt,  size=size, vacuum=vacuum)

        surf_name = M1_name + surf_type_name+"-%d%d%d"%(size[0],size[1],size[2])
        if M2_name != '':
            surf_name = surf_name + '-'+config_bimetal + "-" + M2_name

            if config_bimetal == 'top':
                for atom in surf:
                    if atom.tag == 1 :
                        atom.symbol = M2_name
            elif config_bimetal == 'sub':
                for atom in surf:
                    if atom.tag == 2 :
                        atom.symbol = M2_name
            elif config_bimetal == 'mixtop':
                for atom in surf:
                    if atom.tag == 1 :
                        atom.symbol = M2_name
                        break
            elif config_bimetal == 'mixsub':
                for atom in surf:
                    if atom.tag == 2 :
                        atom.symbol = M2_name
                        break
    else:
        print "ERROR: Not implemented yet for the bimetal configuration %s!!!"%(config) 

    # get the number of layers
    nlayers = max(surf.get_tags())

    mask = [atom.tag > nlayers - fix_layers for atom in surf]
    constr_bottom = FixAtoms(mask=mask)
    surf.set_constraint([constr_bottom])

    return surf, surf_name 

def parse_ads_inp(ads_inp):
    """
    Parse the adsorbate input information  
    """
    tmp_s = ads_inp.split(';') 
    ads_info = []
    for i in range( len(tmp_s)): 
        tmp = tmp_s[i].split(',')
        mol_name = tmp[0]
        site = tmp[1]
        h_ads = 2.0
        rot = 0.0 
        off = 0.0 
        for s in tmp[2:]:
            if 'h=' in s:
                t = s.split('=')
                h_ads = float(t[1]) 

            if 'r=' in s:
                t = s.split('=') 
                rot = float(t[1]) 

            if 'o=' in s:
                t1 = s.split('=') 
                t2 = t1[1].split('_') 
                if len(t2) == 1:
                    off = float(t2[0]) 
                else:
                    off = (float(t2[0]), float(t2[1]))
        ads_info.append([mol_name, site, h_ads, rot, off])

    return ads_info 

def calc_formation_energy(symbl, ene, ref_dict, surf_name='surf' ) : 
    """
    calculate the formation energy in terms of the symbol,  the energy and the reference energy dictory 
       <symbl> can CH2, CH2_s etc. 
    in ref_dict, the pure surface is denoted as 'slab'
    """

    from ase.atoms import string2symbols
    E0 = ene 
    if '_s' in symbl:
        name,site = symbl.split('_') #split key into name/site
        E0 -= ref_dict[surf_name]
    else:
        name = symbl

    #remove - from transition-states
    formula = name.replace('-','')

    #get the composition as a list of atomic species
    composition = string2symbols(formula)

    #for each atomic species, subtract off the reference energy
    for atom in composition:
        E0 -= ref_dict[atom]
        #round to 3 decimals since this is the accuracy of DFT
    E0 = round(E0,3)
    return E0 


def get_formation_energies(energy_dict, ref_dict, surf_name='111'):
    """
    calculate the formation energy based on the energy data given in energy_dict using the reference energies given in ref_dict
    """
    from ase.atoms import string2symbols
    formation_energies = {}
    for key in energy_dict.keys(): #iterate through keys
        E0 = energy_dict[key] #raw energy

        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)
            if site == surf_name :
                E0 -= ref_dict[site] #subtract slab energy if adsorbed

            #remove - from transition-states
            formula = name.replace('-','')

            #get the composition as a list of atomic species
            composition = string2symbols(formula)

            #for each atomic species, subtract off the reference energy
            for atom in composition:
                E0 -= ref_dict[atom]
            #round to 3 decimals since this is the accuracy of DFT
            E0 = round(E0,3)
            formation_energies[key] = E0
    return formation_energies

def make_mkm_input(file_name,surf_name, energy_dict, freq_dict = None):
    """
    Create the mkm input file used for the micro-kinectic modeling with CatMAP 
    """
    #create a header
    header = '\t'.join(['surface_name','site_name',
                        'species_name','formation_energy',
                        'frequencies','reference'])

    lines = [] #list of lines in the output
    for key in energy_dict.keys(): #iterate through keys
        E = energy_dict[key] #raw energy
        name,site = key.split('_') #split key into name/site
        if 'slab' not in name: #do not include empty site energy (0)

            if not freq_dict is None: 
                freq = freq_dict[key]
            else:
                freq = []

            if site == 'gas':
                surface = None
            else:
                surface = surf_name 

            outline = [surface,site,name,E,freq,'Input File Tutorial.']
            line = '\t'.join([str(w) for w in outline])
            lines.append(line)

    lines.sort() #The file is easier to read if sorted (optional)
    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    input = open(file_name,'w') #open the file name in write mode
    input.write(input_file) #write the text
    input.close() #close the file

    print 'Successfully created input file'

def create_pathway_diagram(pathway,energy_dict,surf_name='111'): 
    """
    create the data along a reaction pathway in terms of the formation energies given in energy_dict 
    e.g. pathway = ["CO+0.5*O2","CO-O","CO+O"+"CO2"] 
    """
    ene_diag = [] 
    for i in range(len(pathway)):
        spec_list = pathway[i].split('+') 
        ene = 0.0 
        for spec_in in spec_list:
            if '*' in spec_in:
                tmp = spec_in.split('*') 
                n_spec = int(tmp[0])
                spec = tmp[1] 
            else:
                n_spec = 1
                spec = spec_in

            if '_s' in spec:
                spec_name = spec.strip()[:-2]+'_'+surf_name
            else:
                spec_name = spec.strip()+'_gas'

            ene += n_spec * energy_dict[spec_name] 
        ene_diag.append(ene) 
    return ene_diag

def get_adsorbate_info(mol,surf,surf_mol):
    """
    extract the coordinates of adsorbate on the surface
    """

    nat_surf = len(surf) 
    nat_tot = len(surf_mol)
    nat_mol = len(mol)

    z_surf = 0.0
    for atom in surf_mol:
        if atom.symbol in surf.get_chemical_symbols():
            z_surf= max(atom.z,z_surf)

    # get the scaled position of the adsorbing atom
    pos_internal = surf_mol.get_scaled_positions()
    xyz_C = pos_internal[nat_tot - nat_mol ]
    print "The scaled position of adsorbing atom: %6.3f %6.3f %6.3f"%(xyz_C[0],xyz_C[1],xyz_C[2])

    atom = surf_mol[ nat_tot - nat_mol  ]
    xyz_C = [atom.x,atom.y,atom.z]
    print "The position of adsorbing atom: %6.3f %6.3f %6.3f"%(xyz_C[0],xyz_C[1],xyz_C[2])

    print "The height of the adsorbing atom with respect to the surface:%8.3f"%(xyz_C[2]-z_surf)
    print "The coordidates of %s w.r.t. the adsorbing atom:"%(mol_name)
    for ia in range(nat_tot - nat_mol, nat_tot):
        atom = surf_mol[ia]
        print "%-4s (%8.3f, %8.3f, %8.3f)"%(atom.symbol, atom.x-xyz_C[0], atom.y-xyz_C[1], atom.z-xyz_C[2])

    print "#%20s %12s %12s %12s %12s" % ("surf_mol", 'e_surf', 'e_mol', 'e_surf_mol',"ads_energy")
    print "%20s %12.6f %12.6f %12.6f  %12.3f" % (surf_mol_name, e_surf, e_mol, e_surf_mol, e_ads)
