#!/usr/bin/env python
from math import *
from chem_utils import *

def f_cif_write(name,latt,basis,sym_name='P 1',sym_num=1,title=None):
  """
  Write crystal structure into cif format
  """
  cif_file = name.strip()+'.cif'
  ofile = open(cif_file,'w')

  ofile.write("data_a \n\n")
  ofile.write("_cell_length_a %12.6f\n"%(latt[0]))
  ofile.write("_cell_length_b %12.6f\n"%(latt[1]))
  ofile.write("_cell_length_c %12.6f\n"%(latt[2]))
  ofile.write("_cell_angle_alpha %12.6f\n"%(latt[3]))
  ofile.write("_cell_angle_beta  %12.6f\n"%(latt[4]))
  ofile.write("_cell_angle_gamma %12.6f\n"%(latt[5]))
  ofile.write("_symmetry_space_group_name_H-M '%s'\n"%(sym_name))
  ofile.write("_symmetry_Int_Tables_number %6d\n"%(sym_num))  
  ofile.write("""

loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
  1     'x, y, z'

loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
""")
  for ia in range(len(basis)):
    atom= basis[ia][0]
    (x,y,z) = (basis[ia][1][0],basis[ia][1][1],basis[ia][1][2])
    ofile.write("%5s %5s %12.6f %12.6f %12.6f\n"%(atom,atom,x,y,z))
  ofile.close()

