#!/usr/bin/env python
## this file contains functions used for embedded cluster model (ECM) calculations 
# Abbreviations used in the variable and function names:
#   ECM: embedded cluster model 
#   PCA: point charge array
 
import sys,os,shutil
import commands,string
from math           import *
from numpy          import *
from subprocess     import *
from scipy.optimize import *
from latt_utils     import *
from struct_utils   import *

def f_PCA_PotenField(xyz,pca):
  """
  This function calculate the electrostatic potential at a given point by the PCA
  """
  nc = len(pca) 
  poten = 0.0
  force = [0.0,0.0,0.0]
  for ic in range(nc):
    ch    = pca[ic][0]
    xyz_c = pca[ic][1][:]
    dist = sqrt( (xyz[0]-xyz_c[0])**2 + (xyz[1]-xyz_c[1])**2 + (xyz[2]-xyz_c[2])**2 )  
    if dist < 1.e-10 :
      print "WARNING: the position to calculate V_es is occupied by one of the point charges" 
      print " ic=%d,charge=%10.4f, (x,y,z)=(%12.6f,%12.6f,%12.6f)" %(ic,ch,xyz_c[0],xyz_c[1],xyz_c[2])
      print " -- this point charge is not counted in V_es"
    else: 
      poten += ch/dist
      for i in range(3): force[i] += - ch*(xyz[i]-xyz_c[i])/dist**1.5

  return poten,force

def f_PCA_Dipole(pca):
  """
  This function calculate the total charge and dipole of a give point charge array (PCA)
  """
  nc= len(pca)
  charge = 0.0
  dipole = [0.0,0.0,0.0] 

  center=[0.0, 0.0, 0.0]
  chtot = 0.0
  for ic in range(nc):
    xyz= pca[ic][1][:]
    ch = abs(pca[ic][0])
    for i in range(3): center[i] += xyz[i]*ch
    chtot += ch
  for i in range(3): center[i] /= chtot

  for ic in range(nc):
    ch  = pca[ic][0]
    xyz = pca[ic][1][:]
    charge    += ch
    for i in range(3): dipole[i] += ch*(xyz[i]-center[i])

  return charge,dipole

def f_ECM_Write_View(cluster,pca=None,name=None):
  """
  This function write the cluster and pca to a Gaussian input file 
  format so that it can be viewed by Gaussian
  """
  if name is None: out_file = "ecm-plot.gjf"
  else: out_file = name 

  print """!  write cluster and point charges to a file that can be viewed  by GaussView:
    positive charges represented by He, and negative charges represented by Rn
  """
  ofile = open(out_file, "w")
  ofile.write("#p  B3LYP/6-31G* Test\n\n")
  ofile.write("Cluster model of MgO\n\n")
  ofile.write("0, 1\n")
  for i in range(len(cluster)):
    ofile.write("%-6s %12.6f %12.6f %12.6f\n" %(cluster[i][0],cluster[i][1][0],cluster[i][1][1],cluster[i][1][2]))

  if not pca is None:
    for i in range(len(pca)):
      if pca[i][0] > 0:
        ofile.write("%-6s %12.6f %12.6f %12.6f\n"%("He",pca[i][1][0],pca[i][1][1],pca[i][1][2]))
      else:
        ofile.write("%-6s %12.6f %12.6f %12.6f\n"%("Rn",pca[i][1][0],pca[i][1][1],pca[i][1][2]))

  ofile.close()


def f_ECM_Export(cluster,pca,form='gauv',name=None):
  """ 
  Write out cluster and the embedding charges into a file with a format defined by "form"
  Available options for form: 
    'gauv'   -- Gaussian input file format with positive and negative charges 
             replaced by "He" and "Rn" respectively
    'gjf'  -- Gaussian input file format 
    'aims' -- FHI-aims geometry.in format   
  """
  nat = len(cluster)
  npc = len(pca)
  
  if name is None:  out_file = "ecm-xyz."+form
  else: out_file = name 

  ofile = open(out_file, "w")
    
  if form == 'gauv':
    print "Export cluster+pca geometry in to gjf format to be viewed in GaussView"
    ofile.write("#p  B3LYP/6-31G* Test\n\n")
    ofile.write("Cluster model of MgO\n\n")
    ofile.write("0, 1\n")
    for ia in range(nat):
      xyz = cluster[ia][1]
      atom = cluster[ia][0]
      ofile.write("%-6s %12.6f %12.6f %12.6f\n" %(atom,xyz[0],xyz[1],xyz[2]))
    for ic in range(npc):
      xyz = pca[ic][1]
      if pca[ic][0] > 0: atom = "He"
      else: atom = "Rn"
      ofile.write("%-6s %12.6f %12.6f %12.6f\n" %(atom,xyz[0],xyz[1],xyz[2]))

  elif form == 'aims': 
    print "Export cluster+pca geometry in to FHI-aims geometry.in format"
    for ia in range(nat):
      xyz = cluster[ia][1]
      atom = cluster[ia][0]
      ofile.write("%-12s %12.6f %12.6f %12.6f %6s\n" %("atom",xyz[0],xyz[1],xyz[2],atom))
    for ic in range(npc):
      ch = pca[ic][0]
      xyz = pca[ic][1]
      ofile.write("%-12s %12.6f %12.6f %12.6f 0 %12.6f\n" %("multipole",xyz[0],xyz[1],xyz[2],ch))
  elif form == "gjf":
    print "Export cluster+pca geometry in to Gaussian input file format"
    for ia in range(nat):
      xyz = cluster[ia][1]
      atom = cluster[ia][0]
      ofile.write("%-6s %12.6f %12.6f %12.6f\n" %(atom,xyz[0],xyz[1],xyz[2]))
    ofile.write("\n\n")
    for ic in range(npc):
      xyz = pca[ic][1]
      ch = pca[ic][0]
      ofile.write("%12.6f %12.6f %12.6f %12.6f\n" %(xyz[0],xyz[1],xyz[2],ch))
  else:
    print "ERROR: unsupported export format "+form
    return 1

  return 0 
  

def f_ECM_Write_PCA(pca,name=None):
  if name is None: out_file = "pca.xyz"
  else: out_file = name 

  ofile = open(out_file ,'w')
  for i in range(len(pca)):
    ofile.write("%12.6f %12.6f %12.6f %8.3f\n"%(pca[i][1][0],pca[i][1][1],pca[i][1][2],pca[i][0]))
  ofile.close()

def f_ECM_Write_Cluster(cluster,name=None):
  """
  Write the cluster to external files
  """
  if name is None: out_file = "cluster.xyz"
  else: out_file = name 

  print "\n!  write cluster to "+out_file 
  ofile = open(out_file,'w')
  for i in range(len(cluster)):
    atom = cluster[i][0]
    xyz = cluster[i][1]
    char = cluster[i][2]
    if char == 0.0: 
      ofile.write("%-4s -1 %12.6f %12.6f %12.6f\n" %(atom,xyz[0],xyz[1],xyz[2]))
    else:
      ofile.write("%-6s  %12.6f %12.6f %12.6f\n" %(atom,xyz[0],xyz[1],xyz[2]))
      
  ofile.close()

# make some post analysis 
def f_Post_Analysis(cluster,pca):
  print "\n============ post-analysis of generated point charge array =========="
  n_atoms = len(cluster)
  n_char = len(pca)
  print "  The number of atoms in the cluster:",len(cluster)
  print "  The number of point charges :",len(pca)

  nat_anion = 0
  nat_cation = 0
  char_tot = 0.0 
  for i in range(len(cluster)):
    atom = cluster[i][0]
    xyz = cluster[i][1]
    char =  cluster[i][2]
    print "%6s %12.6f %12.6f %12.6f # %10.2f" %(atom,xyz[0],xyz[1],xyz[2],char)
    if char > 0: nat_cation += 1
    else: nat_anion += 1
    char_tot += char 
  print "\n The number of cations:%6d" %(nat_cation)
  print " The number of anions: %6d" %(nat_anion) 
  print " The total number of charges in the cluster: %10.2f\n" %(char_tot)
  print "\nApparent chemical formula for the cluster:"+f_Check_Formula(cluster)

  if len(pca)>0: 
    charge,dipole = f_PCA_Dipole(pca)
    print "  Total charge of the point charge array:",charge
    print "  Dipole of the point charge array:",dipole
    print "  Electrostatic analysis for each atom in the cluster:"
    for iat in range(len(cluster)):
      ch = cluster[iat][2]
      if ch == 0.0: continue 
      atom = cluster[iat][0]
      xyz_at=cluster[iat][1][:]
      poten,field = f_PCA_PotenField(xyz_at,pca)
      print "  Atom= %2s: V_es = %12.6f E_es=(%12.6f,%12.6f,%12.6f)"%(atom,poten,field[0],field[1],field[2])

def f_Check_Position(iat,r_atom,core,rbuf,buffers):
  """
  Check the position of the current atom
    several possibility: 
     0  -- already included in the core 
     1  -- the buffer region 
     2  -- in the PC region  
  """
  nat_c = len(core)

  dist_min = 1.e10
  for iat_c in range(nat_c):
    rat_in_core = core[iat_c][1][:] 
    dist = f_Distance(rat_in_core,r_atom) 
    if dist < dist_min: dist_min = dist

  if dist_min < zero_dist:
    ipos = 0
  elif dist_min > rbuf: 
    ipos = 2
  elif buffers is None or buffers[iat] is None:
    ipos = 2
  else: 
    ipos = 1
  
  return ipos,dist_min        


def f_Cluster_Embed(latt_type,alatt,core,basis,rmax,\
                    shape =0,     \
                    charges=None, \
                    buffers=None,\
                    rbuf=None,   \
                    buf_charged=True,\
                    origin=None,\
                    dir=None,\
                    depth=0.0,\
                    debug=False): 
  """
  This function is used to embed a cluster into a point charge array with the option to add buffer atoms  
  Atoms not included in "core", but with the "rbuf" distance are replaced by some buffer atoms as defined in "buffers"
  all other atoms satisfying ( r < rmax ) are replaced by point charges according to "charges"
  
    latt_type  (string):  lattice types (following the convention of wien2k)
    alatt[0:5] (float) :  lattice constants 
    basis[0:nat_b-1]   :  basis of the crystal, each element is made of [<atom>,[x,y,z]>]
    core [0:nat_c-1]   :  atoms in the core, each element of core is made of [<atom>,<[x,y,z]>,<charge>]
    rmax (float)       :  the cut-off radius for point charge embedding
  
  Optional: 
    shape              :  shape of the cluster ( 0 -- for sphercial, 1 for rectangular)
    charges[0:nat_b-1] :  charges to replace atoms in "basis" 
    buffers[0:nat_b-1] :  buffer atoms to replace atoms in "basis";
                          for elements with the "None" value, then use the corresponding point 
                          if "buffers" is set None, then no buffer atoms are added
    rbuf  (float)      :  atoms not in the core but within the distance of "rbuf" are replaced by buffer atoms, 
                          if not set, then rbuf is automatically determined as the nearest neighbouring distances 
                          between atoms as determined by "basis" 
    buf_charged (log)  : if true, the atoms in the buffer region are also included in the point charge array 
    dir[0:2]  (float)  :  direction along which to embed the cluster (in the standard [h,k,l] ==[hkl] format)
                          None for the bulk embedding 
    depth (float)      :  used for directional embedding, the distance betweening the origin and the top surface 
    origin[0:2] (float):  the origin with respect to which the embedding is performed 
    debug (logical)    :  print out more detailed information 
  """
  print "\n -------------------------- Start f_Cluster_Embed -----------------------------" 

  nat_b = len(basis)  # the number of atoms in the basis 
  nat_c = len(core)   # the number of atoms in the core 

  # if origin is not defined, then use the charge center as the origin
  if origin is None: origin = [0.0, 0.0, 0.0]
  
  latt_vec= f_Latt_Vectors(latt_type,alatt)

  # set direction vector to extract the cluster
  if not dir is None:  dir_vec = f_Latt_Direction(dir,alatt,latt_type)

  if not buffers is None and rbuf is None:
    rbuf = f_NN_Distance(core)*1.01

  # by default, the basis is given in terms of internal coordinates
  # so it is necessary to convert basis coordinates to cartesian coordinates
  r_basis = f_Latt_Struct_I2C(basis,alatt,0)
  r_origin = f_Intern_to_Cart(origin,alatt,0)  

  i_range=[]
  if shape ==0: 
    for i in range(3): i_range.append(f_Set_Range(rmax*4,alatt[i]))
  else:
    for i in range(3): i_range.append(f_Set_Range(rmax[i]*4,alatt[i]))
 
  if debug:   
    print "Lattice constatns:",alatt
    print "Basis =",basis
    print "R_origin =",origin 
    print "Rbuf =",rbuf
  
  cluster=[]
  pca=[]
  # first put atoms in the core into cluster 
  for ia_c in range(nat_c): 
    atom = core[ia_c][0]
    xyz = core[ia_c][1][:]
    ch = core[ia_c][2]
    cluster.append([atom,xyz,ch])

  for ix in i_range[0]:
    for iy in i_range[1]:
      for iz in i_range[2]:

        for iat in range(nat_b):

          atom = basis[iat][0]
          if charges is None:char = 0.0 
          else: char = charges[iat]

          r_atom = r_basis[iat][1][:]
          for i in range(3): r_atom[i] += ix*latt_vec[0][i]+iy*latt_vec[1][i]+iz*latt_vec[2][i]

          if not dir is None : # along a particular direction 
            proj = 0.0
            for i in range(3): proj += (r_atom[i]-r_origin[i])*dir_vec[i]
            if proj > depth: continue

          dist = f_Distance(r_origin,r_atom)

          l_out_of_range = False 
          if shape==0: 
            if dist > rmax: l_out_of_range = True  
          else:
            for i in range(3): 
              if abs(r_atom[i]-r_origin[i]) > rmax[i]: l_out_of_range = True; break  
          if l_out_of_range: continue  

          # check the position of the current atom
          # several possibility: 
          #    0  -- already included in the core 
          #    1  -- the buffer region 
          #    2  -- in the PC region  
          ipos,dist_min = f_Check_Position(iat,r_atom,core,rbuf,buffers)

          if ipos == 1: cluster.append( [buffers[iat],r_atom,0.0] )

          if not charges is None :
            if ipos == 2 or (ipos==1 and buf_charged): pca.append( [char,r_atom[:],atom] )

  print "\n--------------------------End f_Cluster_Embed  -------------------------" 
  return cluster,pca            

def f_Cluster_Extract(latt_type,alatt,basis,rmax,\
                      shape =0, \
                      dir=None,\
                      charges=None,\
                      depth=0.0,\
                      origin=None,\
                      debug=False): 
  """
  Extract a cluster from a lattice, the cluster is still represented in 
    latt_type        -- lattice type 
    alatt            -- lattice constants
    basis[0:nat-1]   -- atoms and their internal coordinates in the unit cell
    rmax             -- the radius that defines the extraction region 
    dir[0:2]         -- the surface on which a cluster is extracted 
                        = None - bulk cluster
                        = [h,k,l]  extract on (hkl) surface 
    origin[0:2]      -- the origin with respect to which a cluster is extracted  
    
  """

  print "\n------------------------------ Start f_Cluster_Extract ----------------------" 
  nat=len(basis)  # the number of atoms

  # if origin is not defined
  if origin is None: origin = [0.0, 0.0, 0.0]

  latt_vec= f_Latt_Vectors(latt_type,alatt)

  # set direction vector to extract the cluster
  if not dir is None: dir_vec = f_Latt_Direction(dir,alatt,latt_type)

  i_range = []
  if shape==0 :
    for i in range(3): i_range.append(f_Set_Range(rmax*4,alatt[i]))
  else:
    for i in range(3): i_range.append(f_Set_Range(rmax[i]*4,alatt[i]))

  # by default, the basis is given in terms of internal coordinates
  # so it is necessary to convert basis coordinates to cartesian coordinates
  r_basis = f_Latt_Struct_I2C(basis,alatt,0)
  r_origin = f_Intern_to_Cart(origin,alatt,0)  
   
  print "Lattice constatns:",alatt
  print "Lattice _Vectors:", latt_vec
  print "Basis(Cartesian):"
  for iat in range(nat): print "\t",r_basis[iat]
  print "R_origin(Cartesian):",r_origin 
  print "R_cutoff:",rmax

  cluster=[]
  for ix in i_range[0]:
    for iy in i_range[1]:
      for iz in i_range[2]:
        ixyz=[ix,iy,iz]
        for iat in range(nat):
          atom = basis[iat][0]

          if charges is None: char = 0.0 
          else:  char = charges[iat]

          r_atom = r_basis[iat][1][:]
          for i in range(3): r_atom[i] += ix*latt_vec[0][i]+iy*latt_vec[1][i]+iz*latt_vec[2][i]

          dist = f_Distance(r_origin,r_atom)
          l_out_of_range = False 
             
          if shape==0: 
            if dist > rmax: l_out_of_range = True
          else:
            for i in range(3):
              if abs(r_atom[i]-r_origin[i]) > rmax[i]: l_out_of_range = True
          if l_out_of_range: continue

          if not dir is None: 
            proj = 0.0 
            for i in range(3): proj += (r_atom[i]-r_origin[i])*dir_vec[i]
            if debug:
              print "%-4s (%10.4f,%10.4f,%10.4f) proj=%8.2f, dist=%10.4f" %(atom,r_atom[0],r_atom[1],r_atom[2],proj,dist)

            if proj > depth + zero_dist : continue 
          cluster.append([atom,r_atom,char])

  print "\nApparent chemical formula for the extracted cluster:"+f_Check_Formula(cluster)
  print "\n------------------------------ End of f_Cluster_Extract ----------------------" 
  return cluster           


def f_ECM_Outermost(pca,atom=None,origin=None,width=0.0):
  """
  This function finds the indices of outermost layers of point charges with given charge 
    pca[0:npc-1] --- point charge arraygs, each element == [ <char>, <xyz>, <atom> ]
  """
  npc = len(pca)
  if origin is None: origin = [0.0, 0.0, 0.0]

  # find the maximal distance of point charges from the origin 
  rmax = 0.0 
  for ic in range(npc):
    if atom != pca[ic][2]: continue  
    r = f_Distance(pca[ic][1],origin)
    if r > rmax:  rmax = r

  outer_pc =[] 
  for ic in range(npc):
    if atom != pca[ic][2]: continue 
    ch = pca[ic][0]
    xyz = pca[ic][1][:]
    r = f_Distance(xyz,origin)
    if r >= rmax-width and r <= rmax: outer_pc.append(ic)

  return outer_pc

def f_ECM_Total_Charge(cluster,pca):
  """
  This function returns the total charge of the ECM system including the formal charges 
  of the atoms in the cluster
  """
  ch_tot = 0.0 
  for ia in range(len(cluster)): ch_tot += cluster[ia][2]
  for ic in range(len(pca)): ch_tot += pca[ic][0]
  return ch_tot

  
def f_ECM_Neutralize(cluster,pca,scheme=0):
  """ 
  This function neutralizes the cluster plus pca by adjusting the charge of "ad_atom"
  if ad_atom is None, it is automatically adjusted in terms of the sign of total charge 
 
  """
  ch_tot = f_ECM_Total_Charge(cluster,pca) 
  print "Total Charge:",ch_tot

  npc = len(pca)
  
  atom_anion = None
  for ic in range(npc):
    if pca[ic][0] < 0: 
      atom_anion = pca[ic][2] 
      char_anion = pca[ic][0]
      break 
   
  for ic in range(npc):
    if pca[ic][0] > 0: 
      atom_cation = pca[ic][2]
      char_cation = pca[ic][0]
      break 

  outer_anion = f_ECM_Outermost(pca,atom=atom_anion,  width=0.5)
  outer_cation = f_ECM_Outermost(pca,atom=atom_cation, width=0.5)
  n_anion =  len(outer_anion  )
  n_cation = len(outer_cation ) 
  print "The number of outermost cation charges: %d" %(n_cation) 
  print "The number of outermost anion  charges: %d" %(n_anion) 

  dq = 1.0-ch_tot/(n_anion*char_anion+n_cation*char_cation)

  for i in range(n_anion):
    ic = outer_anion[i]
    pca[ic][0] *= dq 

  for i in range(n_cation):
    ic = outer_cation[i]
    pca[ic][0] *= dq 

