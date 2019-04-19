from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate
from ase.io import read, write
from ase.calculators.vasp import Vasp
import sys, os, shutil, glob, os.path


def ase_set_runvasp(nproc):
    # get the current working directory
    wdir = os.getcwd()
    print "wdir=", wdir

    ### create a local run_vasp.py if not existing 
    if not os.path.isfile("run_vasp.py"): 
        ofile = open("run_vasp.py", "w")
        ofile.write("import os\n")
        if nproc <= 1:
            ofile.write("exitcode = os.system('vasp')\n")
        else:
            ofile.write("exitcode = os.system('mpirun -np %d vasp')\n" % (nproc))

        ofile.close()

    os.environ["VASP_SCRIPT"] = wdir.strip() + "/run_vasp.py"
    os.system("echo $VASP_SCRIPT")


def ase_run_vasp(mol,calc,job_dir,init_only=False):
    """
    This is used to set a vasp calculation in an ASE-based script 
      mode = 0: only create input files needed to run a vasp calculation 
      mode 
    """
    wdir = os.getcwd()
    if init_only:
        if not os.path.isdir(job_dir):
            os.mkdir(job_dir)  
        os.chdir(job_dir)
        mol.set_calculator(calc)

        calc.initialize(mol)
        calc.write_input(mol)
        os.chdir(wdir) 
        return 0.0 

    if not os.path.isdir(job_dir):
        os.mkdir(job_dir)
        os.chdir(job_dir)
        mol.set_calculator(calc)
    else:
        os.chdir(job_dir)
        calc = Vasp(restart=True)
        mol = calc.get_atoms()

    try:
        ene = mol.get_potential_energy()
    except:
        print "ERROR: Fail when running vasp in "+job_dir  
        os.system("touch ERROR")
        ene = 0.0 

    os.chdir(wdir) 

    return ene 

