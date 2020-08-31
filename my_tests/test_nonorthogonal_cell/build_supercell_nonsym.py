import quippy
from ase.calculators.calculator import all_changes
import quippy
import ase
import numpy as np
import os
import sys
from collections import OrderedDict
sys.path.insert(0,"/home/es732/programs/function_modules/")
import unwrap_oh2_ch4


this_file_name=sys.argv[1]
s_ind_list=[2,3,2]
atoms = ase.io.read(this_file_name)
atoms_unwrapped = unwrap_oh2_ch4.unwrap_water_methane(atoms)
atoms=ase.Atoms(atoms_unwrapped["phase_mw"] )
atoms = atoms.repeat(s_ind_list)
new_file_name=this_file_name.split(".xyz")[0]+"_"+str(s_ind_list[0])+"x"+str(s_ind_list[1])+"x"+str(s_ind_list[2])+".xyz"

ase.io.write(new_file_name,atoms) 
