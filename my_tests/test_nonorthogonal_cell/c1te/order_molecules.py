import quippy
import quippy
import ase
import numpy as np
import os
import sys
from collections import OrderedDict
sys.path.insert(0,"/home/es732/Eszter/function_modules/")
import unwrap_oh2_ch4


this_file_name=sys.argv[1]
new_file_name=this_file_name.split(".xyz")[0]+"_ordered.xyz"
atoms_in=quippy.Atoms(this_file_name)
atoms=unwrap_oh2_ch4.order_molecules(atoms_in)["ordered_mw"]
ase.io.write(new_file_name,atoms) 
