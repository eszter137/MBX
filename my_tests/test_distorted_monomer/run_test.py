import quippy
import ase
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import numpy as np
import sys
sys.path.insert(0,"/home/es732/Eszter/function_modules")
import mbx_functions
from SuperCellCalcMBX import SuperCellCalcMBX

filename=sys.argv[1]
this_at = ase.io.read(filename,index="-1")

potMBX=SuperCellCalcMBX(min_box_requirement=18.)


    
potMBX.reset()
    
this_at.set_calculator(potMBX)
potMBX.calculate(atoms=this_at,properties=["energy","forces"],system_changes=True)
    
ase.io.write(filename.split(".xyz")[0]+"_calculated.xyz",this_at)
potMBX.reset()    
    
