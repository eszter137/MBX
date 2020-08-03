import quippy
import ase
from ase.calculators.calculator import all_changes

import logging
logging.basicConfig(level=logging.WARNING)

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import os
import sys

sys.path.insert(0,"/home/es732/Eszter/function_modules/")
import mbx_functions

dirpath = os.getcwd()
print dirpath

from collections import OrderedDict



def calculate_numerical_virial_diag(pot, in_atoms, d=1e-2):
    ## modified from the function calculate_numerical_stress() at ase/calculators/calculator.py 
    """Calculate numerical stress using finite difference."""

    atoms=in_atoms.copy()
    stress = np.zeros((3, 3), dtype=float)

    cell = atoms.cell.copy()
    V = atoms.get_volume()
    for i in range(3):
        x = np.eye(3)
        x[i, i] += d
        atoms.set_cell(np.dot(cell, x), scale_atoms=True)
        atoms.set_calculator(pot)
        eplus = atoms.get_potential_energy()#force_consistent=True)
        pot.reset()

        x[i, i] -= 2 * d
        atoms.set_cell(np.dot(cell, x), scale_atoms=True)
        atoms.set_calculator(pot)
        eminus = atoms.get_potential_energy()#force_consistent=True)
        pot.reset()

        stress[i, i] = (eplus - eminus) / (2 * d * V)
        x[i, i] += d
    atoms.set_cell(cell, scale_atoms=True)
    virial= (-stress*atoms.get_volume())
    return virial

d_list=[0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.000000001]

#version_list=["mbx_pbc_no_ww_www_mww","mbx_pbc_no_pips","mbx_pbc"]
version_list=["mbx_pbc_no_2b","mbx_pbc_no_3b","mbx_pbc_no_pips","mbx_pbc"]
version_list=["mbx_pbc_no_pips"]
#version_list=["mbx_pbc_no_3b","mbx_pbc_no_pips","mbx_pbc"]

dict_version_d_virial=OrderedDict([(version,
                                    OrderedDict([(d, [[np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]])
                                                 for d in d_list]))
                                   for version in version_list])
dict_version_model_virial=OrderedDict([(version,[[np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN], [np.NaN, np.NaN, np.NaN]])
                                   for version in version_list])



for this_version in version_list[:2]:
    this_at_tmp=quippy.Atoms("test.xyz")
    dict_pot=mbx_functions.create_mbx_potential_variables(this_at_tmp)
    
    
    pot_MBX = quippy.Potential("IP MBXPBC diagonal_virial=F n_monomers_types={"+
                               dict_pot["n_monomers_types"]+"} nmon="+str(dict_pot["nmon"])+
                               " json_file="+this_version+".json",param_str="")
    
    at = this_at_tmp.copy()
    pot_MBX.calc(at,virial=True)
    dict_version_model_virial[this_version] = np.array( at.params["virial"])
    print dict_version_model_virial[this_version]
    np.savetxt("test_"+this_version+"_virial"+".txt",dict_version_model_virial[this_version])
    
    
    for tmp_d in [0.1,0.01,0.001,0.0001,0.00001, 0.000001, 0.0000001, 0.00000001, 0.000000001]:
        test_at = this_at_tmp.copy()
        test_at.set_calculator(pot_MBX)
        dict_version_d_virial[this_version][tmp_d] =calculate_numerical_virial_diag(pot=pot_MBX,in_atoms=this_at_tmp.copy(),d=tmp_d) 
        print "virial test with ", tmp_d, " done"
        #print dict_version_d_virial[this_version][tmp_d]
        np.savetxt("test_"+this_version+"_fd_virial_d"+str(tmp_d)+".txt",dict_version_d_virial[this_version][tmp_d])
    
    
    
