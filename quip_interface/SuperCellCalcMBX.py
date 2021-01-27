from __future__ import print_function
import numpy as np
from enum import Enum
import ase.units as ase_units
from ase.calculators.calculator import Calculator, all_changes
import logging
import ase
import quippy

import os
import sys
import shutil
dirpath = os.getcwd()

# TO DO: modify for computer...
root_functions=os.path.dirname(__file__)
#from . import mbx_functions
sys.path.insert(0,root_functions)
import mbx_functions
import unwrap_oh2_ch4
shutil.copy(root_functions+"/mbx_pbc.json",dirpath+"/mbx_pbc.json")

class SuperCellCalcMBX(Calculator):
    # 'Properties calculator can handle (energy, forces, ...)'
    implemented_properties = {'energy', 'free_energy', 'forces','virial','stress'}

    # 'Default parameters'
    default_parameters = {}

    def __init__(self, min_box_requirement=24.,my_json="mbx_pbc.json", **kwargs):
        """ Create a potential for a potential and atoms so that cellsize is larger than min_box_requirement
        :potential: the potential that we use for the calculation
        :min_box_requirement: the minimum box length needed
        """
        super(SuperCellCalcMBX, self).__init__(**kwargs)
        self.min_box_requirement = min_box_requirement
        self.my_json = my_json

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        """Do the calculation.

        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        The subclass implementation should first call this
        implementation to set the atoms attribute.
        """

        results_unwrapped = unwrap_oh2_ch4.unwrap_water_methane_not_ordered_check_bondlength(ase.Atoms(atoms.copy()))
        if results_unwrapped["error_bond_lengths"]:
            raise ValueError ("Error with bond lengths in molecules. Check ordering of atoms.")
        atoms_unwrapped = results_unwrapped["phase_mw"]
        
        #if not np.allclose(atoms.cell,atoms.cell*np.eye(3)):
        #    raise ValueError ("Error: MBX does not work with non-diagonal cells yet.")

        lattice = atoms.cell     
        s_i_list =[1,1,1]
        for i in range(3):
            lat_ii=lattice[i][i]
            s_i_list[i]= (int(np.floor(self.min_box_requirement/lat_ii))+1)
        atoms_supercell=atoms_unwrapped.copy()*(s_i_list)
        
        if "virial" in properties and not ("stress" in properties):
            properties.append("stress")
            properties.remove("virial")
            
        dict_pot = mbx_functions.create_mbx_potential_variables(atoms_supercell)
        this_pot_super = quippy.Potential("IP MBXPBC n_monomers_types={"+
                                          dict_pot["n_monomers_types"]+"} nmon="+str(dict_pot["nmon"])+
                                          " json_file="+self.my_json+" diagonal_virial=T",param_str="")
        
        this_pot_super.calculate(atoms_supercell, properties, system_changes)
        results_super =this_pot_super.results
        
        results={}
        getEnergy = True if 'energy' in properties else False
        getForces = True if 'forces' in properties else False
        getStress = True if 'stress' in properties else False
        
        if getEnergy:
            results["energy"] = results_super["energy"]/np.prod(s_i_list)
            results["free_energy"]=results["energy"]
        if getForces:
            results["forces"] = results_super["forces"][:len(atoms)]
        if getStress:
            results["stress"] = results_super["stress"]
            tmp_virial = quippy.fzeros((3,3))
            tmp_virial[:,:] = quippy.stress_matrix(-results['stress']*atoms.get_volume())
            results['virial'] = tmp_virial
             
        self.results = results

        this_pot_super.reset()

    def reset(self):
        """Clear all information from old calculation."""
        super(SuperCellCalcMBX, self).reset()
        self.results={}
