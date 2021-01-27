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
root_functions=os.path.dirname(__file__) #"/home/es732/Eszter/function_modules/"
#from . import mbx_functions
sys.path.insert(0,root_functions)
import mbx_functions
import unwrap_oh2_ch4
shutil.copy(root_functions+"/mbx.json",dirpath+"/mbx.json")

class CalcMBX(Calculator):
    # 'Properties calculator can handle (energy, forces, ...)'
    implemented_properties = {'energy', 'free_energy', 'forces','virial','stress'}

    # 'Default parameters'
    default_parameters = {}

    def __init__(self, my_json="mbx.json", **kwargs):
        """ Create a potential for a potential and atoms
        :potential: the potential that we use for the calculation
        """
        super(CalcMBX, self).__init__(**kwargs)
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


        lattice = atoms.cell     
        dict_pot = mbx_functions.create_mbx_potential_variables(atoms)
        this_pot = quippy.Potential("IP MBX n_monomers_types={"+
                                     dict_pot["n_monomers_types"]+"} nmon="+str(dict_pot["nmon"])+
                                     " json_file="+self.my_json,param_str="")
        
        this_pot.calculate(atoms, properties, system_changes)
        results_tmp =this_pot.results
        
        results={}
        getEnergy = True if 'energy' in properties else False
        getForces = True if 'forces' in properties else False
        getStress = True if 'stress' in properties else False
        
        if getEnergy:
            results["energy"] = results_tmp["energy"]
            results["free_energy"]=results["energy"]
        if getForces:
            results["forces"] = results_tmp["forces"]
        if getStress:
            results["stress"] = results_tmp["stress"]
            tmp_virial = quippy.fzeros((3,3))
            tmp_virial[:,:] = quippy.stress_matrix(-results['stress']*atoms.get_volume())
            results['virial'] = tmp_virial
             
        self.results = results

        this_pot.reset()

    def reset(self):
        """Clear all information from old calculation."""
        super(CalcMBX, self).reset()
        self.results={}
