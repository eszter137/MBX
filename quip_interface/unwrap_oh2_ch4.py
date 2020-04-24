import quippy
from ase.io import read
import quippy
import ase
import numpy as np
import os
import sys

def unwrap_water(atoms):
    atoms = atoms.copy()
    for ind in range(0, len(atoms), 3):

        disp1 = atoms.get_distance(ind, ind + 1, mic=True, vector=True)
        disp2 = atoms.get_distance(ind, ind + 2, mic=True, vector=True)

        if np.linalg.norm(disp1) > 2. or np.linalg.norm(disp1) > 2.:
            print('Please check the order of atoms!')

        atoms[ind + 1].position = atoms[ind].position + disp1
        atoms[ind + 2].position = atoms[ind].position + disp2

    return atoms

def unwrap_water_methane(atoms):

    numO = list(atoms.numbers).count(8)
    numC = list(atoms.numbers).count(6)
    atoms_unwrapped = atoms.copy()
    
    for i in range(numO):
        iO=i*3
        if not (atoms.numbers[iO]==8):
            print "ERROR: check order of atoms"
        for iH in [iO+1,iO+2]:
            if not (atoms.numbers[iH]==1):
                print "ERROR: check order of atoms"
            disp = atoms_unwrapped.get_distance(iO, iH, mic=True, vector=True)
            atoms_unwrapped[iH].position = atoms_unwrapped[iO].position + disp
    
    for i in range(numC):
        iC=numO*3+i*5
        if not (atoms.numbers[iC]==6):
            print "ERROR: check order of atoms"
        for iH in [iC+1,iC+2,iC+3,iC+4]:
            if not (atoms.numbers[iH]==1):
                print "ERROR: check order of atoms"
            disp = atoms_unwrapped.get_distance(iC, iH, mic=True, vector=True)
            atoms_unwrapped[iH].position = atoms_unwrapped[iC].position + disp
            
    atoms_unwrapped_ch4 = atoms_unwrapped[(3*numO):].copy()
    atoms_unwrapped_h2o = atoms_unwrapped[:(3*numO)].copy()
    
    return {"phase_mw":atoms_unwrapped, 
            "phase_m":atoms_unwrapped_ch4,
            "phase_w": atoms_unwrapped_h2o}


def order_molecules(atIn):
    atoms = quippy.Atoms(atIn.copy())
    atoms_ordered = quippy.Atoms(lattice=atoms.lattice)
    for i,at in enumerate(atoms):
        if  (atoms.numbers[i]==8):
            atoms_ordered.append(atoms[i:(i+1)])
            for iH in [i+1,i+2]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                atoms_ordered.append(atoms[iH:(iH+1)])
    for i,at in enumerate(atoms):
        if  (atoms.numbers[i]==6):
            atoms_ordered.append(atoms[i:(i+1)])
            for iH in [i+1,i+2,i+3,i+4]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                atoms_ordered.append(atoms[iH:(iH+1)])
    return {"ordered_mw":atoms_ordered}

