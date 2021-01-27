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

def unwrap_water_methane_not_ordered(atoms):

    atoms_unwrapped = atoms.copy()
    for i,at in enumerate(atoms):
        if (atoms.numbers[i]==8):
            iO=i
            for iH in [iO+1,iO+2]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                disp = atoms_unwrapped.get_distance(iO, iH, mic=True, vector=True)
                atoms_unwrapped[iH].position = atoms_unwrapped[iO].position + disp

        if  (atoms.numbers[i]==6):
            iC=i
            for iH in [iC+1,iC+2,iC+3,iC+4]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                disp = atoms_unwrapped.get_distance(iC, iH, mic=True, vector=True)
                atoms_unwrapped[iH].position = atoms_unwrapped[iC].position + disp


    return {"phase_mw":atoms_unwrapped}

def unwrap_water_methane_not_ordered_check_bondlength(atoms,bond_max=2.5):

    tmp_error=False
    atoms_unwrapped = atoms.copy()
    for i,at in enumerate(atoms):
        if (atoms.numbers[i]==8):
            iO=i
            for iH in [iO+1,iO+2]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                disp = atoms_unwrapped.get_distance(iO, iH, mic=True, vector=True)
                if np.sum(disp**2.)>(bond_max**2.):
                    print "ERROR: problem with OH bond length"
                    tmp_error=True
                atoms_unwrapped[iH].position = atoms_unwrapped[iO].position + disp

        if  (atoms.numbers[i]==6):
            iC=i
            for iH in [iC+1,iC+2,iC+3,iC+4]:
                if not (atoms.numbers[iH]==1):
                    print "ERROR: check order of atoms"
                disp = atoms_unwrapped.get_distance(iC, iH, mic=True, vector=True)
                if np.sum(disp**2.)>(bond_max**2.):
                    print "ERROR: problem with CH bond length"
                    tmp_error=True
                atoms_unwrapped[iH].position = atoms_unwrapped[iC].position + disp

    return {"phase_mw":atoms_unwrapped,"error_bond_lengths": tmp_error}


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

def unbuild_ordered_supercell_h2o_ch4(atoms_ordered_super, s_inds_list):
    num_cells=np.product(s_inds_list)
    n_unit = len(atoms_ordered_super)/num_cells

    numH2O = np.sum(atoms_ordered_super.numbers==8)
    numH2O_unit = numH2O/num_cells
    numCH4 = np.sum(atoms_ordered_super.numbers==6)
    numCH4_unit = numCH4/num_cells

    atoms_unit_cell=ase.Atoms()
    atoms_unit_cell+=(atoms_ordered_super[:(numH2O_unit*3)].copy())
    ind_start_ch4 = numH2O*3
    atoms_unit_cell+=(atoms_ordered_super[ind_start_ch4:(ind_start_ch4+numCH4_unit*5)].copy())

    tmp_cell = np.eye(3)
    for i in range(3):
        tmp_cell[i][i]=atoms_ordered_super.cell[i][i]/s_inds_list[i]
    atoms_unit_cell.set_cell(tmp_cell)
    atoms_unit_cell.pbc=True
    return atoms_unit_cell

def unbuild_supercell_h2o_ch4(atoms_ordered_super, s_inds_list):
    num_cells=np.product(s_inds_list)
    n_unit = len(atoms_ordered_super)/num_cells

    atoms_unit_cell=ase.Atoms()
    atoms_unit_cell+=(atoms_ordered_super[:(n_unit)].copy())

    tmp_cell = np.eye(3)
    for i in range(3):
        tmp_cell[i][i]=atoms_ordered_super.cell[i][i]/s_inds_list[i]
    atoms_unit_cell.set_cell(tmp_cell)
    atoms_unit_cell.pbc=True
    return atoms_unit_cell




def unbuild_genice_supercell_h2o_ch4(orig_atoms_ordered_super, s_inds_list):

    atoms_ordered_super=orig_atoms_ordered_super.copy()
    #atoms_ordered_super.positions=atoms_ordered_super.positions-atoms_ordered_super.positions[0,:]
    num_cells=np.product(s_inds_list)
    numH2O = np.sum(atoms_ordered_super.numbers==8)
    numCH4 = np.sum(atoms_ordered_super.numbers==6)
    numH2O_unit = numH2O/num_cells
    numCH4_unit = numCH4/num_cells

    atoms_unit_cell=ase.Atoms()
    tmp_cell = np.eye(3)
    for i in range(3):
        tmp_cell[i][i]=atoms_ordered_super.cell[i][i]/s_inds_list[i]
    atoms_unit_cell.set_cell(tmp_cell)
    atoms_unit_cell.pbc=True
    
    atoms_unit_cell+=(atoms_ordered_super[:(numH2O_unit*3)].copy())

    
    for i_at,at in enumerate(atoms_ordered_super):
        is_within_unit_cell=True
        if at.number==6:
            for i in range(3):
                if not(at.position[i]<tmp_cell[i][i] and at.position[i] >= 0.):
                    is_within_unit_cell=False
            if is_within_unit_cell:
                atoms_unit_cell+=(atoms_ordered_super[i_at:(i_at+5)].copy())
                print i_at,", ",
    print ""
    new_numH2O=np.sum(atoms_unit_cell.numbers==8)
    new_numCH4=np.sum(atoms_unit_cell.numbers==6)
    new_num_at=len(atoms_unit_cell)
    
    if not(new_num_at==len(atoms_ordered_super)/num_cells):
        print "Error: problem with number of atoms" , new_num_at,len(atoms_ordered_super)/num_cells
        if not(new_numH2O==numH2O/num_cells):
            print "Error: problem with number of O atoms" , new_numH2O,numH2O/num_cells
        if not(new_numCH4==numCH4/num_cells):
            print "Error: problem with number of C atoms" , new_numCH4,numCH4/num_cells
    else:
        return atoms_unit_cell

def unbuild_genice_supercell_h2o_ch4_indC(orig_atoms_ordered_super, s_inds_list,c_indices):

    atoms_ordered_super=orig_atoms_ordered_super.copy()
    #atoms_ordered_super.positions=atoms_ordered_super.positions-atoms_ordered_super.positions[0,:]
    num_cells=np.product(s_inds_list)
    numH2O = np.sum(atoms_ordered_super.numbers==8)
    numCH4 = np.sum(atoms_ordered_super.numbers==6)
    numH2O_unit = numH2O/num_cells
    numCH4_unit = numCH4/num_cells

    atoms_unit_cell=ase.Atoms()
    tmp_cell = np.eye(3)
    for i in range(3):
        tmp_cell[i][i]=atoms_ordered_super.cell[i][i]/s_inds_list[i]
    atoms_unit_cell.set_cell(tmp_cell)
    atoms_unit_cell.pbc=True

    atoms_unit_cell+=(atoms_ordered_super[:(numH2O_unit*3)].copy())


    for i_at in c_indices:
        atoms_unit_cell+=(atoms_ordered_super[i_at:(i_at+5)].copy())
    new_numH2O=np.sum(atoms_unit_cell.numbers==8)
    new_numCH4=np.sum(atoms_unit_cell.numbers==6)
    new_num_at=len(atoms_unit_cell)

    if not(new_num_at==len(atoms_ordered_super)/num_cells):
        print "Error: problem with number of atoms" , new_num_at,len(atoms_ordered_super)/num_cells
        if not(new_numH2O==numH2O/num_cells):
            print "Error: problem with number of O atoms" , new_numH2O,numH2O/num_cells
        if not(new_numCH4==numCH4/num_cells):
            print "Error: problem with number of C atoms" , new_numCH4,numCH4/num_cells
    else:
        return atoms_unit_cell

