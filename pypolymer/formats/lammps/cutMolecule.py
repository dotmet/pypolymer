import numpy as np
from pylmp.File.datafile.utilTools import *
from pylmp.File.datafile.molecule import Molecule

def cutMolecule_1(Mol, nsubMol):
    
    mlist = []
    ini_molid = Mol.id
    total_atoms = Mol.atoms
    subMol_atoms = int(total_atoms/nsubMol)
    atom_type = Mol.atom_type
    
    for num_Mol in range(nsubMol):
        subMol_p = Mol.coord[subMol_atoms*num_Mol:subMol_atoms*(num_Mol+1)]
        subMol_b = genBonds(subMol_atoms, isClosed=False)
        subMol_a = genAngles(subMol_atoms, isClosed=False)
        subMol = Molecule(idx = ini_molid+num_Mol)
        subMol.atom_type = atom_type
        subMol.initialize(Prev_M=mlist[-1], Atoms=subMol_p, Bonds=subMol_b, Angles=subMol_a, atyp=Mol.atom_type, btyp=Mol.bond_type, agtyp=Mol.angle_type, dtyp=Mol.dihedral_type)
        mlist.append(subMol)
    return mlist

def _BondInRange(mol, R):
    distL = []
    for p1,p2 in zip(mol[:-1, 2:5], mol[1:, 2:5]):
        distL.append(np.linalg.norm(p1-p2))
    return True if np.max(np.array(distL))<=R else False

def cutMolecule(Mol, nMol, mlist, nsAtoms):

    ini_molid = Mol.id
    total_atoms = Mol.atoms
    subMol_atoms = int(total_atoms/nMol)
    atom_type = Mol.atom_type
    
    for num_Mol in range(nMol):
        subMol_p = Mol.coord[subMol_atoms*num_Mol:subMol_atoms*(num_Mol+1)]
        subMol_b = genBonds(subMol_atoms, isClosed=False)
        subMol_a = genAngles(subMol_atoms, isClosed=False)
        subMol = Molecule(idx = ini_molid+num_Mol)
        subMol.atom_type = atom_type
        subMol.initialize(Prev_M=mlist[-1], Atoms=subMol_p, Bonds=subMol_b, Angles=subMol_a, atyp=Mol.atom_type, btyp=Mol.bond_type, agtyp=Mol.angle_type, dtyp=Mol.dihedral_type)
        if subMol.atoms==nsAtoms and _BondInRange(subMol_p, 3):
            mlist.append(subMol)
        else:
            continue
    return mlist
