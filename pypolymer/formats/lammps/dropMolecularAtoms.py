import numpy as np
from pylmp.File.datafile.utilTools import *
from pylmp.File.datafile.molecule import Molecule

def dropMolecularAtoms(Mol, atom_id=[]):
    if atom_id == []:
        return Mol
    mol = Mol
    
    newMol = Molecule(idx=Mol.id)
    
    coords = mol.coord
    bonds = mol.bonds
    angles = mol.angles
    dihedrals = mol.dihedrals
    
    bonds1 = bonds
    angles1 = angles
    dihedrals1 = dihedrals
    
    tl = []
    tl.append(Mol.atom_type)
    tl.append(Mol.bond_type)
    tl.append(Mol.dihedral_type)
    
    coords = np.delete(coords, np.array(atom_id)-1, axis=0)
    for aid in atom_id:
        bonds1, bonds = dropAtoms(bonds1, bonds, aid)
        angles1, angles = dropAtoms(angles1, angles, aid)
        dihedrals1, dihedrals = dropAtoms(dihedrals1, dihedrals, aid)
        
    newMol.initialize(Prev_M=Mol.prevM, Atoms=coords, Bonds=bonds, Dihedrals=dihedrals, atyp=tl[0], btyp=tl[1], dtyp=tl[2])
    return newMol
    
def dropAtoms(lst=[], lst1=[], atom_id=0):
    if len(lst) != 0:
        delist = np.where(np.where(lst[:, 1:] == atom_id, 1, 0).sum(axis=1)>0)[0]
        lst = np.delete(lst, delist, axis=0)
        lst1 = np.delete(lst1, delist, axis=0)
        try:
            lst1 = np.where(lst1 >= atom_id, lst1-1, lst1)
        except:
            pass
    return lst, lst1

def NeighList(bonds=[], atom_id=0):
    neighlist = []
    if len(bonds) != 0:
        lst = np.where(bonds[:, 1] == atom_id)[0]
        for i in lst:
            neighlist.append(bonds[i, 2])
        if(len(lst)<6):
            lst = np.where(bonds[:, 2] == atom_id)[0]
            for j in lst:
                neighlist.append(bonds[j, 1])
    return neighlist

def dropMolecularAtomsByRounds(Mol, start_id=1, rounds=0):
    droplist = [start_id]
    roundlist = [start_id]
    for r in range(rounds):
        bonds = Mol.bonds
        ndl = []
        for j in roundlist:
            ndl += NeighList(bonds, j)
        roundlist = ndl
        droplist += roundlist
    droplist=sorted(list(set(droplist)), reverse=True)
    print('Drop atoms:\n', droplist)
    return dropMolecularAtoms(Mol, atom_id=droplist)

#     bond_list = Mol.bonds
#     if rounds == 0:
#         return dropMolecularAtoms(Mol)
#     else:
#         id_list = [1]
#         sublist = [1]
#         for _ in range(rounds):
#             sublist = getnextlist(bond_list, sublist)
#             id_list += sublist
#         id_list = list(set(id_list))
#         print(id_list)
#         return dropMolecularAtoms(Mol, atom_id=id_list)
    
# def getnextlist(bonds, id_l):
#     lst1 = []
#     for ida in id_l:
#         lst2 = []
#         id1 = np.where(bonds[:,1] == ida)
#         for j in id1:
#             for k in j:
#                 lst2.append(bonds[k, 2])
#         lst1 += lst2
#     return list(set(lst1))