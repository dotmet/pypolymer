import numpy as np
import matplotlib.pyplot as plt
import copy
import os

def mol_cross_boundary(mol, bonds, box, bound_info=['p','p','p']):

    # bonds = bonds[np.argsort(bonds[:,0])]
    parsed_index = [0]
    
    for bond in bonds:
        p0 = mol[bond[0]]
        p1 = mol[bond[1]]
        _p1 = copy.deepcopy(p1)
        vec = p1-p0
        absvec = np.abs(vec)
        if bond[1] not in parsed_index:
            for i in range(3):
                if np.abs(vec[i])>box[i]/2 and bound_info[i].lower()[-1]=='p':
                    _p1[i] = -(vec[i]/absvec[i])*box[i]+_p1[i]
            parsed_index.append(bond[1])
        if np.linalg.norm(_p1-p0)>20:
            print('Bond too long in boundary parsing!')
            print(bond[0], p0)
            print(bond[1], _p1)
            print(bond[1], p1)
            print('----')
            break
        mol[bond[1]] = _p1

    return mol
    
def gen_angles(sid=0, atoms=0, isClosed=False):
    angle = []
    for i in range(atoms-2):
        l = [i, i+1, i+2]
        angle.append(l)
    if isClosed:
        angle.append([atoms-2, atoms-1,0])
        angle.append([atoms-1, 0, 1])
    return np.array(angle, dtype=int)+sid

def gen_bonds(sid=0, atoms=100, isClosed=False):
    bond = []
    for i in range(atoms-1):
        l = [i, i+1]
        bond.append(l)
    if isClosed:
        bond.append([atoms-1, 0])
    return np.array(bond, dtype=int)+sid
    
def cut_polymer(topo, length=1, num=0):
    
    if length==1:
        topo.bonds= []
        topo.angles = []
    else:
        bonds = gen_bonds(sid=0, atoms=length)
        angles = gen_angles(sid=0, atoms=length)
        nmol = int(topo.coords.shape[0]/length)
        if num>0 and num<=nmol:
            nmol = num

        for i in range(1, nmol):
            bonds = np.vstack([bonds, gen_bonds(sid=i*length, atoms=length)])
            angles = np.vstack([angles, gen_angles(sid=i*length, atoms=length)])
        topo.coords = topo.coords[:nmol*length]
        topo.bonds = bonds
        topo.angles = angles
        
        print(f'{nmol} molecules generated!')
        
    return topo