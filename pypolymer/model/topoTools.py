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