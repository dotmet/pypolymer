import numpy as np
import matplotlib.pyplot as plt
import sys
import os

####################################
file_name = 'dump.2500sc300b.lammpstrj'
sys_atoms = 5548
chain_length = 2500
rod_length = 10
radius = 15
####################################

def get_data(file_name):
    with open(file_name, 'r') as f:
        data = f.readlines()
        data_mat = list(map(lambda x:x.split(), data[-sys_atoms:]))
        f.close()
    return np.array(data_mat, dtype = float)

def get_chain(data):
    chain = data[np.where(data_mat[:,1]==3.0)[0]]
    return chain[np.argsort(chain[:,0], axis=0)]

def split_chain(chain):
    new_chain = []
    for i in range(0, chain_length, rod_length):
        new_chain.append(chain[i:i+rod_length,:])
    return new_chain

def parse_mol(mol):
    _mol = mol[np.where(mol[:,2]>0)[0]]
    return _mol

def _isOnSmallSphere(mol, R=14):
    x = mol[:, 2]
    if x.max()>R or x.min()<-R:
        return True
    else:
        return False

if __name__ == '__main__':
    data_mat = get_data(file_name)
    chain = get_chain(data_mat)
    new_chain = split_chain(chain)
    for mol in new_chain:
        _mol = parse_mol(mol)
        if len(_mol)==0:
            continue
        elif _isOnSmallSphere(_mol, R=radius):
            continue
        else:
            plt.plot(_mol[:,3], _mol[:,4], c='k', marker='o')
    plt.show()
