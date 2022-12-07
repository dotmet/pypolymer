import numpy as np
import matplotlib.pyplot as plt
import sys
import os

#####################################
file_name = 'dump.2500sc300b.lammpstrj'
mol_length = 10
radius = 14
color = '#808080'
#####################################

def get_data(file_name):
    data_mat = []
    num_atoms = 0
    with open(file_name, 'r') as f:
        data = f.readlines()
        sys_atoms = int(data[3])
        data_mat = list(map(lambda x:x.split(), data[-sys_atoms:]))
        f.close()
    return np.array(data_mat, dtype=float)
    
def get_chain(data):
    chain = data[np.where(data[:,1] == 3.0)[0]]
    return chain[np.argsort(chain[:,0], axis=0)]

def split_chain(chain):
    new_chain = []
    for i in range(0,chain.shape[0], mol_length):
        new_chain.append(chain[i:i+mol_length, :])
    return new_chain

def parse_mol(mol):
    _mol = mol[np.where(mol[:,2]>0)[0]]
    return _mol


def _isCoordOut(mol, axis='x', R=14):
    axiss = ['x', 'y', 'z']
    x = mol[:, axiss.index(axis)+2]
    if x.max()>R or x.min()<-R:
        return True
    else:
        return False

def _isCuted(phy, dist=np.pi/2):
    cut_id = 0
    if phy.max() - phy.min() >= dist:
        for i in range(mol_length-1):
            if abs(phy[i] - phy[i+1]) > dist:
                cut_id = i+1
    return cut_id

if __name__ == '__main__':
    data_mat = get_data(file_name)
    chain = get_chain(data_mat)
    chain = split_chain(chain)
    for mol in chain:
        _mol = parse_mol(mol)
    plt.show()

