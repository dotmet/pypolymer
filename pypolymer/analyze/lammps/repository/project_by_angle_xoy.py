import argparse

import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.pyplot import figure as fig

def parse_args():
    parser = argparse.ArgumentParser(description='Project from Cartesian to Angles')
    parser.add_argument(
        '-f',
        type=str,
        default='dump.2500sc300b.lammpstrj',
        help='file name')
    parser.add_argument(
        '-t',
        type=int,
        default=3,
        help='which type atoms choosed to project')
    parser.add_argument(
        '-n',
        type=int,
        default=0,
        help='number of atoms choosed to project')
    parser.add_argument(
        '-ln',
        type=int,
        default=2500,
        help='length of total chain')
    parser.add_argument(
        '-sn',
        type=int,
        default=10,
        help='length of subchain')
    parser.add_argument(
        '-r',
        type=float,
        default=16.0,
        help='constraint a region')
    parser.add_argument(
        '-o',
        type=str,
        default='out.png',
        help='output file name [picture]')
    parser.add_argument(
        '-c',
        type=str,
        default='#808090',
        help='choose the display color of the picture')
    args = parser.parse_args()
    return args

def get_data(file_name):
    data_mat = []
    with open(file_name, 'r') as f:
        data = f.readlines()
        sys_atoms = int(data[3])
        data_mat = list(map(lambda x:x.split(), data[-sys_atoms:]))
        f.close()
    return np.array(data_mat, dtype=float)
    
def get_chain(data, t):
    chain = data[np.where(data[:,1] == t)[0]]
    return chain[np.argsort(chain[:,0], axis=0)]

def split_chain(chain, mol_length):
    new_chain = []
    for i in range(0, chain.shape[0], mol_length):
        new_chain.append(chain[i:i+mol_length, :])
    return new_chain

def parse_mol(mol):
    pass

def _isCoordOut(mol, axis='x', R=14):
    axiss = ['x', 'y', 'z']
    #x = mol[:, axiss.index(axis)+2]
    #if x.max()>R or x.min()<-R:
    _r = np.array([np.linalg.norm(mol[i,2:5]) for i in range(mol.shape[0])])
    if max(np.abs(_r))>R:
        return True
    else:
        return False

def _isCuted(phy, mol_length, dist=np.pi/2):
    cut_id = 0
    if phy.max() - phy.min() >= dist:
        for i in range(mol_length-1):
            if abs(phy[i] - phy[i+1]) > dist:
                cut_id = i+1
    return cut_id

def rotate(theta, phy, angle):
    new_phy = phy+angle
    cut_theta1 = theta[np.where(phy<np.pi/2+angle)]
    cut_theta2 = theta[np.where(phy>=np.pi/2+angle)]
    new_theta = np.concatenate([cut_theta2, cut_theta1], axis=0)
    return new_theta, new_phy

def main():
    args = parse_args()

    data_matrix = []
    try:
        data_matrix = get_data(args.f)
    except :
        raise ValueError('invalid filename')

    mol_length = args.sn
    radius = args.r
    rotate_angle = np.pi/3
    
    chain = get_chain(data_matrix, args.t)
    chain = split_chain(chain, mol_length)

    figure = fig(figsize=(2*np.pi,np.pi))
    #plot molecule
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    for mol in chain:
        if _isCoordOut(mol, R=radius):
            continue
        _r = np.array([np.linalg.norm(mol[i,2:5]) for i in range(mol.shape[0])])
        theta = np.arcsin(mol[:,4]/_r)
        phy = np.arctan(mol[:,3]/mol[:,2])
        phy = np.where(mol[:, 2]<0, np.pi+phy, phy)

        #theta, phy = rotate(theta, phy, rotate_angle)
      
        c_id = _isCuted(phy, mol_length, dist=np.pi/2)

        a = 1.3
        lw = 0.8
        
        color = args.c
        phy = phy+np.pi/2
        if c_id == 0:
            plt.plot(phy, theta, lw=lw, markersize=a, c=color, marker='o')
        else:
            plt.plot(phy[:c_id], theta[:c_id], lw=lw, markersize=a, c=color, marker='o')
            plt.plot(phy[c_id:], theta[c_id:], lw=lw, markersize=a, c=color, marker='o')
    fig_name = args.o
    plt.savefig(fig_name, dpi=400, orientation='landscape')
    plt.show()

if __name__ == '__main__':
    main()
    
