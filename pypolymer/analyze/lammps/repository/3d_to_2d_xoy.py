import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure as fig
import os

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
        default=5,
        help='length of subchain')
    parser.add_argument(
        '-r',
        type=float,
        default=18.0,
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
    parser.add_argument(
        '-d',
        type=str,
        default='+',
        help="choose the direction '+-'")
    args = parser.parse_args()
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
    
def get_chain(data):
    chain = data[np.where(data[:,1] == 3.0)[0]]
    return chain[np.argsort(chain[:,0], axis=0)]

def split_chain(chain, mol_length):
    new_chain = []
    for i in range(0, chain.shape[0], mol_length):
        new_chain.append(chain[i:i+mol_length, :])
    return new_chain

def parse_mol(mol, direction):
    if direction=='+':
        _mol = mol[np.where(mol[:,3]>0)[0]]
    else:
        _mol = mol[np.where(mol[:,3]<0)[0]]
    return _mol

def _isCoordOut(mol, axis='x', R=14):
    axiss = ['x', 'y', 'z']
    #x = mol[:, axiss.index(axis)+2]
    #if x.max()>R or x.min()<-R:
    mol= mol[:,2:5]
    if mol.max()>R or mol.min()<-R:
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

def main():
    args = parse_args()

    data_matrix = []
    try:
        data_matrix = get_data(args.f)
    except :
        raise ValueError('invalid filename')

    mol_length = args.sn
    radius = args.r
    chain = get_chain(data_matrix)
    chain = split_chain(chain, mol_length)
    
    #plot molecule
    figure = fig(figsize=(6,6))
    a = 2
    lw = 1
    fig_name = args.o
    direction = args.d

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    
    for mol in chain:
        _mol = parse_mol(mol, direction)
        if len(_mol) == 0:
            continue
        elif _isCoordOut(_mol, R=radius):
            continue
        else:
            color = args.c
            plt.plot(_mol[:,4], _mol[:,2], c=color, marker='o', lw=lw, markersize=a)
    plt.savefig(fig_name, dpi=400, orientation='landscape')
    plt.show()

if __name__ == '__main__':
    main()
    
