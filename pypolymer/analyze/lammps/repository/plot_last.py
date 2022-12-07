import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure as fig
import os

from pylmp.Data.LmpTrj.util import *
from pylmp.Data.LmpTrj.trjdata import TrjData

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
        default=21,
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
    
def parse_chain(chain, veca):
    vec_mol = chain[0,2:5]-chain[-1,2:5]
    cosc = np.dot(vec_mol, veca)/(np.linalg.norm(vec_mol)*np.linalg.norm(veca))
    return cosc
        
def main():
    args = parse_args()
    fn = args.f
    mol_length = args.sn
    radius = args.r
    fig_name = args.o
    color = args.c

    trjdata = TrjData(fn, LAST_STEP=True)
    #data_mat = trjdata.split_by_steps(step_length=10)
    data_mat = [trjdata.coords_matrix]
    print(len(data_mat))
    dist = 5
    figure = fig(figsize=(6,6))
    a=2
    lw=1

    #plt.subplot(121)

    groups_list = []
    for wholemol in data_mat:
        chain = atom_type_filter(wholemol, typ=3.0)
        chain = split_to_molecule(chain, mol_length)
        final_chain = []
        chain_vec = []
        for mol,step in zip(chain, range(1,len(chain)+1)):
            if not ismolInSphere(mol, R=radius, center=[0,0,0]):
                continue
            if iswholeMolInRange(mol, dist, center=[0,0,0], mol_length=mol_length):
                continue
            _mol = project_to_2d(mol, axis='x', direction='+')
            if len(_mol)<10:
                continue
            elif not ismolHeadInRange(_mol, R=dist, center=[0,0,0]):
                continue
            else:
                final_chain .append(mol)
                mol_vec = get_mol_vec(_mol, center=[0,0,0], direction='out')
                chain_vec.append(mol_vec)
                
        chain_vec = np.array(chain_vec)
        theta = np.arctan(chain_vec[:,2]/chain_vec[:,1])
        theta = np.where(chain_vec[:,1]<0, theta+np.pi, theta)
        groups,_,_ = count_group_1dim(theta, gap=np.pi/10)
        groups_list.append(groups)
        plt.plot([i for i in range(len(theta))], np.sort(theta)+np.pi/2, marker='o')
    #plt.subplot(122)
    #plt.plot([i for i in range(len(theta))], np.cos(np.sort(theta)), marker='o')
    plt.savefig(fig_name, dpi=400, orientation='landscape')
    plt.show()

if __name__ == '__main__':
    main()
    
