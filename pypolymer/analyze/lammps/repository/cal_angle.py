import numpy as np
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Caculate cosine')
    parser.add_argument(
        '-f',
        type=str,
        default='dump.2500sc300b.lammpstrj',
        help='file name')
    parser.add_argument(
        '-sn',
        type=int,
        default=10,
        help='length of total chain')

    args = parser.parse_args()
    return args
    
################ Get Data ################
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

def chain_xyz(mol):
    return mol[:, 2:5]
###########################################

################ Tools ###################
def _isInRange(center, mol, dist):
    chain_xyz = mol[:, 2:5]
    for point in chain_xyz:
        _dist = cal_distance(center, point)
        if _dist <= dist:
            return point
        else:
            return False
    
def cal_distance(a,b):
    return np.linalg.norm(a-b)

def cal_cosine(direction, vector):
    d_norm = np.linalg.norm(direction)
    v_norm = np.linalg.norm(vector)
    return np.dot(direction, vector)/(d_norm*v_norm)

def cal_direction(R, point, axis=np.array([1,0,0])):
    n = 2*point
    ps = R*axis
    pe = -ps
    p = np.cross((ps-point), (pe-point))
    return np.cross(n,p)
###########################################

################ Parse Data ###############
def chain_filter(center, chain, dist):
    _chain = []
    for mol in chain:
        point = _isInRange(center, mol, dist)
        if point:
            _chain.append(mol)
        else:
            pass
    return _chain
############################################

def main():
    args = parse_args()

    mol_length = args.sn
    data_mat = []
    try:
        data_mat = get_data(args.f)
    except:
        raise ValueError('Invalid filename')
        
    chain = get_chain(data_mat)
    chain = split_chain(chain, mol_length)

    radius = 13.8
    mol_cosine = []  
    for mol in chain:
        mol_vec = mol[0, 2:5] - mol[-1, 2:5]
        axis = np.array([1,0,0])
        direction = cal_direction(axis=axis, R=radius, point=mol[0, 2:5])
        mol_cosine.append(cal_cosine(direction, mol_vec))

    print(len(mol_cosine))
    mean_cosine = np.mean(np.abs(mol_cosine))
    print(mean_cosine)
    
    center = [0,0,0]
    dist = 5
    chain = chain_filter(center, chain, dist)
    

if __name__ == '__main__':
    main()
