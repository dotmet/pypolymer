import numpy as np
import matplotlib.pyplot as plt

### Some small usefull func ###############
def toArray(p):
    return np.array(p, dtype=float)

def cal_distance(p, q):
    dist = 0
    try:
        dist = np.linalg.norm(p-q)
    except:
        try:
            p, q = toArray(p), toArray(q)
            dist = np.linalg.norm(p-q)
        except:
            raise ValueError('Invalid input points')
    return dist

def cal_cosine(vec1, vec2):
    try:
        d_norm = np.linalg.norm(vec1)
        v_norm = np.linalg.norm(vec2)
    except:
        raise ValueError('Invalid vectors!')
    return np.dot(vec1, vec2)/(d_norm*v_norm)

def cal_theta(vec):
    vec = np.array(vec, dtype=float)
    vec = np.array([vec], dtype=float) if len(vec.shape)==1 else vec
    theta = np.arctan(vec[:,-1]/vec[:,-2])
    theta = np.where(vec[:,-2]<0, theta+np.pi, theta)
    return theta

def cal_tangent_on_sphere(point, R, axis='x', center=[0,0,0]):
    d_id = ['x', 'y', 'z'].index(axis)
    p1 = [[1,0,0],[0,1,0],[0,0,1]][d_id]
    direct = R*(toArray(p1)-toArray(center))
    pass
    
def count_group_1dim(data_list, gap, adjust=False):
    dl = np.array(data_list, dtype=float)
    dl = np.sort(dl)
    group = 1
    _list, __list, k = [], [], 0
    for a,b,i in zip(dl[:-1], dl[1:],range(1, dl.shape[0])):
        if np.linalg.norm(a-b) >= gap:
            group += 1
            _list.append(i-k)
            __list.append(dl[k:i])
            k=i
        if i == len(dl)-1:
            if i-k>=0:
                _list.append(i-k+1)
                __list.append(dl[k:])
        
    return group,_list,__list
###########################################

### Some functions used to determinate ###
def isInRange(point, R, center=[0,0,0]):
    STATUS = True
    if np.linalg.norm(np.array(point)-np.array(center)) > R:
        STATUS = False
    return STATUS

def isCuted(theta, mol_length, dist=np.pi/2):
    cut_id = 0
    theta = toArray(theta)
    if theta.max() - theta.min() > dist:
        for i in range(mol_length-1):
            if abs(theta[i]-theta[i+1])>dist:
                cut_id = i+1
    return cut_id

def ismolInSphere(mol, R, center=[0,0,0]):
    STATUS = True
    for x,y,z in zip(mol[:,2], mol[:,3], mol[:,4]):
        if cal_distance([x,y,z], center)>R:
            STATUS = False
    return STATUS

def ismolHeadInRange(mol, R, center=[0,0,0]):
    STATUS = False
    p1,p2 = mol[0,2:5], mol[-1,2:5]
    if isInRange(p1, R, center) or isInRange(p2, R, center):
        STATUS = True
    return STATUS

def iswholeMolInRange(mol, R, center=[0,0,0], mol_length=5, k=1):
    '''
        Use it to eliminate chain that position
    around small sphere.
    '''
    STATUS = False
    times = 0
    for x,y,z in zip(mol[:,2], mol[:,3], mol[:,4]):
        if cal_distance([x,y,z], center)<R:
            times += 1
    if times >= mol_length-k:
        STATUS = True
    return STATUS
##########################################

### Some functions used to filter data ###
def atom_type_filter(data, typ=3.0):
    _new_data = data[np.where(data[:,1]==typ)[0]]
    return _new_data[np.argsort(_new_data[:,0], axis=0)]

def split_to_molecule(data, mol_length=1):
    new_data = []
    for i in range(0, data.shape[0], mol_length):
        new_data.append(data[i:i+mol_length, :])
    return new_data

def get_mol_vec(mol, center, direction='out'):
    ps, pe=mol[0,2:5], mol[-1,2:5]
    vec = toArray(ps) - toArray(pe)
    ds = cal_distance(center, ps)
    de = cal_distance(center, pe)
    if direction=='out':
        vec = -vec if ds<de else vec
    elif direction=='in':
        vec = vec if ds<de else -vec
    return vec

def sort_mol_on_circle(chain, vec, clockwise='+'):
    pass
def mol_filter(mol):
    pass

#########################################

### Two functions used to convert coordinates ###
def mapping_to_angle_coord(mol, radius):
    phy = np.arccos(mol[:,4]/radius) # plane angle
    theta = np.arctan(mol[:,3]/mol[:,2])
    theta = np.where(mol[:,2]<0, np.pi+theta, theta) # 3d angle
    return theta, phy

def project_to_2d(mol, axis='x', direction='+'):
    column = ['x', 'y', 'z'].index(axis)+2
    if direction=='+':
        mol = mol[np.where(mol[:,column]>=0)[0]]
    elif direction=='-':
        mol = mol[np.where(mol[:,column]<=0)[0]]
    else:
        raise ValueError('Invlid direction')
    mol[:, column] = 0
    return mol
################################################ 
