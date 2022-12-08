import numpy as np
from .parsePosition import *

import os

### Read data from files####################
def read_sphere(atoms):
    path = 'F:\\Lammps\\SphereData\\'
    position = read_position(path+'sphere'+str(atoms)+'\\a'+str(atoms)+'xyz.out')
    bonds = read_Data(path+'sphere'+str(atoms)+'\\bond.data')
    dihedrals = read_Data(path+'sphere'+str(atoms)+'\\dihedral.data')
    xyz_out = path+'s'+str(atoms)+'\\'+str(atoms)+'.xyz'
    if not os.path.exists(path+'s'+str(atoms)):
        os.system('mkdir '+ path + 's' + str(atoms))
    data_file = path + 's'+str(atoms)+'\\data.sphere'
    return position, bonds, dihedrals, data_file

### Get radius ###########
def get_radius(atoms):
    position, _, _, _ = read_sphere(atoms)
    radius = max(max(abs(np.array(position)).tolist()))
    return radius

### Move molecule############################ 
def move_molecule(position, direction):
    return np.array(position)+np.array(direction)

### Change atoms id#########################
def change_serial(lst1, atoms):
    return np.array(lst1)+atoms

####Generate rod's position, bonds, and angles####
def init_rod(initposi, direction=[1,0,0], num=11):
    position = []
    bonds = []
    angles = []
    step = np.linalg.norm(direction)
    for i in range(num):
        posi = np.array(list(map(lambda x:step*i*x, direction)))+np.array(initposi)
        position.append(posi.tolist())
        bonds.append([i+1, i+1, i+2])
        angles.append([i+1, i+1, i+2, i+3])
    return np.array(position), np.array(bonds[:-1]), np.array(angles[:-2])

####Generate rod as a molecular###########
# def generateRod(Prev_M, i_d=1, p0=[0,0,0], direction=[0,1,0], num=11):
#     # rod = Molecule(idx=i_d)
#     rod.atom_type = 3
#     p,b,a = init_rod(p0, direction, num)
#     b = change_serial(b, Prev_M.atom_id)
#     a = change_serial(a, Prev_M.atom_id)
#     rod.initialize(Prev_M=Prev_M, Atoms=p, Bonds=b, Angles=a, atyp=3, btyp=2, agtyp=1)
#     return rod

def genAngles(atoms=100, isClosed=True):
    angle = []
    for i in range(atoms-2):
        l = [i+1, i+1, i+2, i+3]
        angle.append(l)
    if isClosed:
        angle.append([atoms-1,atoms-1,atoms,1])
        angle.append([atoms, atoms, 1, 2])
    return np.array(angle)

def genBonds(atoms=100, isClosed=True):
    bond = []
    for i in range(atoms-1):
        l = [i+1, i+1, i+2]
        bond.append(l)
    if isClosed:
        bond.append([atoms, atoms, 1])
    return np.array(bond)

def genDihedrals(faces):
    dihedrals = []
    faces = np.array(faces, dtype=int)
    for nf in range(faces.shape[0]):
        for _nf in range(nf+1, faces.shape[0]):
            f1, f2 = faces[nf, :], faces[_nf, :]
            _a = np.hstack([f1, f2])
            info, idxs, cons = np.unique(_a, return_index=True, return_counts=True)
            if len(info) == 2:
                dihe = [-1]*4
                a,b = info
                ida1, idb1 = np.where(f1==a)[0][0], np.where(f1==b)[0][0]
                ida2, idb2 = np.where(f2==a)[0][0], np.where(f2==b)[0][0]
                if (idb1-ida1)*(idb2-ida2)<0:
                    f2 = f2[:,:,-1]
                    ida2, idb2 = idb2, ida2
                id1, id2 = ida1-1, idb2+1
                id1 = id1 if id1>=0 else len(f1)
                id2 = id2 if id2<len(f2) else 0
                dihe[0] = f1[id1]
                dihe[1] = a
                dihe[2] = b
                dihe[3] = f2[id2]
                dihedrals.append(dihe)
    return dihedrals
    
def divideCircle(radius=1, parts=2, flat_vector=[0,0,1]):
    coordinates = []
    inc_rad = 2*np.pi/parts
    for i in range(parts):
        x = radius*np.cos(i*inc_rad)
        y = radius*np.sin(i*inc_rad)
        z = 0
        coordinates.append([x, y, z])
    return coordinates