from pypolymer.model import Polymer
from pypolymer.model.linear import *

import numpy as np

class Sphere(Polymer):

    def __init__(self, name='Sphere', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
        self.r = 1
        self.haveCenter=False
    
    def create_topo(self, atoms=1, r=1):
        self.topoInfo = self.new_topo()
        self.topoInfo.coords = self.calCoords(atoms)
        self.topoInfo.bonds = []
        if self.haveCenter:
            atoms = atoms+1
            self.topoInfo.atoms = atoms
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)

        return self.topoInfo

    def calCoords(self, npoints):
        
        coords = [[0,0,0]] if self.haveCenter else []
        phi = np.pi * (3. - np.sqrt(5.))  # golden angle in radians

        for i in range(npoints):
            y = 1 - (i / float(npoints - 1)) * 2  # y goes from 1 to -1
            radius = np.sqrt(1 - y * y)  # radius at y

            theta = phi * i  # golden angle increment

            x = np.cos(theta) * radius
            z = np.sin(theta) * radius

            coords.append((x, y, z))

        return np.array(coords)*self.r

    ##
    def createAngles(self):
        pass

    def createBonds(self):
        pass

    def createDihedrals(self):
        pass
    ##
class Plane(Polymer):

    def __init__(self, name='plane', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
        self.r = 1
        self.haveCenter=False
        
    def create_topo(self, init_p=[0,0,0], atoms=0, la=10, lb=10, bond_length=1, direction='x', plane='xoy'):
        
        self.topoInfo = self.new_topo()
        
        if direction not in plane:
            raise ValueError(f'Wrong growth direction {direction} for plane {plane}')
            
        na, nb = int(la/bond_length), int(lb/bond_length)
        shift_vec = np.zeros(3)

        if atoms == 0:
            atoms = int(la/bond_length)*int(lb/bond_length)
        
        if direction==plane[0]:
            la, lb = la-1, 1
            nrect = int(nb/2)
            rdir = plane[-1]
            ratoms = 2*na
        else:
            la, lb = 1, lb-1
            nrect = int(na/2)
            rdir = plane[0]
            ratoms = 2*nb
            
        shift_vec[['x','y','z'].index(rdir)] = 2*bond_length
        
        rect = Rectangle()
        rect.closed = False
        topo = rect.create_topo(init_p, ratoms, la, lb, bond_length, plane, direction)
        coords = topo.coords

        for i in range(nrect):
            if i==0:
                continue
            else:
                _vec = i*shift_vec
                coords = np.vstack([coords, topo.coords+_vec])

        self.topoInfo.coords = coords[:atoms]
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)
        self.topoInfo.bonds = self.createBonds(sid=0, eid=atoms)
        self.topoInfo.angles = self.createAngles(sid=0, eid=atoms)
                
        return self.topoInfo
        
class Lattice(Polymer):

    def __init__(self, name='Lattice', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
        self.full_bonds = True
        
    def create_topo(self, init_p=[0,0,0], atoms=0, bond_length=1, box=[100, 100, 100]):
            
        self.topoInfo = self.new_topo()
        
        if len(box)==3:
            L = box
            move_vec = 0
        else:
            L = [box[1]-box[0], box[3]-box[2], box[5]-box[4]]
            move_vec = np.array([box[1]/2+box[0]/2, box[2]/2+box[3]/2, box[4]/2+box[5]/2])
        
        Lx, Ly, Lz = L
        move_vec += -np.array(L)/2
        
        plane = Plane()
        top = plane.create_topo(la=Lx, lb=Ly, bond_length=bond_length)
        coords = top.coords
        
        Nz = int(Lz/bond_length)
        
        if atoms==0:
            atoms = top.atoms*Nz
            
        for i in range(1, Nz):
            vec = np.array([0, 0, bond_length])*i
            _coords = top.coords+vec
            if i%2==1:
                _coords = _coords[::-1]
            coords = np.vstack([coords, _coords])
        
        self.topoInfo.coords = coords + move_vec
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)
        if self.full_bonds:
            self.topoInfo.bonds = self.createBonds(sid=0, eid=atoms)
            self.topoInfo.angles = self.createAngles(sid=0, eid=atoms)
        else:
            self.topoInfo.bonds = []
            self.topoInfo.angles = []
            
        return self.topoInfo