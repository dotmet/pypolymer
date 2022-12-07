from pypolymer.model import Polymer
from pypolymer.model.polymerTopo import PolymerTopoInfo

import numpy as np
from numpy import array as _arr
import matplotlib.pyplot as plt

class Point(Polymer):
    
    def __init__(self, name='Point', idx=0):
        super().__init__()
        self.topoInfo = PolymerTopoInfo()
        self.topoInfo.name = name
        self.topoInfo.idx = idx
 
    def createPoint(self, x, y, z):

        def f(*args):
            return [[args[0], args[1], args[2]]]
        
        self.topoInfo.coords = _arr(self.calCoords(f, x, y, z))
        self.topoInfo.atoms = 1
        self.topoInfo.angles = []
        self.topoInfo.bonds = []
        self.atom_ids = _arr([0])
        return self.topoInfo
        

class Rod(Polymer):

    def __init__(self, name='Rod', idx=0):
        super().__init__()
        self.topoInfo = PolymerTopoInfo()
        self.topoInfo.name = name
        self.topoInfo.idx = idx

    def createRod(self, init_point=[0,0,0], direct=[0,0,1], bond_length=1, atoms=1):
        
        coords = []
        def f(p0, dir, bl):
            p0, dir = _arr(p0), _arr(dir)
            return p0 + dir*bl

        p0 = init_point
        for i in range(atoms):
            coords.append(p0)
            p0 = self.calCoords(f, p0, direct, bond_length)
        
        self.topoInfo.atoms = atoms
        self.topoInfo.coords = _arr(coords)
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)
        self.topoInfo.bonds = self.createBonds(sid=0, eid=atoms)
        self.topoInfo.angles = self.createAngles(sid=0, eid=atoms)

        return self.topoInfo

class Ring(Polymer):

    def __init__(self, name='Ring', idx=0):
        super().__init__()
        self.topoInfo = PolymerTopoInfo()
        self.topoInfo.name = name
        self.topoInfo.idx = idx
        self.closed = True

    def createRing(self, radius=1, bond_length=1, atoms=0, plane='xoy'):

        def f(t, r, d):
            if d=='xoy':
                return [r*np.cos(t), r*np.sin(t), 0]
            elif d=='xoz':
                return [r*np.cos(t), 0, r*np.sin(t)]
            elif d=='yoz':
                return [0, r*np.cos(t), r*np.sin(t)]

        if atoms == 0:
            atoms = round((2*np.pi*radius)/bond_length)
        if radius == 0:
            radius = (bond_length*atoms)/(2*np.pi)
        rad_inc = 2*np.pi/atoms

        coords = []
        for i in range(atoms):
            coords.append(self.calCoords(f, i*rad_inc, radius, plane))

        self.topoInfo.atoms = atoms
        self.topoInfo.coords = _arr(coords)
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)
        self.topoInfo.bonds = self.createBonds(sid=0, eid=atoms, closed=self.closed)
        self.topoInfo.angles = self.createAngles(sid=0, eid=atoms, closed=self.closed)

        return self.topoInfo


if __name__ == '__main__':
    p = Point()
    ptopo = p.createPoint(1,2,3)
     
    R = Ring()
    rtopo = R.createRing(10, 1, 31)

    ptopo.coords 
    ptopo.angles
    ptopo.bonds

    rtopo.bonds
    rtopo.angles
    rtopo.coords
    
