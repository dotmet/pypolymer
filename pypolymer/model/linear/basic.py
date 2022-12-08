from pypolymer.model import Polymer
from pypolymer.model.polymerTopo import PolymerTopoInfo

import numpy as np
from numpy import array as _arr
import matplotlib.pyplot as plt

class Point(Polymer):
    
    def __init__(self, name='Point', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
 
    def create_topo(self, x, y, z):
        
        self.topoInfo = self.new_topo()
    
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
        self.topoInfo = None
        self.name = name
        self.idx = idx

    def create_topo(self, init_point=[0,0,0], direction=[0,0,1], bond_length=1, atoms=1):
        
        self.topoInfo = self.new_topo()
        coords = []
        def f(p0, dir, bl):
            p0, dir = _arr(p0), _arr(dir)
            return p0 + dir*bl

        p0 = init_point
        for i in range(atoms):
            coords.append(p0)
            p0 = self.calCoords(f, p0, direction, bond_length)
        
        self.topoInfo.atoms = atoms
        self.topoInfo.coords = _arr(coords)
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=atoms)
        self.topoInfo.bonds = self.createBonds(sid=0, eid=atoms)
        self.topoInfo.angles = self.createAngles(sid=0, eid=atoms)

        return self.topoInfo

class Ring(Polymer):

    def __init__(self, name='Ring', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
        self.closed = True

    def create_topo(self, radius=1, bond_length=1, atoms=0, plane='xoy'):
        
        self.topoInfo = self.new_topo()

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
        
class Rectangle(Polymer):
    

    def __init__(self, name='Rectangle', idx=0):
        super().__init__()
        self.topoInfo = None
        self.name = name
        self.idx = idx
        self.closed = False
        
    def create_topo(self, init_point=[0,0,0], atoms=0, la=2, lb=2, bond_length=1, plane='xoy', direction='x'):
        
        self.topoInfo = self.new_topo()

        if (la+lb)*2 < atoms*bond_length:
            raise ValueError('atoms*bond_length i larger than the arc length of rectangle')
        if la<bond_length:
            la = bond_length
        if lb<bond_length:
            lb = bond_length

        plane_ax = [plane[0], plane[-1]]
        if direction not in plane:
            raise ValueError(f'Wrong direction {direction} for plane {plane}')
        else:
            dveca = np.zeros(3)
            dveca[['x', 'y', 'z'].index(direction)] = 1
            
            plane_ax.remove(direction)
            dir2 = plane_ax[-1]
            dvecb = np.zeros(3)
            dvecb[['x', 'y', 'z'].index(dir2)] = 1
            
        rod = Rod()
        na_atoms = int(la/bond_length)+1
        nb_atoms = int(lb/bond_length)+1
        _atoms = atoms
        atoms = (na_atoms+nb_atoms)*2-4
        if _atoms <= 1:
            _atoms = atoms
        coords = np.zeros((atoms, 3))
        init_point = np.array(init_point, dtype=float)

        sid, eid = 0, na_atoms
        _rod1 = rod.create_topo(init_point, dveca, bond_length, na_atoms)
        coords[:eid] = _rod1.coords

        sid, eid = eid, eid+nb_atoms-1
        _rod2 = rod.create_topo(coords[sid-1], dvecb, bond_length, nb_atoms)
        coords[sid:eid] = _rod2.coords[1:]

        sid, eid = eid, eid+na_atoms-1
        _rod3 = rod.create_topo(coords[sid-1], -dveca, bond_length, na_atoms)
        coords[sid:eid] = _rod3.coords[1:]


        sid, eid = eid, eid+nb_atoms-2
        if nb_atoms>2:
            _rod4 = rod.create_topo(coords[sid-1], -dvecb, bond_length, nb_atoms)
            coords[sid:] = _rod4.coords[1:-1]
            
        self.topoInfo.atoms = _atoms
        self.topoInfo.coords = _arr(coords[:_atoms])
        self.topoInfo.atom_ids = self.createAtomIds(sid=0, eid=_atoms)
        self.topoInfo.bonds = self.createBonds(sid=0, eid=_atoms, closed=self.closed)
        self.topoInfo.angles = self.createAngles(sid=0, eid=_atoms, closed=self.closed)

        return self.topoInfo
