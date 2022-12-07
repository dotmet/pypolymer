from pypolymer.model import Polymer

import numpy as np

class Sphere(Polymer):

    def __init__(self, name='Sphere', idx=0):
        super().__init__()
        self.topoInfo = PolymerTopoInfo()
        self.topoInfo.name = name
        self.topoInfo.idx = idx
        self.r = 1
        self.haveCenter=False
    
    def createSphere(self, atoms=1, r=1):
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