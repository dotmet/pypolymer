import numpy as np
import matplotlib.pyplot as plt

from pylmp.Model.polymerTopo import PolymerTopoParser

class Polymer(object):

    def __init__(self, topoInfo=None):
        
        self.name = 'Polymer'
        self.type = 'Linear'
        self.atom_sid = 0

    def calCoords(self, f, *args):
        '''
        A method used to generate coordinates of atoms.
        '''
        return f(*args)

    def createPolymer(self):
        return PolymerTopoParser().init_TopoInfo()

    def printInfo(self):
        print('-------'+ self.name +' Info------')
        print('ID         |', self.id)
        print('Atoms      |', self.atoms)
        print('Bond Length|', self.bond_length)
        print('Radius     |', self.radius)
        
    def snap3D(self):
        figure = plt.figure()
        c = np.array(self.topoInfo.coords)
        ax = figure.gca(projection = '3d')
        ax.plot(c[:,0], c[:,1], c[:,2], marker='o')
        plt.show()
    
    def snap2D(self):
        figure = plt.figure()
        c = np.array(self.topoInfo.coords)
        ax = figure.gca(projection = '3d')
        ax.plot(c[:,0], c[:,1], [0]*c.shape[0], marker='o')
        plt.show()

    '''
        These two default methods are applied to linear and contineuous 
        chain polymer. For complex structures, you must reload these
        methods for your purpose.
    '''
    def createAtomIds(self, sid, eid):
        atom_ids = [i for i in range(sid, eid)]
        return np.array(atom_ids) + self.atom_sid

    def createBonds(self, sid, eid, closed=False):
        bonds = [[i, i+1] for i in range(sid, eid-1)]
        if closed:
            bonds.append([eid-1, sid])
        return np.array(bonds)+self.atom_sid

    def createAngles(self, sid, eid, closed=False):
        angles = [[i, i+1, i+2] for i in range(sid, eid-2)]
        if closed:
            angles += [[eid-2, eid-1, 0], [eid-1, 0, 1]]
        return np.array(angles)+self.atom_sid

    '''
    This method cannot be universal. It must be reload by subclass 
    before usage.
    '''
    def createDihedrals(self, sid, eid):
        return None


        
