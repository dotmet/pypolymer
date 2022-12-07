import numpy as np
from numpy import array as _arr

import matplotlib.pyplot as plt

from pypolymer.model import Polymer
from pypolymer.formats.lammps.utilTools import *
from pypolymer.model.linear import *
from pypolymer.model.shape import *

from pypolymer.model.polymerTopo import PolymerTopoInfo, PolymerTopoParser

# Store topological info
class BranchedTopoInfo(PolymerTopoInfo):
    
    def __init__(self):

        super().__init__()

        self.backbone_type = 'Rod'
        self.backbone_coords = []
        self.arms_types = []
        self.arms_ids = []
        self.arms_coords = {}

# Analysis topological structure
class BranchedTopoParser(PolymerTopoParser):

    def __init__(self):

        super().__init__()
    

# Constructure topological structure
class BranchedPolymer(Polymer):

    def __init__(self, name='Branched', idx=0):
        '''
            Star polymer is a sepcific type of Branched polymers.
            The backbone 
        '''
        self.topoInfo = BranchedTopoInfo()
        self.topoInfo.name = name
        self.topoInfo.idx = idx
        
        self.topoInfo.IS_BASE_POLYMER = False

        ### An example set of params for evenly distributed Linear Branched polymer.
        self.backparams = {
            'bond length': 1,

            'type': 'Rod',
            'pstart': [0,0,0], 
            'atoms': 20,
            'direction': [0,0,1],
            'bond type': 'B1',
            'angle type': 'A1'
        }
        self.armparams = {
            'bond length': [1],

            'narms': 20,
            'back to arm ids': [],
            'types': ['Rod'],
            'atoms': [20], 
            'directions': [[1,0,0]],
            'diff side': True,
            'even dist': True,
            'arm gap': 1,
            'bond types': ['B2'],
            'angle types': ['A2']
        }

        # print(self.backparams)
        # print(self.armparams)

        self.BranchedPolymerType = 'Linear'
        self.StarPolymer = False
        self.SP_r = 1.0
    
    def create_backbone(self):
        
        pas = self.backparams
        bl = pas['bond length']
        btyp = pas['type']
        atoms = pas['atoms']
        if self.StarPolymer:
            atoms=self.armparams['narms']
        bdir = pas['direction']
        bstart = pas['pstart']

        btyp = 'sphere' if self.StarPolymer==True else btyp

        btopo = None

        if btyp.lower() in ['rod', 'l', 'linear', 'line']:
            bc = Rod()
            btopo = bc.createRod(bstart, bdir, bl, atoms)

        elif btyp.lower() in ['c', 'c1', 'c2', 'circle', 'closed circle', 'open circle']:
            radius = atoms/(2*np.pi)
            bc = Ring()
            bc.atom_sid = 0
            if btyp.lower() in ['c2', 'open circle']:
                bc.closed = False
            plane = bdir
            if bdir==[0,0,1]:
                plane = 'yoz'
            elif bdir==[0,1,0]:
                plane = 'yoz'
            elif bdir==[1,0,0]:
                plane = 'xoz'
            btopo = bc.createRing(radius, bl, atoms, plane)

        elif btyp.lower() == 'sphere':
            r = np.sqrt(atoms/(4*np.pi))
            bc = Sphere()
            bc.haveCenter = True
            if self.SP_r<=1.0:
                bc.r = r+2
            else:
                bc.r = self.SP_r
            btopo = bc.createSphere(atoms)
            btopo.bonds = _arr([[0, i] for i in range(1, atoms+1)])

        return btopo

    def create_arms(self, backtopo):

        btoms = backtopo.atoms
        bcoords = backtopo.coords
        bbonds = backtopo.bonds
        bangles = backtopo.angles

        apas = self.armparams
        narms = apas['narms']

        bls = apas['bond length']
        if type(bls) == float:
            bls = [bls]*narms
        elif len(bls) < narms:
            bls = bls+[bls[-1]]*(narms-len(bls))

        back_ids = apas['back to arm ids']
        sids = back_ids
        if type(back_ids)==int or len(back_ids)==0:
            gap = apas['arm gap']
            if apas['even dist']:
                gap = int(btoms/narms)
            gap = 1 if gap==0 else gap
            sids = [i*gap for i in range(np.min([int(btoms/gap), narms]))]
        elif len(back_ids) != narms:
            raise ValueError('Wrong backbone atoms ids for arms!')
        sids = np.array(sids)
        if self.StarPolymer:
            sids = np.array([i for i in range(1,narms+1)])
        spoints = bcoords[sids]

        self.armparams.update({'back to arm ids':sids})

        atyps = apas['types']
        if type(atyps) == str:
            atyps = [atyps]*narms
        elif len(atyps)<narms:
            atyps = atyps+[atyps[-1]]*(narms-len(atyps))

        dirs = apas['directions']
        if type(dirs[0])==float or type(dirs[0])==int:
            dirs = [dirs]*narms
        elif len(dirs)<narms:
            dirs = dirs + [dirs[-1]]*(narms-len(dirs))
        if apas['diff side']:
            dirs = np.array(dirs)
            idirs = _arr([i for i in range(len(dirs))])
            dirs = [np.power(-1, i)*d for i,d in zip(idirs, dirs)]

        if self.StarPolymer:
            dirs = bcoords[1:,:]

        arms_atoms = apas['atoms']
        if type(arms_atoms)==int:
            arms_atoms = [arms_atoms]*narms
        elif len(arms_atoms)<narms:
            arms_atoms = arms_atoms + [arms_atoms[-1]]*(narms-len(arms_atoms))

        atopos = {}
        for i in range(narms):

            atyp = atyps[i]
            adir = np.array(dirs[i])
            adir = adir/np.linalg.norm(adir)
            bl = bls[i]
            atoms = arms_atoms[i]
            
            
            if self.StarPolymer:
                spoint = spoints[i]
                btoms = 1
            else:
                spoint = spoints[i] + bl*adir

            atopo = None
            if atyp.lower() in ['rod', 'l', 'linear', 'line']:
                ac = Rod()
                ac.atom_sid = btoms+i*atoms
                atopo = ac.createRod(spoint, adir, bl, atoms)
                atopo.connect_bonds = np.array([[sids[i], ac.atom_sid]])

            atopos.update({f'arm{i}':atopo})
            
        return atopos
        
    def combine(self, backTopo, armTopos):

        bT, aTs = backTopo, armTopos
        
        atoms_list = self.armparams['atoms']
        atoms_list = [atoms_list] if type(atoms_list)!=list else atoms_list
        if len(atoms_list)==0:
            self.topoInfo = bT
            return bT
        elif np.max(atoms_list)==0:
            self.topoInfo = bT
            return bT

        sids = self.armparams['back to arm ids']

        self.topoInfo.sub_polymers.update({'back':bT})
        self.topoInfo.coords = bT.coords
        self.topoInfo.bonds = bT.bonds
        self.topoInfo.angles = bT.angles
        self.topoInfo.sub_polymers_ids.update({'back':bT.atom_ids})

        self.topoInfo.atoms_typeid.update()
        self.topoInfo.bonds_typeid.update()
        self.topoInfo.angles_typeid.update()

        for key, val, sid in zip(aTs.keys(), aTs.values(), sids):
        
            self.topoInfo.sub_polymers.update({key:val})
            self.topoInfo.sub_polymers_ids.update({key:val.atom_ids})

            self.topoInfo.coords = np.concatenate([self.topoInfo.coords, val.coords], axis=0)

            inter_bond = _arr([[sid, val.atom_ids[0]]])
            if self.StarPolymer:
                inter_bond = _arr([[0, val.atom_ids[0]]])
            if len(self.topoInfo.bonds)==0:
                if len(val.bonds)==0:
                    self.topoInfo.bonds = np.array(val.bonds)
                else:
                    self.topoInfo.bonds = np.concatenate([inter_bond, val.bonds], axis=0)
            else:
                if len(val.bonds)==0:
                    self.topoInfo.bonds = np.concatenate([self.topoInfo.bonds, inter_bond], axis=0)
                else:
                    self.topoInfo.bonds = np.concatenate([self.topoInfo.bonds, inter_bond, val.bonds], axis=0)

            if len(val.atom_ids)<=1:
                continue
                
            inter_ags = _arr([[sid, val.atom_ids[0], val.atom_ids[1]]])
            if self.StarPolymer:
                inter_ags = _arr([[0, val.atom_ids[0], val.atom_ids[1]]])
            if len(self.topoInfo.angles) == 0:
                if len(val.angles)==0:
                    self.topoInfo.angles = np.array(inter_ags)
                else:
                    self.topoInfo.angles = np.concatenate([inter_ags, val.angles], axis=0)
            else:
                if len(val.angles)==0:
                    self.topoInfo.angles = np.concatenate([self.topoInfo.angles, inter_ags], axis=0)
                else:
                    self.topoInfo.angles = np.concatenate([self.topoInfo.angles, inter_ags,val.angles], axis=0)

        return self.topoInfo
   
    def createBranchedPolymer(self, dumpTopo=True, topoFname=None):

        btopo = self.create_backbone()
        atopos = self.create_arms(btopo)

        if self.StarPolymer:
            P = Point()
            btopo = P.createPoint(0,0,0)

        BTopo = self.combine(btopo, atopos)
        
        if dumpTopo:
            BTP = BranchedTopoParser()
            if topoFname is None:
                b_l, a_l = self.backparams['atoms'], self.armparams['narms']
                topoFname=f'b{b_l}a{a_l}.conf'
            BTP.dump(BTopo, topoFname)

        return self.topoInfo
        
if __name__ == '__main__':
    BP = BranchedPolymer()
    bpTopo = BP.createBranchedPolymer(dumpTopo=False)
    data = bpTopo.coords
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(data[:,0], data[:,1], data[:,2])
    plt.show()

    Sp = BranchedPolymer()
    Sp.StarPolymer=True
    Sp.armparams.update({'narms':30, 'atoms':[10]})
    spTopo = Sp.createBranchedPolymer(dumpTopo=False)
    data = spTopo.coords
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(data[:,0], data[:,1], data[:,2])
    plt.show()

