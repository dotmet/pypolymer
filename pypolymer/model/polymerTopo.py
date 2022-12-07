import pickle
import os
import copy
import time
import gsd
import gsd.hoomd

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from .topoTools import mol_cross_boundary
from pypolymer.analyze.lammps.trjdata import TrjData
from pypolymer.util import parse_boundary, compute_gyration, compute_rcm

class PolymerTopoInfo(object):

    def __init__(self):

        self.name = 'Polymer'
        self.idx = 0
        self.atom_ids = [0]
        self.bounds_info = {'xlo':0, 'xhi':0, 'ylo':0, 'yhi':0, 'zlo':0, 'zhi':0}
        self.bounds_rule = []
        
        # If this polymer have sub polymers.
        self.IS_BASE_POLYMER = True

        # Basic properties
        self.atoms = 1
        self.coords = [[0,0,0]] # dynamic properties
        self.velocity = [[0,0,0]]
        self.raw_coords = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        self.Masses = [1.0]
        self.box = []
        self.connect_bonds = []

        # Typeid
        self.atoms_typeid = {}
        self.bonds_typeid = {}
        self.angles_typeid = {}
        self.dihedrals_typeid = {}

        # Sub structures
        '''
            These attributes set for specific sub classes.
            In example:
                BranchedTopoInfo(PolymerTopoInfo)
        '''
        self.sub_polymers = {}
        self.sub_polymers_ids = {}

        # Structural properties at a specific time step
        self.rcm = [0, 0, 0]

        self.gyration = [[0, 0, 0],
                        [0, 0, 0],
                        [0, 0, 0]]
        self.Rg = 0
        self.tan2XG = 0
        self.XG = 0
        
        self.xw = 0 # Polymer configuration span along x
        self.yw = 0 # Polymer configuration span along y
        self.zw = 0 # Polymer configuration span along z

        self.torsional_angle_submol = {}
        
        self.ete_vec = []


class PolymerTopoParser(object):

    def __init__(self):
    
        self.topo = None
        self.bounds_info = {'xlo':0, 'xhi':0, 'ylo':0, 'yhi':0, 'zlo':0, 'zhi':0}
        self.bounds_rule = []
        self.PARSE_BOUND = True
        self.box = []

    def init_TopoInfo(self):
        return PolymerTopoInfo()

    def dump(self, topo=None, fname=None):

        if fname is None:
            fname = 'branched.conf'
        if topo is None:
            topo = self.topo
        
        with open(fname, 'wb') as f:
            pickle.dump(topo, f)

    def load(self, fname=None, topo=None):

        try:
            os.path.exists(fname)
            with open(fname, 'rb') as f:
                topo = pickle.load(f)
        except:
            if topo is None:
                raise ValueError('Wrong configuration file specified!')
        
        self.topo = topo
        
        return topo

    ##
    def match(self, topo=None, data=[], velo=[], analysis=False):
    
        if topo is None:
            topo = self.topo
        
        topo.coords = data
        topo.velocity = velo

        for key,val in topo.sub_polymers_ids.items():
        
            coords = data[np.array(val, dtype=int)]
            velos = velo[np.array(val, dtype=int)]
            topo.sub_polymers[key].coords = coords
            topo.sub_polymers[key].velocity = velos
        
        if analysis:

            topo = self.structure_analysis(topo=topo)
            topo = self.sub_structure_analysis(topo=topo)
            
        self.topo=topo
        return topo
    
    def match_lmp(self, lmpfile=None, topo=None, step=0, analysis=False):
    
        trj = TrjData(lmpfile, step=step, MUTE=True)
        
        data = trj.coords
        self.bounds_info = trj.Header.bounds_info
        self.bounds_rule = trj.Header.bounds_rule

        return self.match(topo=topo, data=data, analysis=analysis)
        
    def match_gsd(self, gsdfile=None, topo=None, step=0, analysis=False):
    
        snaps = gsd.hoomd.open(gsdfile)
        
        data = snaps[step].particles.position
        velo = snaps[step].particles.velocity
        self.box = snaps[step].configuration.box[:3]
        bonds = snaps[step].bonds.group[:]

        if self.PARSE_BOUND:
            data = parse_boundary(data, bonds, self.box, data.shape[0], bonds.shape[0])

        return self.match(topo=topo, data=data, velo=velo, analysis=analysis)

    def structure_analysis(self, topo=None):
    
        if topo is None:
            topo = self.topo
        if len(topo.coords)==0:
            raise RuntimeError('No topological info found to analysis!')
        
        # compute center of mass
        coords = topo.coords
        mass = topo.Masses
        rcm = np.mean(coords, axis=0)
        topo.rcm = rcm
        
        # compute gyration
        r2cm = coords-rcm
        gyra = compute_gyration(r2cm, r2cm.shape[0])
        topo.gyration = gyra
        topo.Rg = np.sum([gyra[0,0], gyra[1,1], gyra[2,2]])
        topo.tan2XG = 2*gyra[0,2]/(gyra[0,0]-gyra[2,2])
        topo.XG = 0.5*np.arctan(topo.tan2XG)
        
        # compute end-to-end vector
        topo.ete_vec = coords[-1]-coords[0]

        # compute torsional angle that connection bond
        # make with near bonds on backbone
        topo = self.polymer_torsional(topo=topo)
        
        # Compute span along x,y,z
        xcrds, ycrds, zcrds = coords[:,0], coords[:,1], coords[:,2]
        topo.xw = np.max(xcrds) - np.min(xcrds)
        topo.yw = np.max(ycrds) - np.min(ycrds)
        topo.zw = np.max(zcrds) - np.min(zcrds)
        
        self.topo = topo
        
        return topo
        
    def sub_structure_analysis(self, topo=None):
        
        if topo is None:
            topo = self.topo
        if len(topo.coords)==0:
            raise RuntimeError('No topological info found to analysis!')
    
        if len(topo.sub_polymers)>=1:
        
            for pname, polymer in topo.sub_polymers.items():

                _polymer = self.structure_analysis(polymer)
                topo.sub_polymers.update({pname:_polymer})
                
        return topo

    def polymer_torsional(self, topo):

        if len(topo.sub_polymers)==1:
            return None
            
        mbonds = None
        for key, val in topo.sub_polymers.items():
            torsions = []
            if mbonds is None:
                mbonds = val.bonds
                mcoords = val.coords
                torsions = [10086]
            else:
                cbonds = val.connect_bonds
                for idx in np.unique(mbonds):
                    res = np.where(cbonds==idx)
                    if len(res[0])==0:
                        continue
                    else:
                        subids = self.find_bonded_atoms(idx, cbonds)[0]
                        # ptor1 = topo.coords[subids]
                        ptor1 = val.coords[-1]

                        ptor0 = mcoords[idx]

                        mids = self.find_bonded_atoms(idx, mbonds)
                        if len(mids)<4:
                            torsions = [10086]
                            # raise Warning('The bonds near connection bead are less than 3, please check')
                        else:
                            p1 = mcoords[mids[0]]
                            p2 = mcoords[mids[-1]]
                            torsions.append(self.compute_torsional_angle(p1,p2,ptor0,ptor1))
            topo.torsional_angle_submol.update({key:np.mean(torsions)}) 
            
        return topo      

    def find_bonded_atoms(self, atomid, bonds):
        
        atomids = []
        res = np.where(bonds==atomid)
        for idx1 in res[0]:
            for idx2 in res[1]:
                idx2 = 1 if idx2==0 else 0
                atomids.append(bonds[idx1, idx2])
        
        return np.array(atomids)

    def compute_torsional_angle(self, p1, p2, ptor0, ptor1):

        v1 = p1-ptor0
        v2 = p2-ptor0
        vtor = ptor1-ptor0

        nvec = np.cross(v1, v2)
        cosval = np.dot(nvec, vtor)/(np.linalg.norm(nvec)*np.linalg.norm(vtor))
        
        return np.pi/2-np.arccos(np.abs(cosval))

    def compute_rcm(self, coords, mass):
    
        if type(mass)==float or type(mass)==int:
            mass = [mass]*coords.shape[0]
        if len(mass)<coords.shape[0]:
            mass = mass + [1.0]*(coords.shape[0]-len(mass))
        
        mcoords = coords*np.array([mass]).T
        rcm = (1/np.sum(mass))*compute_rcm(mcoords, mcoords.shape[0])
        
        return rcm

    def parseBoundary(self, coords, bonds, box=[]):

        '''
            Parse periodic boundary conditions to keep molecular integrity. 
        '''
        if len(box) == 0:
            bf = self.bounds_info
            box.append(bf['xhi']-bf['xlo'])
            box.append(bf['yhi']-bf['ylo'])
            box.append(bf['zhi']-bf['zlo'])
            self.box = box

        box = self.box

        return parse_boundary(coords, bonds, box, coords.shape[0], bonds.shape[0])

    def image(self, topo=None, plane='xoz', dim3=False, show=True):
    
        if topo is None:
            topo = self.topo
        if dim3 or plane.lower()=='xyz':
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            i = 0
            for top in topo.sub_polymers.values():
                coords = topo.coords
                ax.scatter(coords[:,0], coords[:,1], coords[:,2])
                if i==0:
                    ax.plot(coords[:0], coords[:,1], coords[:,2])
                i = i+1
            if show:
                fig.show()
            else:
                return fig, ax
        else:
        
            axis1, axis2 = 0, 2
            if plane.lower() in ['xoz', 'xz']:
                axis1, axis2 = 0, 2
            elif plane.lower() in ['yoz', 'yz']:
                axis1, axis2 = 1, 2
            elif plane.lower() in ['xoy', 'xy']:
                axis1, axis2 = 0, 1
                
            fig, ax = plt.subplots(1,1)
            
            s_topos = topo.sub_polymers.values()

            Ncolors = len(s_topos)
            colormap = plt.cm.viridis# LinearSegmentedColormap
            Ncolors = min(colormap.N,Ncolors)
            mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]

            i = 0
            for top, c in zip(s_topos, mapcolors):
                cds = top.coords
                ax.plot(cds[:,axis1], cds[:,axis2], marker='o', c=c, markersize=5, lw=1, markerfacecolor='w',alpha=0.8)
                if i==0:
                    ax.plot(cds[:,axis1], cds[:,axis2], color='red', lw=3, zorder=Ncolors+1)
                    ax.plot(cds[[0,-1],axis1], cds[[0,-1],axis2], color='red', lw=3, zorder=Ncolors+1)
                i+=1
                
            if show:
                plt.show()
            else:
                return fig, ax
                
    def kde(self, topo=None, plane='xoz', submol=None):
    
        if topo is None:
            topo = self.topo
            
        axis1, axis2 = 0, 2
        if plane.lower() in ['xoz', 'xz']:
            axis1, axis2 = 0, 2
        elif plane.lower() in ['yoz', 'yz']:
            axis1, axis2 = 1, 2
        elif plane.lower() in ['xoy', 'xy']:
            axis1, axis2 = 0, 1
        
        coords = topo.coords
        
        if submol is None:
            pass
        else:
            coords = topo.sub_polymers[submol].coords
            
        fig, ax = plt.subplots(1,1)
        cmap = sns.cubehelix_palette(start=1.8, light=1, as_cmap=True)
        sns.kdeplot(x=coords[:,axis1], 
            y= coords[:,axis2], 
            levels=20, 
            cmap = cmap,
            thresh=.1,
            ax = ax,
            fill=True)
        
        return fig, ax