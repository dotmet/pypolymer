import numpy as np
import time
import copy

from topo_tools_cpp import parse_boundary, compute_gyration, compute_rcm


class TrjAnalysis(object):

    def __init__(self):
    
        self.topos = []
        self.base_topo = None
        self.baseTopo = None
        
        self.options = [
            'msd', 'msd_submol', 'gyration', 'Rg'
        ]
        self.results = {}
        self.compute_items = {}
        self.compute_items_params = {}
        self.no_cal_computes = {}
        
        self.custom_computes = {}
        self.custom_computes_params = {}
        self.custom_results = {}
        
        self._default_funcs()

    def loadtopo(self, topo=None, topoFile=None, topoParser=None):
        PTP = topoParser
        if topo is None:
            topo = PTP.load(topoFile)
        self.base_topo = topo
        bonds = self.base_topo.bonds
        self.base_topo.bonds = bonds[np.argsort(bonds[:,0])]
        return topo
        
    def loadlmp(self, lmpfile, topoParser=None, steps=[]):
    
        if type(steps)==int:
            steps = [steps]
        
        topos = []
        for step in steps:
            PTP = topoParser
            PTP.load(topo=copy.deepcopy(self.base_topo))
            topo = PTP.match_lmp(lmpfile=lmpfile, step=step, analysis=True)
            topos.append(topo)
        self.topos = topos
        return topos
            
    def loadgsd(self, gsdfile, topoParser=None, steps=[]):
    
        if type(steps)==int:
            steps = [steps]
        
        topos = []
        for step in steps:
            PTP = topoParser
            PTP.load(topo=copy.deepcopy(self.base_topo))
            PTP.match_gsd(gsdfile=gsdfile, step=step, analysis=True)
            topos.append(PTP.topo)
            
        self.topos = topos
        return topos
    
    ################
    def add_compute_func(self, name, func, params=[]):
        try:
            func(self.trj[0], self.trj[1], params)
        except:
            warn_message = '''Custom functions must set 'PolymerTopoInfo' at t=0 as 
                            first parameter, and at t=t1 as another, Other parameter 
                            should be stored in a list.'''
            raise Warning(warn_message)
        self.custom_computes.update({name:func})
        self.custom_computes_params.update({name:params})
        
    def remove_compute_func(self, name):
        self.custom_computes.pop(name)
        self.custom_computes_params.pop(name)
        
    def update_compute_func(self, name, params=[]):
        self.custom_computes_params.update({name:params})
    ################
    
    def compute1(self, func, params=[], useT0=True):
        results = []
        for topo in self.topos:
            if useT0:
                results.append(func(self.topos[0], topo, params))
            else:
                results.append(func(topo, params))
        return results
        
    def analysis(self):
    
        # Part 1
        # No more parsing needed, just collect data from static calculation.
        gyras = []
        R_gyra = []
        tan2XG = []
        rcm = []
        XG = []
        Span = []
        for topo in self.topos:
            gyras.append(topo.gyration.flatten())
            # _coords = topo.coords
            # _rcm = np.mean(_coords, axis=0)
            # c2rcm = _coords-_rcm
            # gyration = compute_gyration(c2rcm, c2rcm.shape[0])
            R_gyra.append(topo.Rg)
            tan2XG.append(topo.tan2XG)
            XG.append(topo.XG)
            rcm.append(topo.rcm)
            Span.append([topo.xw, topo.yw, topo.zw])
        self.results.update({'gyration':np.array(gyras)})  
        self.results.update({'Rg':np.array(R_gyra)})
        self.results.update({'tan2XG':tan2XG})
        self.results.update({'XG':XG})
        self.results.update({'rcm':rcm})
        self.results.update({'Configuration span': Span})
        
        gyras = np.array(gyras, dtype=float)
        self.results.update({'gxx':gyras[:,0]})
        self.results.update({'gyy':gyras[:,4]})
        self.results.update({'gzz':gyras[:,8]})

        # Collecting data of submol
        gyra_sub = {}
        R_gyra_sub = {}
        tan2XG_sub = {}
        XG_sub = {}
        torsional_angle_sub = {}
        rcm_sub = {}
        Span_sub = {}
        for topo in self.topos:
            for key,val in topo.sub_polymers.items():
                if len(gyra_sub)<len(topo.sub_polymers):
                    gyra_sub.update({key:[val.gyration.flatten()]})
                    R_gyra_sub.update({key: [val.Rg]})
                    tan2XG_sub.update({key: [val.tan2XG]})
                    XG_sub.update({key:[val.XG]})
                    torsional_angle_sub.update({key:[topo.torsional_angle_submol[key]]})
                    rcm_sub.update({key:[val.rcm]})
                    Span_sub.update({key:[[val.xw, val.yw, val.zw]]})
                else:
                    gyra_sub[key] += [val.gyration.flatten()]
                    R_gyra_sub[key] += [val.Rg]
                    tan2XG_sub[key] += [val.tan2XG]
                    XG_sub[key] += [val.XG]
                    torsional_angle_sub[key] += [topo.torsional_angle_submol[key]]
                    rcm_sub[key] += [val.rcm]
                    Span_sub[key] += [[val.xw, val.yw, val.zw]]

        for key,val in self.topos[0].sub_polymers.items():
            gyra_sub[key] = np.array(gyra_sub[key])
            R_gyra_sub[key] = np.array(R_gyra_sub[key])
            tan2XG_sub[key] = np.array(tan2XG_sub[key])
            XG_sub[key] = np.array(XG_sub[key])
            torsional_angle_sub[key] = np.array(torsional_angle_sub[key])
            rcm_sub[key] = np.array(rcm_sub[key])
            Span_sub[key] = np.array(Span_sub[key])

        self.results.update({'gyration_submol':gyra_sub})
        self.results.update({'Rg_submol': R_gyra_sub})
        self.results.update({'tan2XG_submol': tan2XG_sub})
        self.results.update({'XG_submol':XG_sub})
        self.results.update({'torsional_angles':torsional_angle_sub})
        self.results.update({'rcm_sub':rcm_sub})
        self.results.update({'Configuration span submol':Span_sub})

        # Part 2
        # Properties need to calculate from t=0 and t=t time step
        for item,func in self.compute_items.items():
            results = []

            for i in range(len(self.topos)):
                topo = self.topos[i]
                params = []
                if item in self.compute_items_params:
                    params = self.compute_items_params[item]
                res = func(self.topos[0], topo, params=params)
                if type(res) == dict:
                    if i==0:
                        results = res
                    else:
                        for key,val in res.items():
                            results[key] = np.vstack([results[key], val])
                else:
                    results.append(res)
            if type(results)==dict:
                for key,val in results.items():
                    results[key] = np.squeeze(val)
                self.results.update({item: results})
            else:
                self.results.update({item: np.array(results)})
        
    def _default_funcs(self):
        
        # Calculate Mean Square displacement for time t
        def MSD_self(Tp0, Tpt, params=[]):
            return np.dot(Tpt.rcm-Tp0.rcm, Tpt.rcm-Tp0.rcm)
        self.compute_items.update({'msd':MSD_self})
        self.compute_items_params.update({'msd':[]})
        
        def MSD_submol(Tp0, Tpt, params=[]):
            msd_arms = []
            for poly0, polyt in zip(Tp0.sub_polymers.values(), Tpt.sub_polymers.values()):
                msd_arms.append(np.dot(polyt.rcm-poly0.rcm, polyt.rcm-poly0.rcm))
            return np.mean(msd_arms)
        self.compute_items.update({'msd_submol':MSD_submol})
        self.compute_items_params.update({'msd_submol':[]})

        def self_correlation(Tp0, Tpt, params=[]):
            return np.dot(Tp0.rcm, Tpt.rcm)/np.dot(Tp0.rcm, Tp0.rcm)
        self.compute_items.update({'correlation':self_correlation})
        self.compute_items_params.update({'correlation':[]})

        def self_correlation_submol(Tp0, Tpt, params=[{'sid':1}]):
            crsup = []
            crsdown = []
            for poly0, polyt in zip(Tp0.sub_polymers.values(), Tpt.sub_polymers.values()):
                crsup.append(np.dot(poly0.ete_vec, polyt.ete_vec))
                crsdown.append(np.dot(poly0.ete_vec, poly0.ete_vec))
            return np.mean(crsup[params[0]['sid']:])/np.mean(crsdown[params[0]['sid']:])
        self.compute_items.update({'correlation_submol':self_correlation_submol})
        self.compute_items_params.update({'correlation_submol':[{'sid':1}]})

        def self_correlation_submol_1(Tp0, Tpt, params=[{'sid':1}]):
            crsup1 = []
            crsup2 = []
            crsdown1 = []
            crsdown2 = []
            for poly0, polyt in zip(Tp0.sub_polymers.values(), Tpt.sub_polymers.values()):
                e0,et = poly0.ete_vec, polyt.ete_vec
                crsup1.append(np.dot(e0, et))
                crsup2.append(np.linalg.norm(e0))
                crsdown1.append(np.dot(e0, e0))
                crsdown2.append(np.linalg.norm(e0))
            up = np.mean(crsup1[params[0]['sid']:])-np.power(np.mean(crsup2[params[0]['sid']:]),2)
            dwn = np.mean(crsdown1[params[0]['sid']:])-np.power(np.mean(crsdown2[params[0]['sid']:]),2)
            return up/dwn
        self.compute_items.update({'correlation_submol_1':self_correlation_submol_1})
        self.compute_items_params.update({'correlation_submol_1':[{'sid':1}]})

        def TACF(Tp0, Tpt, params=[]):
            phits = Tpt.torsional_angle_submol
            phi0s = Tp0.torsional_angle_submol
            taup1 = []
            taup2 = []
            tadown1 = []
            tadown2 = []
            for phit, phi0 in zip(phits.values(), phi0s.values()):
                if phi0==10086:
                    continue
                else:
                    taup1.append(np.cos(phit)*np.cos(phi0))
                    taup2.append(np.cos(phi0))
                    tadown1.append(np.cos(phi0)*np.cos(phi0))
                    tadown2.append(np.cos(phi0))
            taup = np.mean(taup1)-np.power(np.mean(taup2),2)
            tadn = np.mean(tadown1)-np.power(np.mean(tadown2),2)
            return taup/tadn
        self.compute_items.update({'TACF':TACF})
        self.compute_items_params.update({'TACF':[]})

        def velo_correlation_submol(Tp0, Tpt, params=[]):
            vcrup = []
            vcrdown = []
            for poly0, polyt in zip(Tp0.sub_polymers.values(), Tpt.sub_polymers.values()):
                pass
        
        def cross_correlation(Tp0, Tpt, params=[]):
            pass

        def stress_tensor():
            pass
            
    def structure_analysis(self, topo=None):
    
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
        # topo = self.polymer_torsional(topo=topo)
        
        # Compute span along x,y,z
        xcrds, ycrds, zcrds = coords[:,0], coords[:,1], coords[:,2]
        topo.xw = np.max(xcrds) - np.min(xcrds)
        topo.yw = np.max(ycrds) - np.min(ycrds)
        topo.zw = np.max(zcrds) - np.min(zcrds)
        
        
        return topo
        
    def sub_structure_analysis(self, topo=None):
        
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
                            raise Warning('The bonds near connection bead are less than 3, please check')
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
        