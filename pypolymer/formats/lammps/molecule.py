import numpy as np
from .utilTools import *

class Molecule(object):
    
    def __init__(self, idx=0, exist_attribute=True):
        
        self.attriDict={
                        "atoms":0, "bonds":0, "angles":0, "dihedrals":0, 
                        "impropers":0, "block":'\n', "atom types":0, "bond types":0, 
                        "angle types":0, "dihedral types":0, "improper types":0
                        }
        
        self.dataPart={
                        "Masses":'', "Nonbond Coeffs":'', "Bond Coeffs":'',
                       "Angle Coeffs":'', "Dihedral Coeffs":'', "Atoms":'',  
                       "Bonds":'',  "Angles":'',  "Dihedrals":'', "Impropers":''
                        }
        
        self.dataGroup = [
                            ['Masses'], ["Nonbond Coeffs", "Bond Coeffs",
                           "Angle Coeffs", "Dihedral Coeffs"], ["Atoms"], 
                          ["Bonds",  "Angles",  "Dihedrals", "Impropers"]
                            ]
        self.id = idx
        self.atom_id = 0
        self.bond_id = 0
        self.angle_id = 0
        self.dihedral_id = 0

        self.coords = [] ####
        
        self.atoms = []
        self.bonds = []
        self.angles = []
        self.dihedrals = []
        
        self._radium = 0
        
        self.atom_type = 1
        self.bond_type = 1
        self.angle_type = 1
        self.dihedral_type = 1
        
        self.exist_attribute = exist_attribute
        self.prevM = 0
        self.CHANGE_SERIAL = True
        self.HAVE_INDEX = True
        
    # Special initialize function
    def read_data(self, Prev_M=None, Masses=[], Coords=[], Bonds=[], Angles=[], Dihedrals=[], atyp=1, btyp=1, agtyp=1, dtyp=1):
        
        dic = self.dataPart
        ############## Set Params ###############
        self.prevM = Prev_M
        self.coords = Coords
        self.bonds = Bonds
        self.angles = Angles
        self.dihedrals = Dihedrals
        self.atom_type = atyp
        self.bond_type = btyp
        self.angle_type = agtyp
        self.dihedral_type = dtyp
        ############## End ######################
        ############## Change serial ############
        _id = 0 if (Prev_M is None) else Prev_M.atom_id
        if self.CHANGE_SERIAL:
            Bonds = change_serial(Bonds, _id)
            Angles = change_serial(Angles, _id)
            Dihedrals = change_serial(Dihedrals, _id)
        else:
            Bonds = np.array(Bonds, dtype=object)
            Angles = np.array(Angles, dtype=object)
            Dihedrals = np.array(Dihedrals, dtype=float)
        ############## End ################
        
        self.atom_id = 0 if (Prev_M is None) else Prev_M.atom_id
        self.bond_id = 0 if (Prev_M is None) else Prev_M.bond_id
        self.angle_id = 0 if (Prev_M is None) else Prev_M.angle_id
        self.dihedral_id = 0 if (Prev_M is None) else Prev_M.dihedral_id
        atoms = len(Coords)
        bonds = len(Bonds)
        angles = len(Angles)
        dhls = len(Dihedrals)
        
        self.dataPart['Masses'] = Masses
        self.dataPart['Atoms'] = Coords
        self.dataPart['Bonds'] = Bonds
        self.dataPart['Angles'] = Angles
        self.dataPart['Dihedrals'] = Dihedrals
        
        kd = self.dataGroup[2]+self.dataGroup[3]
        for key in kd:
            self.attriDict[key.lower()] = len(self.dataPart[key])
#         self.attriDict['atom types'] = atyp
#         self.attriDict['bond types'] = btyp
#         self.attriDict['angle types'] = agtyp
#         self.attriDict['dihedral types'] = dtyp
            
        for key,v in self.dataPart.items():
            if key in self.dataGroup[0]:
                dic[key]=self._MASS(self.dataPart[key])
            if key in self.dataGroup[1] and len(v)>0:
                dic[key]=self._COEF(self.dataPart[key])
            if key in self.dataGroup[2]:
                start_id = self.atom_id+1
                dic[key]=self._ATOM(start_id, self.dataPart[key], atyp)
            if key == 'Bonds' and len(v)>0:
                start_id = self.bond_id + 1
                dic[key] = self._BADI(start_id, self.dataPart[key], btyp)
            if key == 'Angles' and len(v)>0:
                start_id = self.angle_id + 1
                dic[key] = self._BADI(start_id, self.dataPart[key], agtyp)
            if key == 'Dihedrals' and len(v)>0:
                start_id = self.dihedral_id + 1
                dic[key] = self._BADI(start_id, self.dataPart[key], dtyp)
        self.dataPart = dic
        self.atoms = self.attriDict['atoms']
        
        self.atom_id = self.atom_id + atoms
        self.bond_id = self.bond_id + bonds
        self.angle_id = self.angle_id + angles
        self.dihedral_id = self.dihedral_id + dhls
        
        return self
    
    def correct_id(self, prev_m, atoms, bonds=[], angles=[], dihedrals=[]):
        patoms = prev_m.atoms
        pbonds = prev_m.bonds
        pangles = prev_m.angles
        pdihedrals= prev_m.dihedrals
        try:
            b1 = np.linalg.norm(atoms[1]-atoms[0])
        except:
            b1 = 1
        try:
            b2 = np.linalg.norm(patoms[1]-patoms[0])
        except:
            b2=1
        for i in range(atoms.shape[0]):
            for j in range(patoms.shape[0]):
                if np.linalg.norm(atoms[i, 2:5]-patoms[j, 2:5])<np.min([b1, b2])/2:
                    atoms = np.delete(atoms, i, axis=0)
                    _b = bonds[:,2:]
                    _a = angles[:,2:]
                    _d = dihedrals[:,2:]
                    bonds[:,2:] = np.where(_b==atoms[i,0], patoms[j,0], atoms[i,0])
                    angles[:,2:] = np.where(_a==atoms[i,0], patoms[j,0], atoms[i,0])
                    dihedrals[:,2:] = np.where(_d==atoms[i,0], patoms[j,0], atoms[i,0])
        return atoms, bonds, angles, dihedrals
    # Print all informations of current molecule
    def printM(self):
        for k,v in self.attriDict.items():
            if v and k!='block':
                print(k + ':\n '+str(v))
        for k,v in self.dataPart.items():
            if v:
                print(k + ": \n"+v)
    
    def _setParameters(self):
        print("Set attribute of current molecular:")
        for key in self.attriDict:
            if key != "block":
                self.attriDict[key] = int(input(key+" :"))
    
    # MASS part, which can be defined in 'in' file, of the lammps data file
    def _MASS(self, data):
        return ""
    
    # Coefficient of potential function
    def _COEF(self, data):
        return ''
    
    # Parse current molecule's trajectory data
    # 'Atoms' part of lammps' data file
    def _ATOM(self, start_id, data, atom_type=1):
        num=len(data)
        if num==0:
            return ''
        index = self._gen1DArray(start_id, num, dtype=object)
        posi= np.around(np.array(data, dtype=float), decimals=5)
        atom_type = self._toArray([self.atom_type] * num, dtype=object)
        m_id = self._toArray([self.id]*num, dtype=object)
        b_array = np.concatenate([index, m_id, atom_type, posi], axis=1)
        return self._arraytoStr(b_array.tolist())
    
    # Data of bonds, angles, and dihedrals of current molecule
    # 'Bonds', 'Angles', and 'Dihedrals' parts of lammps file
    def _BADI(self, start_id, data, typ=1):
        num = len(data)
        if num==0:
            return ''
        typ_ = self._toArray([typ]*(num))
        index = self._gen1DArray(start_id, num)
        data = data[:,1:] if self.HAVE_INDEX else data
        b_array = np.concatenate([index, typ_, data], axis=1)
        return self._arraytoStr(b_array.tolist())

    # A function used to generate a specified 1D vector
    def _gen1DArray(self, start, num, dtype=int):
        index = [i for i in range(start, start+num)]
        return self._toArray(index, dtype=dtype)
    
    # A function used to convert 1D list to 1*n matrix
    def _toArray(self, l, dtype=int):
        l = np.array(l, dtype=dtype)
        return l.reshape(l.shape[0],1)
    
    # A function used to convert a matrix(a multidimentional numpy array)
    # to a string which can be used to output directly 
    def _arraytoStr(self, l):
        t=[]
        for i in l:
            s=list(map(lambda x:str(x), i))
            t.append("\t".join(s))
        return "\n".join(t)
