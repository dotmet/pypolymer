import numpy as np
import matplotlib.pyplot as plt
from .molecule import Molecule

class LmpDataFile(object):
    
    def __init__(self, box=(100,100,100), molecules=[], atom_typs=0, bond_typs=0, angle_typs=0, dihedral_typs=0, fileName="data.dat"):
    
        self.molecules = molecules
        
        self.fileName = fileName
        self.fileData = ""
        self.header = ''
        self.paramDict={"atoms":0, "bonds":0, "angles":0, "dihedrals":0, 
                        "impropers":0, "block":'\n', "atom types":0, "bond types":0, 
                        "angle types":0, "dihedral types":0, "improper types":0}
        self.filePart = {"Header":'', "Masses":'', "Nonbond Coeffs":'', "Bond Coeffs":'', 
                         "Angle Coeffs":'', "Dihedral Coeffs":'', "Atoms":'',  "Bonds":'', 
                         "Angles":'',  "Dihedrals":''}
        self.box = box
        self.boxRange = {"xlo":-100, "xhi":100, "ylo":-100, "yhi":100, "zlo":-100,"zhi":100}
        self.at = atom_typs
        self.bt = bond_typs
        self.agt = angle_typs
        self.dt = dihedral_typs
        
        self.oldVersion = False
        
        self._parse_box()
    
    def _parse_box(self):
        
        box = self.box
        if box is None:
            pass
        else:
            if len(box)==3:
                box = [-box[0]/2, box[0]/2, -box[1]/2, box[1]/2, -box[2]/2, box[2]/2]
                
            for i,key in enumerate(self.boxRange.keys()):
                self.boxRange[key] = box[i]
    
    # A function used to construct system
    def Merge(self, mlist):
        for mo in mlist:
            for k,v in mo.dataPart.items():
                if v:
                    self.filePart[k] += v+'\n'
            for k,v in mo.attriDict.items():
                if v and v!=1:
                    self.paramDict[k] += v
                if v==1:
                    self.paramDict[k] =1
        self.paramDict['atom types'] = self.at
        self.paramDict['bond types'] = self.bt
        self.paramDict['angle types'] = self.agt
        self.paramDict['dihedral types'] = self.dt
        self.paramDict['block'] = '\n'
    
    # Generate lammps data file
    def genFile(self):
        self.Merge(self.molecules)
        self.fileData = self._genHeader(self.paramDict)
        for k,v in self.filePart.items():
            if k == 'Atoms':
                self.fileData += '\n'+k+' #ID mol type x y z\n\n'+v
            elif v and k!='Masses':
                self.fileData += '\n'+k+' #id type atom1 atom2 ...\n\n'+v
        self._saveFile()
        
    # A function used to get all molecules' coordinates
    def get_posi(self):
        return self.__stack_attr('posi')
    
    def get_bonds(self):
        return self.__stack_attr('bonds')
    
    def get_angles(self):
        return self.__stack_attr('angles')
    
    def get_dihedrals(self):
        return self.__stack_attr('dihedrals')
                
    def __stack_attr(self, item='posi'):
        
        data = None
        for mol in self.molecules:
            if item=='posi':
                _it = mol.coords
            elif item=='bonds':
                _it = mol.bonds
            elif item=='angles':
                _it = mol.angles
            elif item=='dihedrals':
                _it = mol.dihedrals
                
            if data is None:
                data = _it
            else:
                data = np.vstack([data, _it])
        
        return data
    
    def setParam(self):
        self.fileName = input("Input file name:")
        print("Input values:")
        for key in self.paramDict:
            if key != "block":
                self.paramDict[key] = int(input(key+" :"))
        print("Set box range:")
        for key in self.boxRange:
            self.boxRange[key] = float(input(key+" :"))
        print("Set trajectory equation:")

    def Visualize(self, mlist=None):
    
        if mlist is None:
            mlist = self.molecules
    
        if isinstance(mlist, Molecule):
            mlist = [mlist]

        plt.subplot(projection='3d')
        for mol in mlist:
            p = np.array(mol.coord)
            plt.plot(p[:, 0], p[:, 1], p[:, 2], marker='o')
        plt.show()
        
    def _genHeader(self, dic):

        self._header = "#"+ self.fileName + "\n\n"
        for key,value in dic.items():
            if value == 0:
                pass
            elif key=="block":
                self._header += value
            else:
                self._header += "\t" + str(value) + '\t' + key + '\n'
        self._header += "\n"
        for i in range(0,5,2):
            l=list(self.boxRange)
            self._addRangeToHead(l[i],l[i+1])
#         self._header += "\nMasses\n\n" + "1 1.0\n" 
        return self._header
    
    def _addRangeToHead(self,k1, k2):
        br = self.boxRange
        if not self.oldVersion:
            self._header += '\t' + str(br[k1]) + "\t" + str(br[k2]) + "\t"+k1+'\t'+k2+'\n'
        else:
            self._header += str(br[k1]) + " " + str(br[k2]) + " "+k1+" "+k2+'\n'
        
    def _saveFile(self):
    
        with open(self.fileName,'w') as f:
            f.write(self.fileData)
            f.close()