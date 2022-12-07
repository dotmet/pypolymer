
import numpy as np
import matplotlib.pyplot as plt

from pylmp.Data.rawdata import RawData
from pylmp.Data.LmpTrj.util import *

class TrjHeader(object):
    
    def __init__(self, header):
        self.time_step = 0
        self.bounds_info = {'xlo':0, 'xhi':0, 'ylo':0, 'yhi':0, 'zlo':0, 'zhi':0}
        self.bounds_rule = []
        self.atoms = 0
        self.header = header
        self._header_parse()
        
    def _header_parse(self):
        h = self.header
        for i in range(len(h)):
            l = h[i]
            if type(h[i])!=list:
                l = l.split()
            if i==1:
                self.time_step = int(l[0])
            elif i==3:
                self.atoms = int(l[0])
            elif i==4:
                self.bounds_rule = l[3:6]
            elif i==5:
                self.bounds_info['xlo'] = float(l[0])
                self.bounds_info['xhi'] = float(l[1])
            elif i==6:
                self.bounds_info['ylo'] = float(l[0])
                self.bounds_info['yhi'] = float(l[1])
            elif i==7:
                self.bounds_info['zlo'] = float(l[0])
                self.bounds_info['zhi'] = float(l[1])
                
###
class TrjData(RawData):
    
    def __init__(self, file_name='', step=None, LAST_STEP=False, Print_detail=False, MUTE=False):
        super().__init__()
        self.total_atoms = 1
        self.total_steps = 1
        self.file_name = file_name
        self.axis = ['x', 'y', 'z']
        self.step_length = 1
        self.time_step = 1
        self.relax_step = 0
        self.time_line_list = []
        self.one_matrix = []
        self.coords_matrix = []
        self.coords_all_step = []
        self.data = []
        
        self.coords = [] # Only for specific step or steps when initially specified
        self.Header = [] # Trajectory header for one step or steps 
        self.velocities = []
        
        self.PRINT_DETAIL = Print_detail
        if MUTE:
            self.PRINT_DETAIL = False
        
        try:
            self.total_atoms = self.get_atoms()
        except:
            raise ValueError('Invalid input file')
            
        if step in [-1]:
            LAST_STEP=True
            
        if step is None or step in ['all', 'All', 'AlL', 'ALL', -1]:
                self.line_list = self.read_file(file_name)
                if LAST_STEP:
                    lmat = self.line_list[-self.total_atoms:]
                    self.f_data_only_coords(data_mat=lmat)
                    lmat = [i.split() for i in lmat]
                    lcoords = np.array(lmat, dtype=float)
                    lcoords = lcoords[np.argsort(lcoords[:,0])]
                    self.coords = lcoords[:,2:5]
                    self.Header = TrjHeader(header=self.line_list[-self.total_atoms-9:])
                    if not MUTE:
                        print(f'Only xyz component of positions for step {step} are stored in "coords" attribute of TrjData')
                else:
                    self.f_data_only_coords()
                self.raw_data = self.get_raw_data(file_name)
                self.coords_matrix = self.get_coords_matrix(self.data_only_coords)
                self.coords_all_step = self.coords_matrix
        else:
            steps = [step] if type(step)==int else step
            mats = []
            velo_mats = []
            Headers = []
            for step in steps:
                raw_mat = []
                header = []
                atoms = self.total_atoms
                span = atoms+9
                max_loop = span*(np.max(step)+1)
                with open(file_name, 'r') as f:
                    for i in range(max_loop):
                        l = f.readline()
                        if i>=span*step+9 and i<span*(step+1):
                            raw_mat.append(l.split())
                        elif i>=span*step and i<span*step+9:
                            header.append(l.split())
                raw_mat = np.array(raw_mat, dtype=float)
                raw_mat = raw_mat[np.argsort(raw_mat[:,0])]
                mats.append(raw_mat[:,2:5])
                Headers.append(TrjHeader(header))
            if len(steps) == 1:
                self.coords = mats[0]
                self.Header = Headers[0]
            else:
                self.coords = mats
                self.Header = Headers
            if not MUTE:
                print(f'Only xyz component of positions for step {step} are stored in "coords" attribute of TrjData')  
                if len(steps)>1:
                    print(f'Multiple steps ("step" parameter is a list) choosed, then the "coords" is a list.') 
        if not MUTE:
            print("ALL DATA HAS BEEN SORTED BY ATOMS' ID !!!")

    def get_data(self, step):
        data = self.data
        natoms = int(data[3].split()[0])
        span = natoms + 9
        data_mat = []
        data_raw = []
        if step<-1:
            data_raw = data[step*span:(step+1)*span]
        elif step==-1:
            data_raw = data[-natoms:]
        else:
            try:
                data_raw = data[step*span:(step+1)*span]
            except:
                raise RuntimeError('Wrong time step choosed!')
        data_mat = np.array([[x for x in l.split()] for l in data_raw][9:])
        return np.array(data_mat, dtype=float)

    # eliminate lines which are not info of atoms (str data)
    def f_data_only_coords(self, data_mat=[]):

        new_mat = []
        data_mat = self.line_list if len(data_mat)==0 else data_mat

        total_atoms = self.total_atoms
        not_atom = 0
        columns = len(data_mat[9].split())
    
        for line in data_mat:
            try:
                _ = np.array(line.split(), dtype=float)
                if len(_) == columns:
                    new_mat.append(line)
                else:
                    not_atom += 1
            except:
                not_atom += 1
        
        total_steps = len(new_mat)/total_atoms
        if self.PRINT_DETAIL:
            print('Lines of file:\n', len(data_mat))
            print('Lines without coordinates:\n', not_atom)
            print('Total lines with coordinates:\n', len(new_mat))
            print('Total time steps of trajectory:\n', total_steps)
            print('Atoms in every step:\n', total_atoms)
        else:
            print('Reading lammps trajectory file ...\n')
        self.data_only_coords = new_mat
        self.total_steps = total_steps
        
        if int(total_steps)!=total_steps:
            print('Total atoms:', total_atoms)
            print('Total steps:', total_steps)
            return new_mat
        else:
            return new_mat

    #use numpy to get the float data matrix of all information (float_data)
    def get_coords_matrix(self, data_mat=[]):
        _data_mat = list(map(lambda x:x.split(), data_mat))
        try:
            _data_mat = np.array(_data_mat, dtype=float)
            _data_mat = _data_mat[np.argsort(_data_mat[:,0])]
        except:
            raise RuntimeError('Invalid data in data matrix')
        return _data_mat

    # split data matrix by everystep and store them in a list
    def split_by_steps(self, data_mat=[], step_length=1, relax_step=0):

        _new_data = []
        _data = []
        relax_lines = int((relax_step)*self.total_atoms)
        _data = data_mat if data_mat!=[] else self.data_only_coords
        _data = np.array(_data)
        if len(_data)==0:
            _data = self.f_data_only_coords(self.line_list)[relax_lines:]
        else:
            try:
                data_mat.shape
            except:
                try:
                    _data = list(map(lambda x:x.split(), _data[relax_lines:]))
                    _data = np.array(_data, dtype=float)
                except:
                    raise ValueError('Invalid data matrix!')
        step_interv = self.total_atoms*step_length
        for i in range(0, _data.shape[0], step_interv):
            _mat = _data[i:i+self.total_atoms, :]
            if _mat.shape[0]!=self.total_atoms:
                raise RuntimeError('Data incomplete!')
            _mat = _mat[np.argsort(_mat[:,0])]
            _new_data.append(_mat)
            
        self.step_length = step_length
        self.relax_step = relax_step
        self.time_line_list = _new_data
        return _new_data
    
    def get_2d_data_matrix(self, data_mat=[], axis='x'):
        
        column = self.axis.index(axis) + 2
        data_mat[:, column] = 0
        return data_mat
    
    def get_timeline_list(self):
        return self.time_line_list
    
    def get_last_step_data(self):
        self.one_matrix = np.array(list(map(lambda x:x.split(), self.data_matrix)))
        return self.one_matrix[-self.total_atoms:, :]
