import numpy as np

class RawData(object):
    
    def __init__(self, file_name=''):
        self.file_name = file_name
        self.raw_data = ''
        self.line_list = []
        self.data_only_coords = []
        
    def read_file(self, file_name=''):
        if file_name == '':
            file_name = self.file_name
        try:
            with open(file_name, 'r') as f:
                data = f.readlines()
                self.line_list = data
                f.close()
        except:
            raise ValueError('Invalid input file')
        return self.line_list
    
    def get_raw_data(self, file_name=''):
        try:
            with open(file_name, 'r') as f:
                self.raw_data = f.read()
        except:
            raise ValueError('Invalid input file')
        return self.raw_data
    
    def get_line_list(self):
        return self.line_list
        
    def get_atoms(self):
        atoms = 0
        with open(self.file_name, 'r') as f:
            i = 0
            for i in range(9):
                l = f.readline()
                if i==3:
                    atoms = int(l.split()[0])
                i+=1
        if atoms==0:
            raise ValueError('Wrong atoms find!')
        return atoms

            
                