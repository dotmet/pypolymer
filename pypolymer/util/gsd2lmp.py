import gsd.hoomd
import argparse
import sys
import numpy as np
from pypolymer.util import parse_boundary
from tqdm import tqdm

class trjconvert:
    
    def __init__(self, infile, outfile, stride=None, parse_bound=False):
        self.infile = infile
        self.outfile = outfile
        self.stride = stride
        self.parse_boundary = parse_bound
        
        self.inext = self.infile.split(".")[-1]
        self.outext = self.outfile.split(".")[-1]
        
        self.lammps_head = self.gen_lammps_head()
        self.lmp_keys = [k for k in self.lammps_head.keys()]
        
        if self.inext == "lammpstrj":
            pass
        elif self.inext == "gsd":
            pass
        
        if self.outext == "lammpstrj":
            pass
        elif self.outext == "gsd":
            pass
    
    def gen_lammps_head(self):
        head = {"ITEM: TIMESTEP": 0,
                "ITEM: NUMBER OF ATOMS":0,
                "ITEM: BOX BOUNDS ":'',
                "boundary rule": ["pp", "pp", "pp"],
                "boundary box": '',
                "ITEM: ATOMS id type x y z":'',}
        return head
    
    def gen_lammps_text(self, head):
        text = ''
        head["ITEM: BOX BOUNDS "] = " ".join(head["boundary rule"])+ "\n" + \
        self.convert_to_str(head["boundary box"])
        for key,val in head.items():
            if "boundary" in key:
                continue
            if "BOX" in key:
                text += key +self.convert_to_str(val) + "\n"
            else:
                text += key + "\n" +self.convert_to_str(val) + "\n"
        return text
    
    def convert_to_str(self, arr):
        arr = np.array(arr).astype(str)
        if arr.ndim == 0:
            return str(arr)
        if arr.ndim == 1:
            return " ".join(arr)
        else:
            return "\n".join([" ".join(l) for l in arr])
        
    def read_lammpstrj(self):
        raw_data = None
        with open(self.infile, 'r') as f:
            raw_data = f.readlines()
            f.close()
            
        trj_data = []
        
        sid = 0
        atoms = int(raw_data[3])
        
        for i in range(100000):
            keys = self.lmp_keys
            head = self.gen_lammps_head()
            head[keys[0]] = int(raw_data[1])
            head[keys[1]] = atoms
            head[keys[3]] = raw_data[4].split("BOUNDS")[-1].split()
            head[keys[4]] = [[float(x) for x in l.split()] for l in raw_data[5:8]]
            head[keys[5]] = np.array(
                [l.split() for l in raw_data[sid+8:sid+8+atoms]]).astype(float)
            atoms = int(raw_data[sid+8+atoms+3])
            trj_data.append(head)
        return trj_data
    
    def read_gsd(self):
        return gsd.hoomd.open(self.infile, 'rb')

    def write_lammpstrj(self):
        snaps = self.read_gsd()
        trj_data = []
        with open(self.outfile, 'w') as f:
            steps = len(snaps)
            for stp in tqdm(range(steps)):
                snap = snaps[stp]
                keys = self.lmp_keys
                head = self.gen_lammps_head()
                head[keys[0]] = snap.configuration.step
                head[keys[1]] = snap.particles.N
                # head[keys[3]] = raw_data[4].split("BOUNDS")[-1].split()
                head[keys[4]] = [[-l/2, l/2] for l in snap.configuration.box[:3]]
                atomid_typ = np.vstack([np.arange(snap.particles.N), 
                                        snap.particles.typeid]).T + 1
                posi = snap.particles.position
                bonds = snap.bonds.group
                box = snap.configuration.box
                if self.parse_boundary:
                    posi = parse_boundary(posi, bonds, box)
                atoms = np.concatenate([atomid_typ, 
                                  posi], axis=1, dtype=object)
                head[keys[5]] = atoms
                # trj_data.append(head)
                    # for head in trj_data[::self.stride]:
                f.write(self.gen_lammps_text(head))
        f.close()
        snaps.close()
    
    def write_gsd(self):
        pass
    
import argparse
def parse_args():
    
    parser = argparse.ArgumentParser(description='Threading')
    
    parser.add_argument('-inp', type=str, default="", help='input file you need to convert. ".gsd"')
    parser.add_argument('-out', type=str, default="", help='output file you need ot convert to. ".lammpstrj"')
    parser.add_argument('-stride', type=int, default=1, help='convert per stride steps.')
    parser.add_argument("--parse_boundary", action="store_true", help="Parse boundary when convert trajectory.")
    
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    
    args = parse_args()
    
    if args.out=="":
        args.out = "".join((args.inp).split(".")[:-1]+[".lammpstrj"])
    
    tjv = trjconvert(args.inp, args.out, args.stride, args.parse_boundary)
    tjv.write_lammpstrj()
    print(f"\n*** Convert [{args.inp}] to [{args.out}] success! ***\n")
    
