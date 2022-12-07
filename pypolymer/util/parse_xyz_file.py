import numpy as np

def extract_chain(file_name='data.xyz', atom_typ=3.0):
    data = []
    raw_data = ''
    with open(file_name, 'r') as f:
        _d = f.readlines()
        data = np.array([l.split() for l in _d[2:]], dtype=float)
        raw_data = _d
        f.close()
    _mat = data[np.where(data[:,0]==3.0)]
    print('Dimension of matrix of chain\n', _mat.shape)
    return _mat, raw_data

def split_chain(chain, mol_length):
    _c = [chain[i:i+mol_length, :] for i in range(chain.shape[0])]
    return _c

def _mol_out(mol, r_ange):
    dist = [np.linalg.norm(l[1:]) for l in mol]
    if max(abs(np.array(dist)))>r_ange:
        return True
    else:
        return False
    
def main():
    fn = 'last.xyz'
    at = 3.0
    R = 15
    chain, ftext = extract_chain(file_name=fn, atom_typ=at)
    _l = []
    for line in ftext:
        ls = line.split()
        if ls[0]!='3':
            _l.append(line)
            continue
        else:
            if np.linalg.norm(ls[1:])<=R:
                _l.append(line)
            else:
                continue
    with open('parsed_'+fn, 'w') as f:
        f.write(''.join(_l))
        f.close()

if __name__ == '__main__':
    main()
