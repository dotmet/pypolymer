import numpy as np

def read_position(file_name):
    f = open(file_name, 'r')
    xyz = []
    for i in f:
        s = list(map(lambda x: float(x), i.split()))
        xyz.append(s)
    f.close()
    xyz = np.array(xyz, dtype=np.float32)
    return xyz

def read_Data(file_name):
    f = open(file_name, 'r')
    data = []
    for i in f:
        s = list(map(lambda x: float(x), i.split()))
        data.append(s)
    f.close()
    data = np.array(data, dtype=int)
    return data

def caculate_distance(xyz):
    distance = []
    k=1
    for i in range(len(xyz)):
        for j in range(i+1, len(xyz)):
            d = np.linalg.norm(xyz[i]-xyz[j])
            if abs(d-1)<0.3:
                d = 'bond_id: '+str(k)+' ('+str(i+1)+', '+str(j+1)+'): '+str(d)
                distance.append(d)
                k=k+1
    return distance

def bond_length_file(file_name, distance):
    with open(file_name, 'w') as f:
        for i in distance:
            f.write(str(i)+'\n')
        f.close()