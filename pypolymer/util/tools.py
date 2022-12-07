import numpy as np

# convert a point to float type numpy array
def _cnf(point):
    _p = []
    try:
        _p = np.array(point, dytpe=float)
    except:
        raise ValueError('Invalid point!')
    return _p

# Caculate distance between two points
def cal_distance(a, b):
    result = 0
    try:
        result = np.linalg.norm(a-b)
    except:
        try:
            a,b = _cnf(a), _cnf(b)
            result = np.linalg.norm(a-b)
        except:
            raise ValueError('Invalid two points!')
    return result

# Judge
def _is_in_range(point, center, radius):
    STATUS = True
    try:
        p, c, r = _cnf(point), _cnf(center), radius
        if cal_distance(p, c)>radius:
            STATUS = False
    except:
        raise RuntimeError('Invalid parameters!')
    return STATUS
        