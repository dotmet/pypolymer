import numpy as np
import matplotlib.pyplot as plt
import copy

def cal_line_eq(p1,p2):

    x1,y1 = p1[0], p1[1]
    x2,y2 = p2[0], p2[1]
    dx, dy = x2-x1, y2-y1
    dx = 0.01 if dx==0 else dx
    dy = 0.01 if dy==0 else dy
    A = 1/dx
    B = -1/dy
    C = y1/dy-x1/dx

    def f(x=None,y=None):
        if x is None and y is None:
            return A, B, C
        elif x is None and A!=0:
            return (C-B*y)/A
        elif y is None and B!=0:
            return (C-A*x)/B
        else:
            return A*x+B*y+C

    return f

def cal_cross_point(bond1, bond2):

    p11, p12 = bond1[0], bond1[1]
    p21, p22 = bond2[0], bond2[1]

    A1,B1,C1 = cal_line_eq(p11, p12)()
    A2,B2,C2 = cal_line_eq(p21, p22)()

    rd = A1*B2-A2*B1
    rx = B1*C2-B2*C1
    ry = A2*C1-A1*C2

    if rd==0:
        if A1*C2-A2*C1==0:
            if 'chong he':
                return 'equal'
            else:
                return np.array([])
        else:
            return np.array([])

    else:
        cp = np.array([rx/rd, ry/rd])
        res1 = np.dot(p11-cp, p12-cp)
        res2 = np.dot(p21-cp, p22-cp)
        if res1>0 or res2>0:
            return np.array([])
        else:
            return cp
    
def _cal_area(ps):
    ''' 坐标必须是按顺时针或逆时针排过序的！'''
    s = 0
#     print(ps)
    try:
        ps = ps.tolist()
    except:
        pass
    if len(ps)<=2:
        Warning('Points less than 3 choosed! please check!')
    for p1,p2 in zip(ps, ps[1:]+ps[0:1]):
#         print(p1, p2)
        s += (p1[0]*p2[1]-p1[1]*p2[0])/2

    return np.abs(s)

def split_area(coords, plot=False):

    coords = np.array(coords)
    nps = len(coords)
    if nps<=3:
        return False

    coords = coords.tolist()
    coords = coords+coords[0:1]

    res = []
    p2coords = []

    for i in range(nps-1):
        res.append(coords[i])
        for j in range(i+1, nps):
            bond1 = [coords[i], coords[i+1]]
            bond2 = [coords[j], coords[j+1]]
            cross = cal_cross_point(bond1, bond2)
            if len(cross)==5:
                print('发生重合了')
                coords[i+1][1]+=0.01
                bond1 = [coords[i], coords[i+1]]
                bond2 = [coords[j], coords[j+1]]
                cross = cal_cross_point(bond1, bond2)
            if len(cross)==2:
                if j-i!=1 and j-i!=nps-1:
                    if j+1 >= nps:
                        res = res + [cross]
                    else:
                        res = res + [cross] + coords[j+1:-1]
                    p2coords = [cross]+coords[i+1:j+1]
                    res = np.array(res)
                    p2coords = np.array(p2coords)
                    if plot:
                        plt.plot(res[:,0], res[:,1], color='red', marker='^')
                        plt.plot([res[0,0], res[-1,0]], [res[0,1], res[-1,1]], color='red', marker='^')
                        plt.plot(p2coords[:,0], p2coords[:,1], color='blue', marker='d')
                        plt.plot([p2coords[0,0], p2coords[-1,0]], [p2coords[0,1], p2coords[-1,1]], color='blue', marker='^')
                        plt.show()
                    return res, p2coords

    return False
#                 return res, p2coords
     
def split_area_2d(coords, status=[1], plot=False, check=False):
    
    if len(coords[0])==2:
        coords = [coords]

    _coords = copy.deepcopy(coords)

    for i in range(len(coords)):
        coord = coords[i]

        if status[i]==1 and len(coord)>3:

            if plot:
                plt.plot(coord[:,0], coord[:,1], marker='o')
                plt.show()
            res = split_area(coord, plot=plot)
            
            if check:
                print("------------")
                print('Raw:')
                print(coord)
                print('Points:', len(coord))
                print('Result:')
                print(res)
                print("------------")

            _res = []
            stats = []

            if res is not False:
                c1, c2 = res[0], res[1]
                if len(c1)>=3:
                    if len(c1)>3:
                        stat=1
                    else:
                        stat=0
                    _res.append(c1)
                    stats.append(stat)
                if len(c2)>=3:
                    if len(c2)>3:
                        stat=1
                    else:
                        stat=0
                    _res.append(c2)
                    stats.append(stat)

            else:
                status[i]=0

            if len(_res)!=0 and i!=len(_coords):
                _coords = _coords[:i] + _res + _coords[i+1:]
                status = status[:i]+stats+status[i+1:]
            elif len(_res)!=0 and i==len(_coords):
                _coords = _coords[:i] + _res
                status = status[:i]+stats
                
    return _coords, status

def analysis_2d_shape(coords, plot=False):

    results = {}
    
    coords = [coords]
    status = [1]

    for i in range(10):
        coords, status = split_area_2d(coords, status, plot=plot)
        for i in range(len(coords)):
            if len(coords[i])<=3:
                status[i]=0
        if np.sum(status)==0:
            break

    s = 0 
    for area in coords:
        s += np.abs(_cal_area(area.tolist()))

    results.update({'Area':s})
    results.update({'Coords sub shape':coords})

    return results