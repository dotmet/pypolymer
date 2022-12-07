import numpy as np
import matplotlib.pyplot as plt

class Triangle():
    def __init__(self, coords):
        
        '''
        The coords[1] are choosed to be the apex(顶点) of the triangle.
        '''
        
        coords = np.array(coords)
        if coords.shape[0]!=3:
            raise ValueError('Wrong triangle defined!')
            
        self.coords = coords
            
        self.p1 = coords[0]
        self.p2 = coords[1] # The apex of one angle
        self.p3 = coords[2]
        
        self.nvec = self.cal_nvec()
        self.area = self.cal_area()
        
        self.center_of_hypotenuse = coords[0]/2+coords[2]/2 #斜边中点
        self.apex_angle = self.cal_angle(self.p1-self.p2, self.p3-self.p2) # 顶角
        
    def cal_nvec(self, coords=None):
        if coords is None:
            coords = self.coords
        nvec = np.cross(coords[0]-coords[1], coords[2]-coords[1])
        return nvec
    
    def cal_area(self, coords=None):
        if coords is None:
            coords = self.coords
        nvec = self.cal_nvec(coords)
        area = np.linalg.norm(nvec)/2
        return area
    
    def cal_angle(self, p1, p2):
        module = np.linalg.norm(p1)*np.linalg.norm(p2)
        if module==0:
            print('The result may be wrong or any value!')
            print(p1, p2)
            return None
        return np.arccos(np.dot(p1,p2)/module)
    
    def apex_vs_rcm(self, rcm=[0,0,0]):
        
        rcm = np.array(rcm)
        d1 = np.linalg.norm(self.p1-rcm)
        d2 = np.linalg.norm(self.p2-rcm)
        d3 = np.linalg.norm(self.p3-rcm)
        
        if d2==d1 and d2==d3:
            return 'line'
        elif d2<np.max([d1,d3]):
            return False
        else:
            return True

def split_by_triangles(coords):
    
    nps = len(coords)
    
    if nps==3:
        return coords, [Triangle(coords)]
    
    elif nps==4:
        return coords, [Triangle(coords[:2, :]), Triangle(coords[[2,3,0], :])]
    
    elif nps%2!=0 :
        # print('Add points')
        coords = np.vstack([coords, coords[0]/2+coords[-1]/2])
        nps += 1
    
    new_apex_coords = []
    triangles = []
    remain_apex = []
    
    # Find next apex of angles
    coords = np.vstack([coords, coords[0]]) # Add first point to end to make chain to be closed
    for i in range(0,nps-1,2):
        tricoords = coords[i:i+3]
        trag = Triangle(tricoords)
        new_apex_coords.append(trag.center_of_hypotenuse)
        triangles.append(trag)

        remain_apex.append(coords[i]) # Record the apex points that doesnot counted
    
    # Store all triangles between the non-apex (odd) points of input coords and 
    # the apex points for output
    for i in range(len(remain_apex)):
        if i==len(remain_apex)-1:
            i=0
        tri_coords = [new_apex_coords[i], remain_apex[i+1], new_apex_coords[i+1]]
        trag = Triangle(np.array(tri_coords))
        triangles.append(trag)
    
    return np.array(new_apex_coords), triangles

def split_area_3d(coords, plot=False):
    ''' The coords must be sorted! '''
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        
    triangles = []
    _coords = [coords]
    for i in range(5*len(coords)):
        
        if plot:
            pcoords = np.vstack([coords, coords[-1]])
            ax.plot(pcoords[:,0], pcoords[:,1], pcoords[:,2])

        res = split_by_triangles(coords)
        coords = res[0]
        triangles += res[1]
        if len(res[1])<=2:
            break
        else:
            _coords.append(coords)
    if plot:
        plt.show()
    return _coords, triangles

def analysis_3d_shape(coords, plot=False, detail_info=False):

    results = {}
    nvec = 0
    area = 0
    _coords, trags = split_area_3d(coords, plot=plot)
    for tag in trags:
        nvec += tag.nvec
        area += tag.area
        
    results.update({'Area':area})
    results.update({'Normal vector':nvec})

    if detail_info:
        results.update({'All apex points':_coords})
        results.update({'All triangles':trags})
        results.update({'All shapes coords': _coords})

    return results