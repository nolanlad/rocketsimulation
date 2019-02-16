import numpy as np
import matplotlib.pyplot as plt

#drag coefficient of a flat plate
# source : https://www.grc.nasa.gov/www/k-12/InteractProgs/index.htm
# source : https://www.grc.nasa.gov/www/k-12/airplane/kiteincl.html
CD = 1.28

def drag_force(rho, velocity, drag_coefficient , area):
    #source : https://en.wikipedia.org/wiki/Drag_equation
    return 0.5*rho*(velocity**2) * drag_coefficient * area


def fin_drag(F,v,rho):
    """drag on a fin (F) from velocity vector (v)"""
    v_normal = normal_plane(v)
    F_normal = finmatmul(v_normal,F)
    A_normal = F_normal.area()
    C_normal = F_normal.centroid()
    drag = drag_force(rho,mag(v),CD,A_normal)
    return C_normal, drag
    
    

'''
the traingle object
all 2D shapes can be represented as a collection of triangles
triangles are chosen to constitute the fin because it is easy to calculate
their geometric properties
'''
class Triangle:
    def __init__(self,p1,p2,p3):
        self.p1 = np.array(p1)
        self.p2 = np.array(p2)
        self.p3 = np.array(p3)
        
    def area(self):
        ''' returns the area of the triangle
        this works because for any two vectors the magnitude of the
        cross product is the area of the parallelogram spanned by the
        two vectors. Therefore a triangle will be half that area
        '''
        return 0.5*mag(np.cross(self.p2-self.p1,self.p3-self.p1))
    
    def plot(self):
        '''
        plot the triangle to a graph'''
        plt.plot([self.p1[0],self.p2[0]],[self.p1[1],self.p2[1]],'b-')
        plt.plot([self.p1[0],self.p3[0]],[self.p1[1],self.p3[1]],'b-')
        plt.plot([self.p2[0],self.p3[0]],[self.p2[1],self.p3[1]],'b-')
        
    def __matmul__(self,value):
        ''' implements the @ operator so that a matrix multiplication
        can be distributed to each point in a Triangle '''
        return Triangle(value @ self.p1, value @ self.p2, value @ self.p3)

    def __mul__(self,other):
        return Triangle(other * self.p1, other * self.p2, other * self.p3)

class Fin(list):
    def __init__(self):
        super().__init__()

    def area(self):
        A = 0.0
        for f in self:
            A+= f.area()
        return A

    def centroid(self):
        centroid = np.zeros(3)
        for f in self:
            centroid += tri_centroid(f)*f.area()
        return centroid

    def plot(self):
        for f in self:
            f.plot()
            
    

def trimatmul(mat,T):
    '''multiply each point in triangle "T" by the matrix "mat" '''
    return Triangle(mat @T.p1, mat @T.p2, mat @T.p3)

def trimul(n,T):
    '''multiply each point in triangle "T" by the matrix "mat" '''
    return Triangle(mat @T.p1, mat @T.p2, mat @T.p3)

def triadd(v,T):
    ''' add vector v to each point in a triangle
    good for translating shapes'''
    return Triangle(np.add(v,T.p1), np.add(v,T.p2), np.add(v,T.p3))
    


def mag(v):
    '''magnitude of a vector '''
    return np.sqrt(np.dot(v,v))

def unit(v):
    ''' returns unit of vector "v" '''
    return np.array(v)/(np.sqrt(np.dot(v,v)))


def normal_plane(v):
    '''
    returns a matrix representing the space of the plane with "v" as its normal
    '''
    a = v[0]
    b = v[1]
    c = v[2]
    if a != 0:
        return np.array(
            [ [-1*b/a,1,0],
                 [-1*c/a,0,1],
                 [0, 0, 0] ])
    elif b != 0:
        return np.array([ [1,-a/b,0],
                 [0,-c/b,1],
                 [0, 0, 0 ]])
    elif c != 0:
        return np.array([[1,0,-a/c],
                [0,1,-b/c],
                [0, 0, 0 ]])
    else:
        print("you moron, are you actually trying to find the plane normal to a null vector?")
        return 0
    

def rot(axis,thet):
    '''
    returns the rotation matrix of an angle theta around the axis "axis"
    '''
    x = axis[0]
    y = axis[1]
    z = axis[2]
    cost = np.cos(thet)
    sint = np.sin(thet)
    xx = x*x
    yy = y*y
    zz = z*z
    xy = x*y
    xz = x*z
    yz = y*z
    return np.array([
        [cost + xx*(1 - cost),xy*(1-cost) - z*sint, xz*(1-cost)+y*sint],
        [xy*(1 - cost) + z*sint,cost+yy*(1 - cost),yz*(1 - cost) - xx*sint],
        [xz*(1 - cost) - y*sint, x*np.sin(thet)+yz*(1 - cost),cost+zz*(1 - cost)]])


def example():
    T = Triangle([0,0,0],[0,1,0],[0,0,1])
    rot45 = rot([0,1,0],np.pi/4.)
    rot90 = rot([0,1,0],np.pi/2.)
    T45 = trimul(rot45,T)
    T90 = trimul(rot90,T)
    T.plot()
    T45.plot()
    T90.plot()
    plt.show()
    down = normal_plane([0,0,1]) # project looking up the z axis
    T45_flat = trimul(down,T45)
    T90_flat = trimul(down,T90)
    T_flat   = trimul(down,T)
    print("the area of T is %f"%(T_flat.area()))
    print("the area of T45 is %f"%(T45_flat.area()))
    print("the area of T90 is %f"%(T90_flat.area()))


def tri_centroid(T):
    '''
    find the centroid of a triangle, useful for finding drag force
    source:
    https://www.mathopenref.com/coordcentroid.html
    '''
    Ox = (T.p1[0] + T.p2[0] + T.p3[0])/3.
    Oy = (T.p1[1] + T.p2[1] + T.p3[1])/3.
    Oz = (T.p1[2] + T.p2[2] + T.p3[2])/3.
    return np.array([Ox, Oy, Oz])


def make_basic_rocket():
    fin1 = Triangle([0,0,0],[0,0,1],[0,1,0])
    fin1 = triadd([0,1,0],fin1)
    rot90 = rot([0,0,1],np.pi/2)
    fin2 = trimatmul(rot90,fin1)
    fin3 = trimatmul(rot90,fin2)
    fin4 = trimatmul(rot90,fin3)
    F = Fin()
    F.append(fin1)
    F.append(fin2)
    F.append(fin3)
    F.append(fin4)
    F = finadd([0,0,-2],F)
    rot45 = rot([0,1,0],np.pi/4)
    
    return F

def finmatmul(v,F):
    F2 = Fin()
    for f in F:
        F2.append(trimatmul(v,f))
    return F2

def finadd(v,F):
    F2 = Fin()
    for f in F:
        F2.append(triadd(v,f))
    return F2
