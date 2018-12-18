import numpy as np
import matplotlib.pyplot as plt
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
        return Triangle(value @ self.p1, value @ self.p2, value @ self.p3)

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
        return
    

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


