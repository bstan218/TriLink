import numpy as np
import math as m

def FindRotationVector(p1,p2):
    #finds (unit vector) rotation vector based on two points
    vec = np.cross(p1,p2)
    mag = np.linalg.norm(vec)
    return vec/mag

def FindTheta(p1,p2):
    #finds theta based on two points
    m1 = np.linalg.norm(p1)
    m2 = np.linalg.norm(p2)
    theta = np.arccos(np.dot(p1,p2)/(m1*m2))
    return theta
                  
def RotatePoint(p, rotationvector, theta):
    #rotate a point around axis v by theta
    q = np.array([0.0,0.0,0.0])

    c = np.cos(theta)
    t = (1 - np.cos(theta))
    s = np.sin(theta)
    X = rotationvector[0]
    Y = rotationvector[1]
    Z = rotationvector[2]

    d11 = t*X**2 + c
    d12 = t*X*Y - s*Z
    d13 = t*X*Z + s*Y
    d21 = t*X*Y + s*Z
    d22 = t*Y**2 + c
    d23 = t*Y*Z - s*X
    d31 = t*X*Z - s*Y
    d32 = t*Y*Z + s*X
    d33 = t*Z**2 + c

    q[0] = round((d11*p[0] + d12*p[1] + d13*p[2]),4)
    q[1] = round((d21*p[0] + d22*p[1] + d23*p[2]),4)
    q[2] = round((d31*p[0] + d32*p[1] + d33*p[2]),4)

    return q



if __name__ == "__main__":
    vec = np.array([1,0,0])
    theta = 5
    p1 = np.array([2,0,0])
    print(RotatePoint(p1,vec, theta))
    
    