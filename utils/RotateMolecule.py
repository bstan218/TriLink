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

def rotate_molecule(ligand):
    
    def helper(sideindex):
        #define inputs for first rotation
        p11 = ligand.centered_atoms[sideindex]
        p1mag = np.linalg.norm(p11)
        p12 = np.array([p1mag,0,0])
        v1 = FindRotationVector(p11,p12)
        t1 = FindTheta(p11,p12)
        #first rotation
        first_rotated_atoms = []
        for atom in ligand.centered_atoms:
            first_rotated_atoms.append(RotatePoint(atom,v1,t1))
        
        #define inputs for second rotation
        centerindex = ligand.centerconnection[3]-1
        centerconnection = first_rotated_atoms[centerindex]
        p21 = np.array([0,centerconnection[1], centerconnection[2]])
        p2mag = np.linalg.norm(p21)
        p22 = np.array([0,0,p2mag])
        v2 = FindRotationVector(p21,p22)
        t2 = FindTheta(p21,p22)
        #second rotation
        final_rotated_atoms = []
        for atom in first_rotated_atoms:
            final_rotated_atoms.append(RotatePoint(atom, v2, t2))
        
        return final_rotated_atoms

    sideindex = ligand.sideconnections[0][3]-1
    altsideindex = ligand.sideconnections[1][3]-1
    rotated_atoms = helper(sideindex)
    alt_rotated_atoms = helper(altsideindex)

    return rotated_atoms, alt_rotated_atoms

def rock_rotate(template, theta):
    coords = template.optimized_coordinates.coords
    vec = np.array([1,0,0])

    return_coords = []
    for atom in coords:
        return_coords.append(RotatePoint(atom, vec, theta))
    
    return return_coords