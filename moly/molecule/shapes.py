"""
Produce basic shapes for molecules

"""

import numpy as np

def rotation_matrix(vec1, vec2):
    """ Find the rotation matrix that aligns vec1 to vec2
    :param vec1: A 3d "source" vector
    :param vec2: A 3d "destination" vector
    :return mat: A transform matrix (3x3) which when applied to vec1, aligns it with vec2.
    """
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))
    return rotation_matrix


def get_single_cylinder(radius, points=100):
    phi = np.linspace(0, 2*np.pi, points)
    x   = radius * np.cos(phi)
    y   = radius * np.sin(phi)
    
    data = []
    for z in np.linspace(0, 1, 2):
        z = np.ones(points) * z
        data.extend(np.vstack([x,y,z]).T)
    data = np.vstack(data)
    
    return data


def get_bonds_cylinders(geometry, bonds, radius=0.3, points=60):
    
    cilinder_a = []
    cilinder_b = []
        
    for idx1, idx2 in bonds:
            
        vec1 = geometry[idx1]
        vec2 = geometry[idx2]
        length = np.linalg.norm(vec2 - vec1)
        
        #mid = (vec1 + vec2)/2
        #r   = 0.3
        #resolution = 50 #n-1 sides
        
        #Create a circle
        phi = np.linspace(0, 2*np.pi, resolution)
        x   = radius * np.cos(phi)
        y   = radius * np.sin(phi)
        
        for bond in range(2):
        
            #Extend circle in z
            #z = np.linspace(-length/2 + (length/2)*bond,  length/2 - (length/2)*(1-bond)  , 2)           
            data = []
            if bond == 0:
                for z in np.linspace(0, length/2, 2):
                    z = np.ones(resolution) * z
                    data.extend(np.vstack([x,y,z]).T)
                
            if bond == 1:
                for z in np.linspace(length/2, length , 2):
                    z = np.ones(resolution) * z
                    data.extend(np.vstack([x,y,z]).T)
                
            data = np.vstack(data)

            #Rotate cilinder
            R = rotation_matrix(np.array([0,0,1]), vec2 - vec1)
            data = R.dot(data.T).T
            data += vec1
        
            if bond == 0:
                cilinder_a.append(data)
            if bond == 1:
                cilinder_b.append(data)
        
    return [cilinder_a, cilinder_b]


def get_sphere(r=1.0, points=20):

    phi   = np.linspace(0,        2*np.pi, 2*points)
    theta = np.linspace(-np.pi/2, np.pi/2, points)

    #Mesh of sphere
    phi, theta = np.meshgrid(phi[1:], theta)

    xsphere = r * (np.cos(theta) * np.sin(phi)).flatten()
    ysphere = r * (np.cos(theta) * np.cos(phi)).flatten()
    zsphere = r * (np.sin(theta)).flatten()
    
    return [xsphere,ysphere,zsphere]

def get_atoms_spheres(symbols, atomic_numbers, cube=False):

    spheres = []

    for atom in range(len(symbols)):
        radius = atomic_numbers[atom]/30 + 0.6
        sphere = get_sphere(r=radius)
        spheres.append(sphere)

    return spheres




