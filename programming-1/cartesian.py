# Create a library of essential functions for 3D point manipulation and Cartesian transformations
import numpy as np
import os

# Read input file
def readInput_Body(filename):
    os.chdir("inputs")
    file = open(filename)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(filename, "r").readline()  # Read the first line from file
    size = f.split(',')[:3]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    os.chdir("..")
    return data, size

def readInput_Readings(filename):
    os.chdir("inputs")
    file = open(filename)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(filename, "r").readline()  # Read the first line from file
    size = f.split(',')[:4]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    os.chdir("..")
    return data, size

def readInput_EmPivot(filename):
    os.chdir("inputs")
    file = open(filename)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(filename, "r").readline()  # Read the first line from file
    size = f.split(',')[:2]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    os.chdir("..")
    return data, size

def readInput_OptPivot(filename):
    os.chdir("inputs")
    file = open(filename)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(filename, "r").readline()  # Read the first line from file
    size = f.split(',')[:3]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    os.chdir("..")
    return data, size


# Define identity matrix
'''
def identity(mat):
    i = np.identity(mat.size)
    return i


# Define rotation matrix - not sure what this is for?
def rotation():
    '''

# Rotation angle and rotation from rotation matrix - not sure how to do this?
def angleAndRotation(rotMatrix):
    return 0

'''
# Product of rotation matrix and vector - use np.matmul
def rotProdMatVec(rotMatrix, vect):
    return 0
'''

# Product of frame and vector
def frameVecProd(frame, vect):
    return np.matmul(frame.get_rot(), vect) + frame.get_vec

'''
# Cross product of two vectors - np.matmul
def crossProd(vec1, vec2):
    return 0

# Dot product of two vectors - np.vdot
def dotProd(vec1, vec2):
    return vdot(vec1, vec2)

# Inverse of rotation matrix - np.inv
def rotInv(rotMatrix):
    return 0
'''
# Inverse of frame
def frameInv(frame):
    '''
        F^-1 = [R^-1, R^-1 dot p]
    '''
    rot = np.inv(frame.get_rot) # is this correct way to ge
    vec = np.matmul(rot, frame.get_vec)
    inv = Frame(rot, vec)
    return inv

'''
# Product of two rotation matrices - np.matmul
def rotProdMatMat(rot1, rot2):
    return 0
'''

# Product of two frames
def frameProd(frame1, frame2):
    '''
        F1 dot F2 =[R1 dot R2, R1p2 + p1]
    '''
    rot = np.matmul(frame1.get_rot(), frame2.get_rot())
    inv = np.matmul(frame1.get_rot(), frame2.get_vec()) + frame1.get_vec()
    prod = Frame(rot, inv)
    return prod

'''
# Least squares solver - scipy
def leastSq():
    return 0
'''

class Frame:
    def __init__(self, rot, vec):
        self.rot = rot
        self.vec = vec
    
    def get_rot(self):
        return self.rot

    def get_vec(self):
        return self.vec
    
    def set_rot(self, x):
        self.rot = x

    def set_vec(self, x):
        self.vec = x