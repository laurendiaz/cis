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

def readOutput(filename):
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



# Product of frame and vector
def frameVecProd(frame, vect):
    return np.matmul(frame.get_rot(), vect) + frame.get_vec


# Inverse of frame
def frameInv(frame):
    '''
        F^-1 = [R^-1, R^-1 dot p]
    '''
    rot = np.inv(frame.get_rot) # is this correct way to ge
    vec = np.matmul(rot, frame.get_vec)
    inv = Frame(rot, vec)
    return inv


# Product of two frames
def frameProd(frame1, frame2):
    '''
        F1 dot F2 =[R1 dot R2, R1p2 + p1]
    '''
    rot = np.matmul(frame1.get_rot(), frame2.get_rot())
    inv = np.matmul(frame1.get_rot(), frame2.get_vec()) + frame1.get_vec()
    prod = Frame(rot, inv)
    return prod


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
