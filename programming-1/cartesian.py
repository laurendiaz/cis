# Create a library of essential functions for 3D point manipulation and Cartesian transformations
import numpy as np
import os

# Read input file
def readInput_Body(filename):
    #fname = os.getcwd() + "/programming-1/inputs/" + filename
    fname = os.getcwd() + "/inputs/" + filename
    file = open(fname)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(fname, "r").readline()  # Read the first line from file
    size = f.split(',')[:3]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    return data, size

def readInput_Readings(filename):
    #fname = os.getcwd() + "/programming-1/inputs/" + filename
    fname = os.getcwd() + "/inputs/" + filename
    file = open(fname)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(fname, "r").readline()  # Read the first line from file
    size = f.split(',')[:4]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    return data, size

def readInput_EmPivot(filename):
    #fname = os.getcwd() + "/programming-1/inputs/" + filename
    fname = os.getcwd() + "/inputs/" + filename
    file = open(fname)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(fname, "r").readline()  # Read the first line from file
    size = f.split(',')[:2]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    return data, size

def readInput_OptPivot(filename):
    #sfname = os.getcwd() + "/programming-1/inputs/" + filename
    fname = os.getcwd() + "/inputs/" + filename
    file = open(fname)
    # Get data
    data = np.loadtxt(file, delimiter=',', skiprows=1, dtype=str)
    # Get size
    f = open(fname, "r").readline()  # Read the first line from file
    size = f.split(',')[:3]  # extract the first three values
    size = list(map(int, size))  # Convert list to Array of int
    return data, size

# Product of frame and vector
def frameVecProd(frame, vect):
    dotprod = np.dot(frame.get_rot(), vect)
    return np.add(dotprod, frame.get_vec())


# Inverse of frame
def frameInv(frame):
    '''
        F^-1 = [R^-1, R^-1 dot p]
    '''
    rot = np.inv(frame.get_rot()) # is this correct way to ge
    vec = np.dot(rot, frame.get_vec())
    inv = Frame(rot, vec)
    return inv


# Product of two frames
def frameProd(frame1, frame2):
    '''
        F1 dot F2 =[R1 dot R2, R1p2 + p1]
    '''
    rot = np.dot(frame1.get_rot(), frame2.get_rot())
    dotprod = np.dot(frame1.get_rot(), frame2.get_vec())
    inv = np.add(dotprod, frame1.get_vec())
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
