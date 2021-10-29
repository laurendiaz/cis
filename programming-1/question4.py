# Question 4
# Solve for distorted calibration

import numpy as np
from cartesian import *
from icp import *
import math

# Get user input for file name:
print('Enter the name of your input data file (e.g. pa1-debug-a): ')
filename = input()
theta = 45
R = np.array([[1, 0, 0],
            [0, math.cos(45), -math.sin(45)],
            [0, math.sin(45), math.cos(45)]])
F0 = Frame(R, [1, 1, 1])
eta0 = 1000000000000000

# Get data from the input files
calBodyData, calBodySize = readInput_Body(filename + '-calbody.txt')
calReadData, calReadSize = readInput_Readings(filename + '-calreadings.txt')
N_D = calBodySize[0]
N_A = calBodySize[1]
N_C = calBodySize[2]
N_framescal = calReadSize[3]

# Position of markers with respect to the calibration body
d = calBodyData[np.arange(1, N_D + 1), :]
d = d.astype(float)
a = calBodyData[np.arange(N_D + 1, N_A + N_D + 1), :]
a = a.astype(float)
c = calBodyData[np.arange(N_D + N_A + 1, len(calBodyData) - 1), :]
c = c.astype(float)

# Position of markers with respect to the trackers
D = np.zeros((N_D, 3, N_framescal))
A = np.zeros((N_A, 3, N_framescal))
C = np.zeros((N_C, 3, N_framescal))

ind = 0
for i in np.arange(1, N_framescal).reshape(-1):
    D[:, :, i] = calReadData[np.arange(ind, ind + N_D), :]
    ind = ind + N_D
    A[:, :, i] = calReadData[np.arange(ind, ind + N_A), :]
    ind = ind + N_A
    C[:, :, i] = calReadData[np.arange(ind, ind + N_C), :]
    ind = ind + N_C

# Part A: Calculate F_D = [R_D, p_D]
# R_D is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
# p_D is a 3XN_frames 2D matrix, each column corresponds to the translation of a frame
R_D = np.zeros((3, 3, N_framescal))
p_D = np.zeros((3, 3, N_framescal))
for i in np.arange(0, N_framescal):
    F = ICP(d, D[:, :, i], F0, eta0)
    R_i = F.get_rot()
    p_i = F.get_vec()
    R_D[:, :, i] = R_i
    print(p_D.shape)
    print(p_i.shape)
    p_D[:, :, i] = p_i

# Part B: Calculate F_A = [R_A, p_A]
# R_A is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
# p_A is a 3XN_frames 2D matrix, each column corresponds to the translation of a frame
R_A = np.zeros((3, 3, N_framescal))
p_A = np.zeros((3, 3, N_framescal))
for i in np.arange(0, N_framescal):
    F = ICP(a, A[:, :, i], F0, eta0)
    R_i = F.get_rot()
    p_i = F.get_vec()
    R_A[:, :, i] = R_i
    p_A[:, :, i] = p_i

# Part C: Compute C_i expected = inv(R_D) * (R_A*ci + p_A - p_D)
C_exp = np.zeros((N_C, 3, N_framescal))
for i in np.arange(0, N_framescal):
    for j in np.arange(0, N_C):
        C_exp[j, :, i] = np.transpose(np.linalg.inv(R_D[:, :, i]) * (R_A[:, :, i] * (c[j]) + p_A[:, :, i] - p_D[:, :, i]))

# Part D: Output C_i expected
# Reshaping C_exp
C_final = []
for i in np.arange(0, N_framescal):
    C_final = np.array([[C_final], [C_exp[:, :, i]]])

# Export to output.txt file
q4_exp = np.transpose(C_final)


