# Question 4
# Solve for distorted calibration

import numpy as np
from cartesian import *
from icp import *

# Get user input for file name:
print('Enter the name of your input data file (e.g. pa1-debug-a): ')
filename = input()

#Get data from the input files
calBodyData,calBodySize = readInput_Body(filename + '-calbody.txt')
calReadData,calReadSize = readInput_Readings(filename + '-calreadings.txt')
N_D = calBodySize(1)
N_A = calBodySize(2)
N_C = calBodySize(3)
N_framescal = calReadSize(4)

# Position of markers with respect to the calibration body
d = calBodyData[np.arange(N_D+1,1),:]
a = calBodyData[np.arange(N_D + 1,N_A + N_D+1),:]
c = calBodyData[np.arange(N_D + N_A + 1,len()-1),:]

# Position of markers with respect to the trackers
D = np.zeros((N_D,3,N_framescal))
A = np.zeros((N_A,3,N_framescal))
C = np.zeros((N_C,3,N_framescal))

ind = 1
for i in np.arange(1,N_framescal+1).reshape(-1):
    D[:,:,i] = calReadData[np.arange(ind,ind + N_D - 1+1),:]
    ind = ind + N_D
    A[:,:,i] = calReadData[np.arange(ind,ind + N_A - 1+1),:]
    ind = ind + N_A
    C[:,:,i] = calReadData[np.arange(ind,ind + N_C - 1+1),:]
    ind = ind + N_C

# Part A: Compute F_D = [R_D, p_D]
# R_D is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
# p_D is a 3XN_frames 2D matrix, each column corresponds to the translation
# of a frame
R_D = np.zeros((3,3,N_framescal))
p_D = np.zeros((3,N_framescal))
for i in np.arange(1,N_framescal+1).reshape(-1):
    R_i,p_i = ICP(d,D[:,:,i])
    R_D[:,:,i] = R_i
    p_D[:,i] = p_i

# Part B: Compute F_A = [R_A, p_A]
# R_A is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
# p_A is a 3XN_frames 2D matrix, each column corresponds to the translation
# of a frame
R_A = np.zeros((3,3,N_framescal))
p_A = np.zeros((3,N_framescal))
for i in np.arange(1,N_framescal+1).reshape(-1):
    R_i,p_i = ICP(a,A[:,:,i])
    R_A[:,:,i] = R_i
    p_A[:,i] = p_i

## Part C: Compute C_i expected = inv(R_D) * (R_A*ci + p_A - p_D)
C_exp = np.zeros((N_C,3,N_framescal))
for i in np.arange(1,N_framescal+1).reshape(-1):
    for j in np.arange(1,N_C+1).reshape(-1):
        C_exp[j,:,i] = np.transpose((np.linalg.inv(R_D[:,:,i]) * (R_A[:,:,i] * np.transpose(c(j,:)) + p_A[:,i] - p_D[:,i])))

## Part D: Prepare for output file
# Reshaping C_exp
C_fin = []
for i in np.arange(1,N_framescal+1).reshape(-1):
    C_fin = np.array([[C_fin],[C_exp[:,:,i]]])

# For export to output.txt file
q4_exp = np.transpose(C_fin)