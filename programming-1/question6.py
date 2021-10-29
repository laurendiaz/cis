# Question 6
# Optical tracking pivot calibration

# Get data from the input files
import numpy as np
from cartesian import *
from icp import *
from pivotCalibration import *
import question4
import question5

filename = input()
calBodyData, calBodySize = readInput_Body(filename + '-calbody.txt')
optPivotData, optPivotSize = readInput_OptPivot(filename + '-optpivot.txt')
N_D = optPivotSize[0]
N_H = optPivotSize[1]
N_frames = optPivotSize[2]
d = calBodyData[np.arange(1, N_D + 1), :]

# Position of markers with respect to the sensor
D = np.zeros((N_D, 3, N_frames))
H = np.zeros((N_H, 3, N_frames))

ind = 0
for i in np.arange(0, N_frames):
    D[:, :, i] = optPivotData[np.arange(ind, ind + N_D), :]
    ind = ind + N_D
    H[:, :, i] = optPivotData[np.arange(ind, ind + N_H), :]
    ind = ind + N_H

# Compute F_D = [R_D, p_D]
# R_D is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
# p_D is a 3XN_frames 2D matrix, each column corresponds to the translation of a frame
R_D = np.zeros((3, 3, N_frames))
p_D = np.zeros((3, N_frames))
for i in np.arange(0, N_frames):
    R_i, p_i = ICP(d, D[:, :, i], F0, eta0)
    R_D[:, :, i] = R_i
    p_D[:, i] = p_i

# Transform the optical tracker beacon positions into EM tracker coordinates
P = np.zeros((N_H, 3, N_frames))

for i in np.arange(0, N_frames):
    for j in np.arange(0, N_H):
        P[j, :, i] = np.transpose((np.linalg.solve(R_D[:, :, i], (np.transpose(H[j, :, i]))) - np.linalg.solve(R_D[:, :, i], p_D[:, i])))

# Define local probe coordinate system and compute p_i using the first frame
P1 = P[:, :, 1]
P_0 = np.mean(P1, 1)

# Translate the observations relative to the midpoint
p = np.zeros((N_H, 3))
for i in np.arange(0, N_H):
    p[i, :] = P1[i, :] - P_0

# Implement pivot calibration
t_P, p_dimple = pivotCalibration(p, P)
q6_exp = p_dimple

# Creating output file to export to output folder
filenameout = (filename + '-output-1.txt')
filepathout = ('Output/' + filenameout)
fileID = open(filepathout, 'w')
fileID.write("{}, {}, {}\n".format(N_C, N_framescal, filenameout))
fileID.write("{}, {}\n".format(q5_exp, q6_exp))
fileID.write("{}\n".format(q4_exp))
fileID.close()
