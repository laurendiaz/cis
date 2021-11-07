# Question 5
# Pivot calibration

import numpy as np
import cartesian
import pivotCalibration

# Get data from the input files
filename = input()
emPivotData, emPivotSize = cartesian.readInput_EmPivot(filename + '-empivot.txt')
N_G = emPivotSize[0]
N_frames = emPivotSize[1]

F0 = cartesian.Frame(45, [1, 1, 1])
eta0 = 1000

# Position of markers with respect to the sensor
G = np.zeros((N_G, 3, N_frames))

ind = 0
for i in np.arange(0, N_frames):
    G[:, :, i] = emPivotData[np.arange(ind, ind + N_G)]
    ind = ind + N_G

# Part A: Define local probe coordinate system and compute g_j using the first frame
# Compute Midpoint
Gj = G[:, :, 1]
G_0 = np.mean(Gj, 1)

# Translate the observations relative to the midpoint
g = np.zeros((N_G, 3))
for i in np.arange(0, N_G):
    g[i, :] = Gj[i, :] - G_0[i]

# Part B and C: Implement Pivot Calibration
t_G, P_dimple = pivotCalibration.pivotCalibration(g, G)

# For export to output.txt file
q5_exp = P_dimple
