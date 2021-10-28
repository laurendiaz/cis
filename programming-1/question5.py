# Question 5: Pivot calibration

# Get data from the input files
import numpy as np
from cartesian import *
from pivotCalibration import *

filename = input()
emPivotData, emPivotSize = readInput_EmPivot(filename + '-empivot.txt')
N_G = emPivotSize(1)
N_frames = emPivotData(2)
# Position of markers with respect to the sensor
G = np.zeros((N_G, 3, N_frames))

ind = 1
for i in np.arange(1, N_frames + 1).reshape(-1):
    G[:, :, i] = emPivotData(np.arange(ind, ind + N_G - 1 + 1))
    ind = ind + N_G

# Part A: Define local probe coordinate system and compute g_i using the first frame
# Compute Midpoint
G1 = G[:, :, 1]

G_0 = np.mean(G1, 1)

# Translate the observations relative to the midpoint
g = np.zeros((N_G, 3))
for i in np.arange(1, N_G + 1).reshape(-1):
    g[i, :] = G1[i, :] - G_0

# Part B&C: Implement Pivot Calibration
t_G, p_dimple = pivotCalibration(g, G)
# For export to output.txt file
q5_exp = p_dimple
