# Computes position of the tool tip and the calibration post for pivot calibration
# Input is the position of the points with respect to the probe and the tracker measurements in N different frames
# Outputs two 3x1 vectors ptip and ppivot

import numpy as np
import scipy
from icp import *
import cartesian


def pivotCalibration(j, J):
    N_frames = J.shape[3 - 1]
    ''''
    Calculate transformation for N frames
    R is a 3x3xN_frames 3D matrix, each page corresponding to the rotation matrix of the frame
    p is a 3XN_frames 2D matrix, each column corresponding to the translation of the frame
    '''
    theta = 45
    r = np.array([[1, 0, 0],
                [0, math.cos(theta), -math.sin(theta)],
                [0, math.sin(theta), math.cos(theta)]])
    F0 = cartesian.Frame(r, [1, 1, 1])
    eta0 = 1000000000000000
    R = np.zeros((3, 3, N_frames))
    p = np.zeros((3, N_frames))
    for i in np.arange(1, N_frames + 1).reshape(-1):
        F = ICP(j, J[:, :, i], F0, eta0)  # apply ICP registration
        R_i = F.get_rot()
        p_i = F.get_vec()
        R[:, :, i] = R_i
        p[:, i] = p_i

    '''
    The least square function needs to be solved:
    [Ri|-I]*[ptip;ppivot] = [-pi]
    '''
    # Build Left matrix [Ri|-I] and Right matrix [-pi] for all frames
    Left = []
    Right = []
    for i in np.arange(1, N_frames + 1).reshape(-1):
        Left = np.array([[Left], [R[:, :, i], - np.eye(3)]])
        Right = np.array([[Right], [- p[:, i]]])

    # Left.shape
    # Right.shape

    # Solving least squares
    x, resnorm, residual, exitflag, output, lambda_ = scipy.linalg.lstsq(Left, Right)
    print(resnorm)
    ptip = x(np.arange(1, 3 + 1))
    ppivot = x(np.arange(4, 6 + 1))

    return ptip, ppivot, R, p
