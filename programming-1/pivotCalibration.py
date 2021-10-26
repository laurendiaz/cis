# Computes position of the tool tip and the calibration post for pivot calibration
# Input is the position of the points with respect to the probe and the tracker measurements in N different frames
# Outputs two 3x1 vectors ptip and ppivot

import numpy as np

def pivotCalibration(j=None, J=None):
    N_frames = J.shape[3 - 1]
    ## Find the transformation for N frames
    # R is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix in a frame
    # p is a 3XN_frames 2D matrix, each column corresponds to the
    # translation in a frame
    R = np.zeros((3, 3, N_frames))
    p = np.zeros((3, N_frames))
    for i in np.arange(1, N_frames + 1).reshape(-1):
        R_i, p_i = registration(j, J(:,:, i))
        R[:, :, i] = R_i
        p[:, i] = p_i

    ## Costructing functions for least squares
    # The least square function needs to be solved:
    # [Ri|-I]*[p_tip;p_pivot] = [-pi]

    # Build LHS matrix [Ri|-I] and RHS [-pi] for all frames
    LHS = []
    RHS = []
    for i in np.arange(1, N_frames + 1).reshape(-1):
        LHS = np.array([[LHS], [R(:,:, i), - np.eye(3)]])
        RHS = np.array([[RHS], [- p(:, i)]])

        LHS.shape
        RHS.shape
        ## Solving least squares
        x, resnorm, residual, exitflag, output, lambda_ = leastSquareSolver(LHS, RHS)
        # MAKE LEAST SQUARES SOLVER
        print(resnorm)
        ptip = x(np.arange(1, 3 + 1))
        ppivot = x(np.arange(4, 6 + 1))
        return ptip, ppivot, R, p

        return ptip, ppivot, R, p