import cartesian
import numpy as np
import icp
import math
import pivotCalibration
import scipy

'''
In this assignment, you will fit polynomials to model distortion and then use these polynomials to “dewarp” 
the EM tracker space. You will use the dewarped EM tracker space to repeat your pivot calibration.
Then you will use the resulting calibration to compute registration to a CT coordinate system. Finally,
given some pointer data frames, you will report corresponding CT coordinates.
    1. Input the body calibration data file and process it to determine the values of C_i^expected [k]
        corresponding to each C_i [k] in each “frame” k of data.
    2. Use the method described in class to produce a suitable distortion correction function.
    3. Use this distortion correction function to repeat your “pivot calibration” for the EM probe.
    4. Using the distortion correction and the improved pivot value, compute the locations b_j of the
        fiducials points with respect to the EM tracker base coordinate system.
    5. Compute the registration frame F_reg.
    6. Apply the distortion correction to all the subsequent values of [G_1, ..., G_N_G], compute the pointer
        tip coordinates with respect to the tracker base, and apply F_reg to compute the tip location with
        respect to the CT image.
'''


def ScaleToBox(q, qmin, qmax):
    return (q - qmin) / (qmax - qmin)


def Bernie(a, b):
    return scipy.special.comb(5, b) * (a ** b) * ((1 - a) ** (5 - b))

def distortionCorrection(p, q):
    '''
        Compute distortion correction function for a distorted
        3D Navigational Sensor
        Construct a "tensor form" interpolation polynomial using 5th degree
        Bernstein polynomials F_ijk(u_x, u_y, u_z) = B_5,i(u_x)B_5,j(u_y)B_5,k(u_z)
        Given: - p = measurement to be corrected
               - q = coeff from dist calibration
    '''
    # 1) determine bounding box to scale q_i values
    # pick upper and lower limits and compute u = ScaleToBox(qs,qmin,qmax)
    upper = 10000
    lower = 0

    mat = np.zeros(1, 3)
    for i in np.arange(start=1, stop=3, step=1):
        mat[i] = ScaleToBox(p[i], lower, upper)

    # 2) set up and solve least squares problem:
    # [F_000(u_s) ... F_555(u_s)] [[c_000^x, c_000^y, c_000^z][... ... ...][c_555^x, c_555^y, c_555^z]] = [p_s^x, p_s^y, p_s^z]
    vectCorr = np.zeros(1, 3)
    count = 0
    for i in np.arange(1, 5, 1):
        for j in np.arange(1, 5, 1):
            for k in np.arange(1, 5, 1):
                c_ijk = q[count, :]
                five = Bernie(vectCorr[0], i) * Bernie(vectCorr[1], j) * Bernie(vectCorr[2], k)
                vectCorr = vectCorr + (c_ijk * five)
                count += 1
    # return Sigma Sigma Sigma c_i,j,k B_5,i(u_x) B_5,j(u_y) B_5,k(u_z)
    return vectCorr


def main():
    # Part 1: Process data the values of C_i^expected [k] corresponding to each C_i [k] in each “frame” k of data.
    print('Enter the name of your input data file (e.g. pa2-debug-a): ')
    filename = input()
    calBodyData, calBodySize = cartesian.readInput_Body(filename + '-calbody.txt')
    calReadData, calReadSize = cartesian.readInput_Readings(filename + '-calreadings.txt')
    R = np.array([[1, 0, 0],
                  [0, math.cos(45), -math.sin(45)],
                  [0, math.sin(45), math.cos(45)]])
    F0 = cartesian.Frame(R, [1, 1, 1])
    eta0 = 1000000000000000

    N_D = calBodySize[0]
    N_A = calBodySize[1]
    N_C = calBodySize[2]
    N_framescal = calReadSize[3]

    # Position of markers with respect to the calibration body
    d = calBodyData[np.arange(0, N_D), :]
    d = d.astype(float)
    a = calBodyData[np.arange(N_D, N_A + N_D), :]
    a = a.astype(float)
    c = calBodyData[np.arange(N_D + N_A, len(calBodyData)), :]
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

    # Calculate F_D = [R_D, p_D]
    # R_D is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
    # p_D is a 3XN_frames 2D matrix, each column corresponds to the translation of a frame
    R_D = np.zeros((3, 3, N_framescal))
    p_D = np.zeros((3, N_framescal))
    for i in np.arange(0, N_framescal):
        F = icp.ICP(d, D[:, :, i], F0, eta0)
        R_i = F.get_rot()
        p_i = F.get_vec()
        R_D[:, :, i] = R_i
        print(p_D.shape)
        print(p_i.shape)
        p_D[:, i] = p_i[:, 1]

    # Calculate F_A = [R_A, p_A]
    # R_A is a 3x3xN_frames 3D matrix, each page corresponds to the rotation matrix of a frame
    # p_A is a 3XN_frames 2D matrix, each column corresponds to the translation of a frame
    R_A = np.zeros((3, 3, N_framescal))
    p_A = np.zeros((3, N_framescal))
    for i in np.arange(0, N_framescal):
        F = icp.ICP(a, A[:, :, i], F0, eta0)
        R_i = F.get_rot()
        p_i = F.get_vec()
        R_A[:, :, i] = R_i
        p_A[:, i] = p_i

    # Compute C_i expected = inv(R_D) * (R_A*ci + p_A - p_D)
    C_i = np.zeros((N_C, 3, N_framescal))
    for i in np.arange(0, N_framescal):
        for j in np.arange(0, N_C):
            C_i[j, :, i] = np.transpose(
                R_D[:, :, i].dot((R_A[:, :, i].dot(np.transpose(c[j, :]) + p_A[:, i] - p_D[:, i]))))

    # Part 2: Distortion Correction (See function)
    # truth = np.zeros(((N_C * N_framescal), 3))
    # measurements = np.zeros(((N_C * N_framescal), 3))
    # for i in np.arange(0, N_framescal):
    #     truth = C_i[:, :, i]
    #     measurements = C[:, :, i]
    #
    # # 1) determine bounding box to scale q_i values
    # # pick upper and lower limits and compute u = ScaleToBox(qs,qmin,qmax)
    # upper = 10000
    # lower = 0
    #
    # mat = np.zeros((1, 3))
    # for i in np.arange(start=1, stop=3, step=1):
    #     mat[i] = ScaleToBox(p[i], lower, upper)
    #
    # # Coefficient from distance using least squares
    # coefficient = np.zeros((216, 3))
    # for i in np.arange(0, 2):
    #     x, resnorm, residual, exitflag, output, lambda_ = scipy.linalg.lstsq(mat, truth[:, i])
    #     coefficient = [x, i]


    # Part 3: EM pivot calibration using distortion correction
    emPivotData, emPivotSize = cartesian.readInput_EmPivot(filename + '-empivot.txt')
    N_G = emPivotSize[0]
    N_framesEM = emPivotSize[1]

    # Find position of markers relative to sensor
    G = np.zeros((N_G, 3, N_framesEM))
    ind = 0
    for i in np.arange(0, N_framesEM):
        G[:, :, i] = emPivotData[np.arange(ind, ind + N_G), :]
        ind = ind + N_G

    # Use distortion correction
    G_correct = np.zeros((N_G, 3, N_framesEM))
    for i in np.arange(0, N_framesEM):
        for j in np.arange(0, N_G):
            G_correct[j, :, i] = distortionCorrection(G[j, :, i], vectCorr)

    # Define and use probe coordinate system to find g
    G2 = G_correct[:, :, 1]
    G_mid = np.mean(G2, 1)

    g = np.zeros(N_G, 3)
    for i in np.arange(N_G, 3):
        g[i, :] = G2[i, :] - G_mid

    # Pivot calibration using distortion correction
    p_tip, p_dimple = pivotCalibration.pivotCalibration(g, G_correct)

    # Part 4: Using the distortion correction and the improved pivot value, compute b_j, the locations of the
    # fiducials points with respect to the EM tracker base coordinate system
    emFiducialsData, emFiducialsSize = cartesian.readInput_EmFiducialss(filename + '-em-fiducialss.txt')
    N_G = emFiducialsSize[0]
    N_B = emFiducialsSize[1]

    # Find position of markers relative to sensor for N_B fiducials
    G = np.zeros((N_G, 3, N_B))
    ind = 0
    for i in np.arange(0, N_B):
        G[:, :, i] = emFiducialsData[np.arange(ind, ind + N_G), :]
        ind = ind + N_G

    # Use distortion correction with fiducial data
    G_correct = np.zeros((N_G, 3, N_framesEM))
    for i in np.arange(0, N_framesEM):
        for j in np.arange(0, N_G):
            G_correct[j, :, i] = distortionCorrection(G[j, :, i], vectCorr)

    # Compute b_j of fiducials using new p_tip
    R_ptr = np.zeros((3, 3, N_B))
    p_ptr = np.zeros((3, N_B))
    B = np.zeros((N_B, 3))

    for i in np.arange(0, N_B):
        F = icp.ICP(g, G[:, :, i], F0, eta0)
        R_i = F.get_rot()
        p_i = F.get_vec()
        R_ptr[:, :, i] = R_i
        p_ptr[:, i] = p_i
        B[i, :] = np.transpose(R_i * p_tip + p_i)

    # Part 5: Compute F_reg
    ctFiducialsData, ctFiducialsSize = cartesian.readInput_CtFiducials(filename + '-ct-fiducials.txt')
    F_reg = icp.ICP(B, ctFiducialsData, F0, eta0)
    R_reg = F_reg.get_rot()
    p_reg = F_reg.get_vec()

    # Part 6: Compute tip location with respect to CT image
    emNavData, emNavSize = cartesian.readInput_EmNav(filename + '-EM-nav.txt')
    N_G = emNavSize[0]
    N_framesEM = emNavSize[1]

    # Position of markers relative to N_frameEM points
    G = np.zeros((N_G, 3, N_framesEM))
    ind = 0
    for i in np.arange(0, N_framesEM):
        G[:, :, i] = emNavData[np.arange(ind, ind + N_G), :]
        ind = ind + N_G

    # Apply distortion correction to G
    G_correct = np.zeros((N_G, 3, N_framesEM))
    for i in np.arange(0, N_framesEM):
        for j in np.arange(0, N_G):
            G_correct[j, :, i] = distortionCorrection(G[j, :, i], vectCorr)

    # Compute pointer tip coordinates wrt tracker base

    R_ptr = np.zeros((3, 3, N_framesEM))
    p_ptr = np.zeros((3, N_framesEM))
    B = np.zeros((N_framesEM, 3))
    for i in np.arange(0, N_framesEM):
        F_ptr = icp.ICP(g, G_correct[:, :, i])
        R_ptr = F_ptr.get_rot()
        p_ptr = F_ptr.get_vec()
        B[i, :] = np.transpose(R_i * p_tip + p_i)

    # Compute b_j of test points with respect to EM tracker base coordinates
    # and apply F_reg to find tip locations
    b_j = np.zeros(N_framesEM, 3)
    for i in np.arange(0, N_framesEM):
        b_j[i, :] = np.transpose(R_reg * np.transpose(B[i, :]) + p_reg)

    # Save and output results
    outname = filename + '-output2.txt'
    outpath = 'outputs/' + outname
    fileID = open(outpath, 'w')
    fileID.write('%3d, %s\n' % N_framesEM, outname)
    fileID.write('%d, %d, %d\n' % b_j[0], b_j[1], b_j[2])

    return 0


main()
