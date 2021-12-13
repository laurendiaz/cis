import os
import numpy as np
import cartesian
from scipy.spatial.transform import Rotation as R
import hyperspy as hs
import math

''' Functions created for PA3 '''

def get_mesh():
    '''
    Function to import mesh from input file and store in appropriate
    data structures
    input: none
    output: verts_arr - array of vertices from mesh
            indices - indices of vertices from mesh
            num_tris - total number of triangles in mesh
    '''
    sur_dict = hs.io_plugins.sur.file_reader(os.getcwd() + '/INPUT/Problem4MeshFile.sur')
    print(sur_dict)
    mesh_data = sur_dict  # .get("data")
    num_verts = mesh_data[0]
    verts_arr = np.transpose(np.reshape(mesh_data[2:num_verts * 3 + 1], np.array(3, num_verts)))
    num_tris = mesh_data[num_verts * 3 + 2]
    indices = np.reshape(mesh_data[num_verts * 3 + 3:-1], np.array(6, num_tris))
    indices = indices[:, 1:3]
    return verts_arr, indices, num_tris


def readInput(name):
    '''
    Function to read input from given file
    input: name - name of the input file
    output: size - number of elements in file
            data - data from input file
    '''
    fp = os.getcwd() + '/INPUT/' + name
    f = open(fp)
    n = np.loadtxt(f, max_rows=1, dtype=str)
    data = np.loadtxt(f, dtype=float)
    size = 1 + int(n[0])
    f.close()
    return size, data


def readInputSampleReadings(name):
    fp = os.getcwd() + '/INPUT/' + name
    f = open(fp)
    n = np.loadtxt(f, max_rows=1, dtype=str, delimiter=',')
    data = np.loadtxt(f, dtype=str, delimiter=',')
    numLEDSInFrame = int(n[0])
    numSampleFrames = int(n[1])
    f.close()
    return numLEDSInFrame, numSampleFrames, data


def ProjOnSeg(c, p, q):
    l = np.dot(c - p, q - p) / np.dot(q - p, q - p)
    if (l < 0):
        l = 0
    elif (l > 1):
        l = 1
    c = p + l * (q - p)
    return c


def FindClosestNeighbor(a, triangle):
    p = triangle[:, 0]
    q = triangle[:, 1]
    r = triangle[:, 2]

    ''' find params for interp using least squares '''
    left = np.array(np.array(q - p, 1), np.array(r - p, 1))
    right = np.array(a - p, 1)
    l, m = np.linalg.lstsq(left, right)
    c = p + l * (q - p) + m * (r - p)
    c_curr = c

    if (l + m > 1):
        c = ProjOnSeg(c_curr, q, r)
    elif (l < 0):
        c = ProjOnSeg(c_curr, r, p)
    elif (m < 0):
        c = ProjOnSeg(c_curr, p, q)

    dist = abs(c - a)
    return c, dist


def Reg(A, B):
    ''' find centroid '''
    a = np.mean(A, axis=0)
    b = np.mean(B, axis=0)

    ''' find position of point relative to centroid '''
    dist_a = np.zeros(A.shape)
    dist_b = np.zeros(B.shape)
    for i in np.arange(0, A.shape[0]):
        dist_a[i, :] = A[i, :] - a
        dist_b[i, :] = B[i, :] - b

    ''' solve for R '''
    # compute H
    H = np.zeros((3, 3))
    for i in np.arange(0, A.shape[0]):
        a_curr = dist_a[i, :]
        b_curr = dist_b[i, :]
        h_curr = np.array([np.array(
            [a_curr[0] * b_curr[0], a_curr[0] * b_curr[1], a_curr[0] * b_curr[2]]
        ), np.array(
            [a_curr[1] * b_curr[0], a_curr[1] * b_curr[1], a_curr[1] * b_curr[2]]
        ), np.array(
            [a_curr[2] * b_curr[0], a_curr[2] * b_curr[1], a_curr[2] * b_curr[2]]
        )])
        H = np.add(H, h_curr)

    print(H)

    # compute G
    delta = np.array([H[1, 2] - H[2, 1],
                      H[2, 0] - H[0, 2],
                      H[0, 1] - H[1, 0]])
    delta = np.reshape(delta, (3, 1))
    test = np.trace(H)
    # this doesn't work
    '''
    g = np.vstack((np.append((np.trace(H), np.transpose(delta))),
                   np.append((delta, H + np.transpose(H) - np.trace(H) * np.eye(3)))))
                   '''
    g = 0

    # perform eigenvalue decomposition
    v, d = np.linalg.eig(g)
    dummy, max_i = np.max(np.array(
        d[0, 0], d[1, 1], d[2, 2], d[3, 3]
    ))
    quaternion = v[:, max_i[0]]
    rot = R.from_quat(np.transpose(quaternion))

    ''' solve for p '''
    pos = np.transpose(b) - R * np.transpose(a)

    return rot, pos


''' Functions created for PA4 '''

def termTest(sigma, Emax, Eavg):
    '''
    Test for termination of ICP
    input: sigma
           Emax
           Eavg
    output: ret - True or False corresponding to whether or not the termination
                  test has determined it fit for the ICP algorithm to stop
                  iteration
    '''
    i = sigma.shape[0]
    ret = True
    if sigma[i] <= 0.01 or Emax[i] <= 0.01:
        print('Condition 1')
    else:
        if i > 1 and (Eavg[i] / Eavg[i - 1]) <= 1 and (Eavg[i] / Eavg[i - 1]) >= 0.95:
            print('Condition 2')
        else:
            ret = False
    return ret



def ICP(numSamples, d):
    '''
    Implementation of Iterative Closest Point Algorithm 
    input: numSamples - number of samples total
           d - 
    output: Rn - rotation matrix of result of ICP
            pn - position vector of result of ICP
    '''
    verts_arr, indices, num_tris = get_mesh()

    ''' init transformation '''
    Rn = np.eye(3)
    pn = np.array([0, 0, 0])
    D = np.zeros(numSamples)[0]
    thresh = 100000
    c = np.zeros(numSamples)[2]

    sigma = []
    Emax = []
    Eavg = []

    iter = True

    while (iter):
        # Find point pairs w bounding box
        A = np.zeros(numSamples, 3)
        B = np.zeros(numSamples, 3)
        E = np.zeros(numSamples, 1)
        D = np.zeros(numSamples, 1)

        for i in np.arange(1, numSamples):
            s = np.transpose(Rn * np.transpose(d[i, :]) + pn)
            dist = np.linalg.norm(s - verts_arr[indices[0, 0] + 1, :])

            for triangle in np.arange(1, num_tris):
                # get vert coords
                tri = indices[triangle, :]
                p = verts_arr[tri[0] + 1, :]
                q = verts_arr[tri[1] + 1, :]
                r = verts_arr[tri[2] + 1, :]

                c = 0

                # calc lower and upper bounds
                lower_x = np.min(np.array([p[0], q[0], r[0]]))
                lower_y = np.min(np.array([p[1], q[1], r[1]]))
                lower_z = np.min(np.array([p[2], q[2], r[2]]))

                upper_x = np.max(np.array([p[0], q[0], r[0]]))
                upper_y = np.max(np.array([p[0], q[0], r[0]]))
                upper_z = np.max(np.array([p[0], q[0], r[0]]))

                if (lower_x - dist <= s[0] <= upper_x + dist <= s[1] <= upper_y + dist and
                        lower_z - dist <= s[2] <= upper_z + dist):
                    hi, disti = FindClosestNeighbor(np.transpose(s),
                                                       np.array(np.transpose(p), np.transpose(q), np.transpose(r)))
                    if np.linalg.norm(disti) < dist:
                        dist = np.linalg.norm(disti)
                        c = hi
            if (dist < thresh):
                A[i, :] = s
                B[i, :] = np.transpose(c)
                D[i] = dist

        # update transformation
        Rn, pn = Reg(A, B)
        edot = np.zeros(numSamples, 1)
        for i in np.arange(1, numSamples):
            s = np.transpose(np.transpose(Rn * d[i, :]) + pn)
            e = B[i, :] - s
            edot[i] = np.dot(e, e)
            A[i, :] = s
            E[i] = np.linalg.norm(e)

        sigma = np.array([sigma, np.mean(edot)])
        Emax = np.array([Emax, np.max(math.sqrt(edot))])
        Eavg = np.array([Eavg, np.mean(math.sqrt(edot))])

        # adjust
        thresh = 5 * np.mean(E[i])

        # term condition
        if (termTest(sigma, Emax, Eavg)):
            iter = False

    return Rn, pn