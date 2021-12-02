import pa
import numpy as np
import cartesian
import math
import os


def termTest(sigma, Emax, Eavg):
    i = sigma.shape[0]
    ret = True
    if sigma[i] <= 0.01 or Emax[i] <= 0.01:
        print('Condition 1')
    else:
        if i > 1 and (Eavg[i] / Eavg[i - 1]) <= 1 and (Eavg[i] / Eavg[i - 1]) >= 0.98:
            print('Condition 2')
        else:
            ret = False
    return ret


'''
Full ICP
'''


def ICP(numSamples, d, filename):
    verts_arr, indices, num_tris = pa.get_mesh()

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
                    hi, disti = pa.FindClosestNeighbor(np.transpose(s),
                                                       np.array(np.transpose(p), np.transpose(q), np.transpose(r)))
                    if np.linalg.norm(disti) < dist:
                        dist = np.linalg.norm(disti)
                        c = hi
            if (dist < thresh):
                A[i, :] = s
                B[i, :] = np.transpose(c)
                D[i] = dist

        # update transformation
        Rn, pn = pa.Reg(A, B)
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

    # generate output file
    exp = np.array([A, B, E])
    avgerr = np.mean(E)
    maxerr = np.max(E)
    stderr = np.std(E)

    outname = filename + '-Output.txt'
    os.makedirs('OUTPUT', mode=0o777, exist_ok=False)
    outpath = 'OUTPUT/' + outname
    fileID = open(outpath, 'w+')
    fileID.write('%3d, %s\n' % (numSamples, outname))

    print(avgerr, maxerr, stderr)

    for i in np.arange(0, numSamples):
        fileID.write('%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n', exp[i][0], exp[i][1], exp[i][2], exp[i][0],
                     exp[i][1], exp[i][2])
    fileID.close()


def main():
    ''' import data '''
    print('Enter the name of your input data file (e.g. PA4-A-Debug): ')
    filename = input()
    size_A, BAData = pa.readInput('Problem4-BodyA.txt')
    size_B, BBData = pa.readInput('Problem4-BodyB.txt')

    size_S, numSamples, readingData = pa.readInputSampleReadings(filename + '-SampleReadingsTest.txt')

    size_D = size_S - size_A - size_B

    body_A = BAData[1:size_A, :]  # marker pos wrt ba coords
    body_B = BBData[1:size_B, :]  # marker pos wrt bb coords

    tip_A = BAData[size_A - 1, :]  # tip a pos wrt ba coords
    tip_B = BBData[size_B - 1, :]  # tip b pos wrt bb coords

    ''' marker pos wrt tracker '''
    read_A = np.zeros((size_A, 3, numSamples))
    read_B = np.zeros((size_B, 3, numSamples))

    count = 0
    for i in np.arange(0, numSamples):
        read_A[:, :, i] = readingData[count:count + size_A, :]
        count += size_A
        read_B[:, :, i] = readingData[count:count + size_B, :]
        count += size_B + size_D

    ''' compute frame A '''
    rot_A = np.zeros((3, 3, numSamples))
    pos_A = np.zeros((3, numSamples))
    for i in np.arange(0, numSamples):
        rot, pos = pa.Reg(body_A, read_A[:, :, i])
        rot_A[:, :, i] = rot
        pos_A[:, i] = pos
    F_A = cartesian.Frame(rot_A, pos_A)

    ''' compute frame B '''
    rot_B = np.zeros((3, 3, numSamples))
    pos_B = np.zeros((3, numSamples))
    for i in np.arange(0, numSamples):
        rot, pos = pa.Reg(body_B, read_B[:, :, i])
        rot_B[:, :, i] = rot
        pos_B[:, i] = pos
    F_B = cartesian.Frame(rot_B, pos_B)

    ''' compute dk '''
    '''d = np.zeros(numSamples, 3)
    for i in numSamples:
        d[i,:] = '''
    d = cartesian.frameVecProd(cartesian.frameProd(np.inv(F_B), F_A), tip_A)

    # find points c_k on surface mesh that are closest to s_k
    verts_arr, indices, num_tris = pa.get_mesh()
    c = []
    distances = []
    for j in np.arange(0, numSamples):
        for i in np.arange(0, num_tris):
            triangle = indices[i, :]
            p = verts_arr[triangle[0] + 1, :]
            q = verts_arr[triangle[1] + 1, :]
            r = verts_arr[triangle[2] + 1, :]
            c_j, dist = pa.FindClosestNeighbor(np.transpose(d[j, :]),
                                               np.array(np.transpose(p),
                                                        np.transpose(q),
                                                        np.transpose(r)))
            c[j] = c_j
            distances[j] = dist

    ICP(numSamples, d, filename)


main()
