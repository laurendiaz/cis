import pa
import numpy as np
import cartesian
import math
import os
from ICP_Reg_util import *


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
