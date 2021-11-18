import os
import numpy as np
import cartesian
from scipy.spatial.transform import Rotation as R

#implement matching - pa3

def readInput(name):
    fp = '../INPUT/' + name
    f = open(fp)
    data = np.loadtxt(f, skiprows=1, dtype=str)
    n = np.loadtxt(f, max_rows = 1, dtype=str)
    size = 1 + n[0][0]
    return size, data


def ICP(A, B):
    ''' find centroid '''
    a = np.mean(A, axis=0)
    b = np.mean(B, axis=0)

    ''' find position of point relative to centroid '''
    dist_a = np.zeros(A.shape)
    dist_b = np.zeros(B.shape)
    for i in np.arange(0, A.shape[0]):
        dist_a[i,:] = A[i,:] - a
        dist_b[i,:] = B[i,:] - b

    ''' solve for R '''
    # compute H
    H = np.zeros(3)
    for i in np.arange(0, A.shape[0]):
        a_curr = dist_a[i,:]
        b_curr = dist_b[i,:]
        h_curr = np.array(np.array(
            a_curr[0] * b_curr[0], a_curr[0] * b_curr[1], a_curr[0] * b_curr[2]
        ), np.array(
            a_curr[1] * b_curr[0], a_curr[1] * b_curr[1], a_curr[1] * b_curr[2]
        ), np.array(
            a_curr[2] * b_curr[0], a_curr[2] * b_curr[1], a_curr[2] * b_curr[2]
        ))
        H = np.concatenate(H, h_curr)

    # compute G
    delta = np.array(H[1,2] - H[2,1], 
                H[2,0] - H[0,2],
                H[0,1] - H[1,0])
    g = np.array(np.concatenate(np.trace(H), np.transpose(delta)),
        np.concatenate(delta, H + np.transpose(H) - np.trace(H) * np.eye(3)))
    
    # perform eigenvalue decomp
    v,d = np.linalg.eig(g)
    dummy, max_i = np.max(np.array(
        d[0,0], d[1,1], d[2,2], d[3,3]
    ))
    quaternion = v[:, max_i[0]]
    rot = R.from_quat(np.transpose(quaternion))

    ''' solve for p '''
    pos = np.transpose(b) - R * np.transpose(a)

    return rot, pos

def main():
    print('Enter the name of your input data file (e.g. PA3-A-Debug): ')
    filename = input()
    size_A, BAData = readInput('Problem3-BodyA.txt')
    size_B, BBData = readInput('Problem3-BodyB.txt')

    readingSize, readingData = readInput(filename + '-SampleReadingsTest.txt')

    size_S = readingSize[0]
    size_D = size_S - size_A - size_B
    numSamples = readingSize[1]

    body_A = BAData[1:size_A,:] # marker pos wrt ba coords
    body_B = BBData[1:size_B,:] # marker pos wrt bb coords

    tip_A = BAData[size_A + 1,:] # tip a pos wrt ba coords
    tip_B = BBData[size_B + 1,:] # tip b pos wrt bb coords

    ''' marker pos wrt tracker '''
    read_A = np.zeros(size_A, 3, numSamples)
    read_B = np.zeros(size_B, 3, numSamples)

    count = 0
    for i in np.arange(0, numSamples):
        read_A[:,:,i] = readingData[count:count + size_A - 1,:]
        count += size_A
        read_B[:,:,i] = readingData[count:count + size_B - 1,:]
        count += size_B + size_D
        
    
    ''' compute frame A '''
    rot_A = np.zeros(3,3,numSamples)
    pos_A = np.zeros(3,numSamples)
    for i in np.arange(0, numSamples):
        rot, pos = ICP(body_A, read_A[:,:,i])
        rot_A[:,:,i] = rot
        pos_A[:,i] = pos
    F_A = cartesian.Frame(rot_A, pos_A)


    ''' compute frame B '''
    rot_B = np.zeros(3,3,numSamples)
    pos_B = np.zeros(3,numSamples)
    for i in np.arange(0, numSamples):
        rot, pos = ICP(body_B, read_B[:,:,i])
        rot_B[:,:,i] = rot
        pos_B[:,i] = pos
    F_B = cartesian.Frame(rot_B, pos_B)

    ''' compute dk '''
    '''d = np.zeros(numSamples, 3)
    for i in numSamples:
        d[i,:] = '''
    d = cartesian.frameVecProd(cartesian.frameProd(np.inv(F_B), F_A), tip_A)
    
    '''output the CT coordinates ck corresponding
        to each sample taken (ck = F dot dk) 
        for now, assume F = I => d is our answer'''

    '''save and output results'''
    out = str(output of CT coordinates)
    outname = filename + '-Output.txt'
    os.makedirs('OUTPUT', mode=0o777, exist_ok=False)
    outpath = 'OUTPUT/' + outname
    fileID = open(outpath, 'w+')
    fileID.write('%d, %s\n' % (numSamples, outname))
    fileID.write(out) # This should be whatever the output is 'output = [A, B, norm of distances]
                      # should be 15x7 or 20x7
    fileID.close()

    return 0

main()