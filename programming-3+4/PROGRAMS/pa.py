import numpy as np
import cartesian

#implement matching - pa3

def readInput(name):
    fp = '../INPUT/' + name
    f = open(fp)
    data = np.loadtxt(f, skiprows=1, dtype=str)
    n = np.loadtxt(f, max_rows = 1, dtype=str)
    size = 1 + n[0][0]
    return size, data


def ICP(body, reading):
    return 0, 0

def main():
    print('Enter the name of your input data file (e.g. pa3-a-debug): ')
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
    for i in numSamples:
        read_A[:,:,i] = readingData[count:count + size_A - 1,:]
        count += size_A
        read_B[:,:,i] = readingData[count:count + size_B - 1,:]
        count += size_B + size_D
        
    
    ''' compute frame A '''
    rot_A = np.zeros(3,3,numSamples)
    pos_A = np.zeros(3,numSamples)
    for i in numSamples:
        rot, pos = ICP(body_A, read_A[:,:,i])
        rot_A[:,:,i] = rot
        pos_A[:,i] = pos
    F_A = cartesian.Frame(rot_A, pos_A)


    ''' compute frame B '''
    rot_B = np.zeros(3,3,numSamples)
    pos_B = np.zeros(3,numSamples)
    for i in numSamples:
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
    return 0

main()