from scipy.spatial import KDTree
import numpy as np
import math
import cartesian as cart


def FindClosestPoint(fq, bnd, T):
    '''
    Given: - fq: Fn dot Qk
           - bnd: bndk
           - T: T
    Return: - c: closest points on M to Q
            - i: index of triangle mik corresponding to c
            - d: distance; || ck - Fn dot qk ||
    '''
    d, i = T.query(x=fq, k=1, distance_upper_bound=bnd)
    return T.self.data[i], i, d

def FindBestRigidTransformation(A, B):
    '''
    Method Based Off "Least Squares Rigid Motion Using SVD
    by Olga Sorkine-Hornung & Michael Rabinovich @ ETH Zurich
    https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    Given: - A
           - B
           - F
    Return: - F: Frame data type from cartesian.py containing
              best rigid transformation calculated in the method
    '''
    # compute centroids of both point sets
    a = np.sum(A) / A.shape[1]
    b = np.sum(B) / B.shape[1]

    # compute centered vectors
    X = np.empty(1)
    Y = np.empty(1)

    x = 0
    y = 0
    begin = True
    for i in A:
        x = i - a
        if begin == True:
            X = np.array([x])
            begin = False
        else:
            np.append(X, x)
    
    begin = True
    for i in B:
        y = i - b
        if begin == True:
            Y = np.array([y])
            begin = False
        else:
            np.append(Y, y)

    # compute d x d covariance matrix S = XWY^T
    W = np.ones(X.shape[1])
    xw = np.dot(X, W)
    t = np.transpose(Y)
    S = np.dot(xw, t)

    # compute the svd (gives rotation)
    u,s,vh = np.linalg.svd(S) #A = USV^H

    ye = np.eye(Y.shape[1])
    ye[Y.shape[1] - 1] = np.linalg.det(np.dot(np.transpose(vh), np.transpose(u)))
    R = np.dot(np.transpose(vh), ye)
    R = np.dot(R, np.transpose(vh))

    # compute optimal translation as t = q - Rp
    t = b - np.dot(R, a)

    F = cart.Frame(R, t)
    return F

def TerminationTest(Sigma, Epsilonmax, Epsilon):
    '''
    Given: - Sigma: array of sigma values
           - Epsilonmax: array of max epsilon values
           - Epsilon: array of epsilon values
    Return: - False (so that iterations continue) if within proper bounds
            - Otherwise, True (to end iterations)
    '''
    if 0.95 <= Epsilon[-1] <= 1 and Sigma >=0 and Epsilonmax >= 0:
        return False
    return True

'''
Iterative Closest Point registration algorithm
Given: - surface model M consisting of triangles {mi} (2d array?)
       - set of points Q = {q1, ..., qN} known to be on M
       - initial guess Fo for transformation F0 s.t. the
         points F dot qk lie on M
       - initial threshol eta0 for match closeness
Return: - c: closest point on M to Q; np.array
        - i: indices of triangles mik corresponding to c
        - d: distance; || ck - Fn dot qk ||
'''
def ICP(M, Q, F0, eta0):
    '''
    Temporary Variables: 
        - n iteration number
        - Fn = [R, p] current estimate of transformation
        - etan current match distance threshold
        - C  = {...,ck,...} - closest points on M to Q
        - D = {...,dk,...} - distances dk = || ck - Fn dot qk || (not vec)
        - I = {...,ik,...} - indices of triangles mik coreesp to ck (not vec)
        - A = {...,ak,...} - subset of Q with valid matches
        - B = {...,bk,...} - points on M corresponding to A
        - E = {...,ek,...} - residual errors bk - F dot ak
    '''

    ''' step 0 - initialization '''
    T = KDTree(M)
    n = 0
    I = np.empty(1) #
    C = np.empty(1) #
    D = 0 
    E = 0

    F = F0

    Sigma = np.empty(1)
    Epsilonmax = np.empty(1)
    Epsilon = np.empty(1)
    Eta = np.array([eta0])
    eta = eta0

    ''' step 1 - matching '''
    A = 0 
    B = 0
    
    
    terminate = False
    while terminate == False:
        k = 0
        while k < Q.shape[1]:
            bnd = 0
            c = 0
            i = 0
            d = 0
            e = 0

            if k == 0:
                # pick the distance to any point in the other cloud
                bnd = np.linalg.norm(Q[0] - M[0], 1)  # i think this is right? double check me
            else:
                bnd = np.linalg.norm(C[k - 1] - cart.frameVecProd(F, Q[k - 1]), 1)
            # develop first w simple search, later make more sophisticated using T

            [c, i, d] = FindClosestPoint(cart.frameVecProd(F, Q[k - 1]), bnd, T) 
            
            if k == 0:
                C = np.array([c])
                I = np.array([i])
                D = np.array([d])
            else:
                np.append(C, c)
                np.append(I, i)
                np.append(D, d)

            if d < eta:
                if(k == 0):
                    A = np.array(Q[k])
                    B = np.array(C[k])
                else:
                    np.append(A, Q[k])
                    np.append(B, C[k])

                # Ek - residual errors bk - F dot ak
                if k == 0:
                    E = np.array(C[k] - cart.frameVecProd(F, Q[k]))
                else:
                    np.append(E, C[k] - cart.frameVecProd(F, Q[k]))
                
            k += 1

        ''' step 2 - transformation update '''
        n += 1
        F = FindBestRigidTransformation(A, B)
        
        Edot = np.empty(1)

        ekdotsum = 0
        count = 0
        while count < E.shape[1]:
            e = np.dot(E[count], E[count])
            ekdotsum += np.dot(E[count], E[count])
            if count == 0:
                Edot = np.array(math.sqrt([e]))
            else:
                np.append(Edot, math.sqrt(e))
            
            count += 1

        # sigma =  (sqrt (sum (ek dot ek)))/ numelements(E)
        sigma = math.sqrt(ekdotsum) / Edot.shape[1]

        # epsilonmax = max(sqrt(ek dot ek))
        epsilonmax = np.amax(Edot)

        # epsilon = (sum(sqrt(ek dot ek)))/numelmenets(E)
        epsilon = np.sum(Edot) / Edot.shape[1]

        if n == 0:
            Sigma = np.array([sigma])
            Epsilonmax = np.array([epsilonmax])
            Epsilon = np.array([epsilon])
        else:
            np.append(Sigma, sigma)
            np.append(Epsilonmax, epsilonmax)
            np.append(Epsilon, epsilon)

        ''' step 3 - adjustment '''
        # compute etan from {eta0, ..., etan-1}

        # threshold can be used to restrict the influence of clearly wrong
        # matches on the computation of Fn

        # generally it should start at a fairly large value and then
        # decrease after a few iterations. one not unreasonable val 
        # might be something like 3 epsilon
        neweta = eta - (3 * epsilon)
        # if the num valid matches begins to fall significantly one
        # can increase it adaptively. too tight a bound may encourage
        # false minima
        if (3 * epsilon) > (eta/5):
            neweta = -(eta + epsilon)
        eta -= neweta
        np.append(Eta, eta)

        # also, if the mesh is incomplete, it may be adventageous to
        # exclude any matches with triangles at the edge of the mesh

        # update threshold etan and dist parameters phi

        ''' step 4 - iteration '''

        # termination test:
        terminate = TerminationTest(Sigma, Epsilonmax, Epsilon)

        n += 1
    return F
        
        









