from scipy.spatial import KDTree
import numpy as np
import math
import cartesian as cart
from collections.abc import Iterable


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
    if(isinstance(d, Iterable)):
        d = d[0]
    return T.data[i - 1], i, d


def FindBestRigidTransformation(A, B):
    '''
    Method Based Off "Least Squares Rigid Motion Using SVD
    by Olga Sorkine-Hornung & Michael Rabinovich @ ETH Zurich
    https://igl.ethz.ch/projects/ARAP/svd_rot.pdf
    Given: - A
           - B
    Return: - F: Frame data type from cartesian.py containing
              best rigid transformation calculated in the method
    '''
    # inputs should be of size Nx3
    # print(A.shape)
    # assert(A.shape[1] == 3)
    # assert(B.shape[1] == 3)
    # assert(len(A) == len(B))

    # calc size
    N = A.shape[0]

    # compute centroids of both point sets
    a = np.mean(A, axis=0)
    b = np.mean(B, axis=0)

    # compute centered vectors
    X = A - np.tile(a, (N, 1))
    Y = B - np.tile(b, (N, 1))

    # compute d x d covariance matrix S = XWY^T
    # X and Y should be dxn where n and d are both 3
    cov = np.dot(np.transpose(X), Y)

    # compute the svd
    u, s, vh = np.linalg.svd(cov)  # A = USV^H

    # calculate rotation
    R = vh * u

    if np.linalg.det(R) < 0:
        vh[2, :] *= -1 # fix if neg
        R = vh.T * u.T

    # compute optimal translation 
    r1 = (a[0]*R[0][0]) + (a[1]*R[0][1]) + (a[2]*R[0][2])
    r2 = (a[0]*R[1][0]) + (a[1]*R[1][1]) + (a[2]*R[1][2])
    r3 = (a[0]*R[2][0]) + (a[1]*R[2][1]) + (a[2]*R[2][2])
    sub = np.array([r1, r2, r3])
    t = b - sub

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
    if 0.95 <= Epsilon[-1] <= 1 and Sigma[-1] >= 0 and Epsilonmax[-1] >= 0:
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
    I = np.empty(1)  #
    C = np.empty(1)  #
    D = np.empty(1)
    E = np.empty(1)

    F = F0

    Sigma = np.empty(1)
    Epsilonmax = np.empty(1)
    Epsilon = np.empty(1)
    Eta = np.array([eta0])
    eta = eta0

    ''' step 1 - matching '''
    A = np.empty(1)
    B = np.empty(1)

    c = 0
    i = 0
    d = 0
    e = 0
    terminate = False
    while terminate == False:
        k = 0
        done = False
        while k < Q.shape[0]:
            bnd = 0
            if k == 0:
                # pick the distance to any point in the other cloud
                M[0] = np.array(list(map(np.float, M[0])))
                Q = Q.astype(np.float)
                M = M.astype(np.float)
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
                C = np.concatenate((C, [c]), axis=0)
                I = np.concatenate((I, [i]), axis=0)
                D = np.concatenate((D, [d]), axis=0)
            if d < eta:
                if (done == False):
                    A = np.array([Q[k]])
                    B = np.array([C[k]])
                    done = True
                else:
                    A = np.concatenate((A, [Q[k]]), axis=0)
                    B = np.concatenate((B, [C[k]]), axis=0)

                # Ek - residual errors bk - F dot ak
                if k == 0:
                    E = np.array(C[k] - cart.frameVecProd(F, Q[k]))
                else:
                    E = np.concatenate((E, C[k] - cart.frameVecProd(F, Q[k])), axis=0)

            k += 1

        ''' step 2 - transformation update '''
        n += 1
        F = FindBestRigidTransformation(A, B)

        Edot = np.empty(1)

        ekdotsum = 0
        count = 0
        while count < E.shape[0]:
            e = np.dot(E[count], E[count])
            ekdotsum += np.dot(E[count], E[count])
            if count == 0:
                Edot = np.array([math.sqrt(e)])
            else:
                earr = np.array([math.sqrt(e)])
                Edot = np.concatenate((Edot, earr), axis=0)

            count += 1

        # sigma =  (sqrt (sum (ek dot ek)))/ numelements(E)
        sigma = math.sqrt(ekdotsum) / Edot.shape[0]

        # epsilonmax = max(sqrt(ek dot ek))
        epsilonmax = np.amax(Edot)

        # epsilon = (sum(sqrt(ek dot ek)))/numelmenets(E)
        epsilon = np.sum(Edot) / Edot.shape[0]

        if n == 0:
            Sigma = np.array([sigma])
            Epsilonmax = np.array([epsilonmax])
            Epsilon = np.array([epsilon])
        else:
            sigarr = np.array([sigma])
            epsmaxarr = np.array([epsilonmax])
            epsarr = np.array([epsilon])

            Sigma = np.concatenate((Sigma, sigarr), axis=0)
            Epsilonmax = np.concatenate((Epsilonmax, epsmaxarr), axis=0)
            Epsilon = np.concatenate((Epsilon, epsarr), axis=0)

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
        if (3 * epsilon) > (eta / 5):
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
