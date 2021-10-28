# 3D point set to 3D point set registration algorithm
# let's do ICP
from scipy.spatial import KDTree
import numpy as np
import cartesian as cart

    # Fn = [R, p] - current estimate of transformation
    # current match distance threshold
    # C = {...,ck,...} - closest points on M to Q
    # D = {...,dk,...} - distances dk = || ck - Fn dot qk || (not vec)
    # I = {...,ik,...} - indices of triangles mik coreesp to ck (not vec)
    # A = {...,ak,...} - subset of Q with valid matches
    # B = {...,bk,...} - points on M corresponding to A
    # E = {...,ek,...} - residual errors bk - F dot ak

def ICP(M, Q, F0, eta0):
    '''
    Iterative Closest Point registration algorithm
    input: M - surface model consisting of triangles {mi}
           Q - set of points {q1, ..., qN} known to be on M
           F0 - initial guess for transformation F0 s.t. the
                points F dot qk lie on M
           eta0 - initial threshold for match closeness
    output: probability density function of the relative pose between maps
    '''

    '''
    step 0: initialization
        - input surface model M and points Q

        - build an appropriate data structure to facilitate
          finding the closest point matching search
    '''
    T = KDTree(M)

    n = 0 #iteration counter
    N = Q.size() # max itrs

    F = F0

    eta = 100000000000000000000000000 #large number, effectively infinity
    
    I = np.ones(N) # {...,1,...}
    C = np.zeros(N) # {...,point on m1,...}
    D = np.zeros(N) # {...,||ck - F0 dot qk||,...}

    '''
    step 1: matching
    '''
    A = np.zeros(N)
    acount = 0
    B = np.zeros(N)
    bcount = 0

    k = 1
    while(k < N):
        bndvect = cart.frameVecProd(F, Q[k]) - C[k - 1] #??
        bnd = np.linalg.norm(bndvect, 1) # length of the vector [Fn dot qk - ck]
        c, i, d = FindClosestPoint(cart.frameVecProd(F, Q[k]),C[k],I[k],T) # come back and fill in
        if(d < eta):
            A[acount] = Q[k]
            B[bcount] = c
            acount +=1
            bcount += 1
        else:
            #test if prob(qk ~ ck) > etan
            k=k
        k += 1
    
    '''
    step 2: transformation update
    '''
    n += 1




    return 0

def FindClosestPoint(Fdotq,c,i,T): 
    return 0, 0, 0 #come back and fill in

