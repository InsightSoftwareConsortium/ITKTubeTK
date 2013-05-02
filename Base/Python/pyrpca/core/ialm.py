"""ialm.py

Implements the inexact Lagrangian multiplier approach (IALM) to solve
the matrix deconvolution problem

\min_{P,C} ||P||_* + \gamma ||C||_1, s.t. ||M-P-C||_{fro} < eps

that was proposed as an approach to solve Candes et al.'s robust
PCA formulation, cf.,

[1] Candes et al., "Robust Principal Component Analysis?",
    In: Journal of the ACM, Vol. 58, No. 3, 2011

Reference for the IALM algorithm:

[2] Lin et al., "The Augmented Lagrangian Multiplier Method
    for Exact Recovery of Corrupted Low-Rank Matrices", 2011
    (available online on archivX).
"""


__license__ = "Apache License, Version 2.0"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import os
import sys
import time
import numpy as np
from optparse import OptionParser
import warnings


def main(argv=None):
    """ Functionality to test the module from the commane line.
    """

    if argv is None:
        argv = sys.argv

    parser = OptionParser()
    parser.add_option("-i", "--dat", help="data file")
    parser.add_option("-g", "--gam", help="gamma value", type="float")
    parser.add_option("-s", "--sav", help="save to")

    (options, args) = parser.parse_args()

    t0 = time.clock()
    X = np.genfromtxt(options.dat)
    t1 = time.clock()
    print "data loading took %.2g [sec]" % (t1-t0)

    m, n = X.shape
    print "data (%d x %d) loaded in %.2g [sec]" % (m, n, (t1-t0))
    low_rank, sparse, n_iter = recover(X, options.gam)

    if not options.sav is None:
        # save rounded low-rank result - usefull for images
        np.savetxt(options.sav, np.round(low_rank), delimiter=' ')


def recover(D, gamma=None):
    """Recover low-rank and sparse part.

    Paramters
    ---------

    D : numpy ndarray, shape (N, D)
        Input data matrix.

    gamma : float, default = None
        Weight on sparse component. If 'None', then gamma = 1/sqrt(max(D,N))
        as shown in [1] to be the optimal choice under a set of suitable
        assumptions.

    Returns
    -------

    LL : numpy array, shape (N,D)
        Low-rank part of data

    SP : numpy array, shape (N,D)
        Sparse part of data

    n_iter : int
        Number of iterations until convergence.
    """

    m, n = D.shape
    Y=D

    if gamma is None:
        gamma = 1/np.sqrt(np.amax([m, n]))

    tol = 1e-07
    max_iter = 1000

    # l2n is the (matrix) 2-norm, i.e., the max singular value
    l2n = np.linalg.norm(Y, ord=2)
    # l2i is the 'inf' norm, i.e., the max. abs. value of Y
    l2i = np.linalg.norm(np.asarray(Y).ravel(), ord=np.inf)

    # computes J(Y) = max(||Y||_2,\gamma^{-1}||Y||_{\inf}), cf. [2], eq. (10)
    dual_norm = np.amax([l2n, l2i/gamma])
    # computes line 1 in Alg. 5 of [2]
    Y = Y/dual_norm

    A_hat = np.zeros(D.shape)
    E_hat = np.zeros(D.shape)

    # cf. section "Choosing Parameters" of [2]
    m = 1.25/l2n
    m_b = m*1e07
    rho = 1.5
    sv = 10
    # Frobenius norm of the original matrix
    D_fro = np.linalg.norm(D, ord='fro')

    k=0
    total_svd=0
    converged=False
    while not converged:
        # start timer for k-th iteration
        t0 = time.clock()

        # part of line 4 of Alg. 5
        tmp = D-A_hat+(1/m)*Y
        E_hat = np.maximum(np.asmatrix(tmp)-gamma/m,0)
        E_hat = E_hat + np.minimum(np.asmatrix(tmp)+gamma/m,0)

        # line 4 of Alg. 5
        U, S, V = np.linalg.svd(D-E_hat+(1/m)*Y, full_matrices=False)

        svp = len(np.where(S>1/m)[0])
        if svp < sv:
            sv = np.amin([svp+1, n])
        else:
            sv = np.amin([svp+np.round(0.05*n), n])

        A_hat = np.mat(U[:,0:svp]) * np.diag(S[0:svp]-1/m) * np.mat(V[0:svp,:])
        total_svd = total_svd+1

        # line 8 of Alg. 5
        Z = D-A_hat-E_hat
        Y = Y+m*Z

        # lno. (9) of Alg. 5 in [2]
        m = np.amin([m*rho, m_b])

        t1 = time.clock()

        # eval. stopping crit.
        conv_crit = np.linalg.norm(Z, ord='fro')/D_fro
        if conv_crit < tol:
            converged = True

        if k % 10 == 0:
            rank_est = np.linalg.matrix_rank(A_hat)
            non_zero = len(np.where(np.asarray(np.abs(E_hat)).ravel()>0)[0])
            print ("[iter: %.4d]: rank(P) = %.4d, |C|_0 = %.4d, crit=%.10f" %
                    (k, rank_est, non_zero, conv_crit))

        # handle non-convergence
        k = k+1
        if not converged and k>max_iter:
            warnings.warn("terminate after max. iter.", UserWarning)
            converged = True

    return (A_hat, E_hat, k)


if __name__ == "__main__":
    main()
