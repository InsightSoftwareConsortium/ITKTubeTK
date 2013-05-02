"""ealm.py

Implements the exact Lagrangian multiplier approach (EALM) to solve
Candes et al.'s approach to robust PCA, i.e.,

min_{P,C} ||P||_* + gamma ||C||_1, s.t. ||M-P-C||_{fro} < eps

[1] Candes et al., "Robust Principal Component Analysis?", In:
    Journal of the ACM, Vol. 58, No. 3, 2011

The EALM algorithm was originally proposed in

[2] Lin et al., "The Augmented Lagrangian Multiplier Method for Exact Recovery
    of Corrupted Low-Rank Matrices", 2011 (available online on archivX).

This code implements Algorithm 4 of [2].
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
    parser.add_option("-s", "--sav", help="save to (rounds low-rank values)")

    # parse CLI args
    options, _ = parser.parse_args()

    # load the data
    D = np.genfromtxt(options.dat)

    # run low rank + sparse recovery
    n, m = D.shape
    low_rank, _, _ = recover(D, options.gam)

    if not options.sav is None:
        np.savetxt(options.sav, np.round(low_rank), delimiter=' ')


def recover(D, gamma=None):
    """Recover low-rank and sparse part, using Alg. 4 of [2].

    Note: gamma is lambda in Alg. 4.

    Parameters
    ---------

    D : numpy ndarray, shape (N, D)
        Input data matrix.

    gamma : float, default = None
        Weight on sparse component. If 'None', then gamma = 1/sqrt(max(D, N))
        as shown in [1] to be the optimal choice under a set of suitable
        assumptions.

    Returns
    -------

    LL : numpy array, shape (N, D)
        Low-rank part of data

    SP : numpy array, shape (N, D)
        Sparse part of data

    n_iter : int
        Number of iterations until convergence.
    """

    n, m = D.shape

    if gamma is None:
        gamma = 1/np.sqrt(np.amax([n, m]))

    # the following lines implement line 1 of Alg. 4
    Y = np.sign(D)
    l2n = np.linalg.norm(Y, ord=2)
    l2i = np.linalg.norm(np.asarray(Y).ravel(), ord=np.inf)
    dual_norm = np.amax([l2n, l2i])
    Y = Y/dual_norm

    # line 4 of Alg. 4
    A_hat = np.zeros(D.shape)
    E_hat = np.zeros(D.shape)
    D_fro = np.linalg.norm(D, ord='fro')

    # cf. section "Choosing Parameters" of [2]
    proj_tol = 1e-06*D_fro
    term_tol = 1e-07
    iter_max = 1e+03

    num_svd = 0  # track # of SVD calls
    m = 0.5/l2n  # \mu in Alg. 4
    r = 6  # \rho in Alg. 4
    sv = 5
    svp = sv

    k = 0
    converged = False
    while not converged:
        primal_converged = False
        sv = sv+np.round(n*0.1)
        primal_iter = 0

        while not primal_converged:
            # implement line 10 in Alg. 4
            T_tmp = D-A_hat+1/m*Y
            E_tmp = (np.maximum(T_tmp-gamma/m, 0) +
                     np.minimum(T_tmp+gamma/m, 0))

            # line 7 of Alg. 4
            U, S, V = np.linalg.svd(D-E_tmp+1/m*Y, full_matrices=False)
            # line 8 of Alg. 4
            svp = len(np.where(S > 1/m)[0])
            if svp < sv:
                sv = np.amin([svp+1, n])
            else:
                sv = np.amin([svp + np.round(0.05*n), n])
            A_tmp = (np.mat(U[:,0:svp]) *
                     np.diag(S[0:svp]-1/m) *
                     np.mat(V[0:svp,:]))

            # check convergence of inner optimization
            if (np.linalg.norm(A_hat-A_tmp, ord='fro') < proj_tol and
                np.linalg.norm(E_hat-E_tmp, ord='fro') < proj_tol):
                primal_converged = True

            A_hat = A_tmp
            E_hat = E_tmp
            primal_iter = primal_iter+1
            num_svd = num_svd+1

        # line 13 of Alg. 4
        Z = D-A_hat-E_hat
        Y = Y+m*Z
        m = r*m

        # evaluate stopping criteria
        stop_crit = np.linalg.norm(Z,'fro')/D_fro
        if stop_crit < term_tol:
            converged = True

        # some information about the iteration
        non_zero = len(np.where(np.asarray(np.abs(E_hat)).ravel()>0)[0])
        message = ["[iter: %.4d]" % k,
                   "#svd=%.4d" % num_svd,
                   "rank(P)=%.4d" % svp,
                   "|C|_0=%.4d" % non_zero,
                   "crit=%.4g" % stop_crit]
        print ' '.join(message)

        k = k+1
        # handle non-convergence
        if not converged and k > iter_max:
            warnings.warn("terminate after max. iter.", UserWarning)
            converged = True

    return (A_hat, E_hat, k)


if __name__ == "__main__":
    sys.exit(main())
