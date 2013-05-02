"""op.py

This module implements

Xu et al., "Robust PCA via Outlier Pursuit", In: IEEE Trans.
on Information Theory", Vol. 58, No. 5, 2012

The problem that is solved is the following matrix
deconvolution problem:

\min_{P,C} ||P||_* + \gamma ||C||_{1,2}
s.t. || X - P + C ||_{fro} < tol

where ||C||_{1,2} is the sum of l2 norms of the columns
of C.
"""


__license__ = "Apache License, Version 2.0"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import os
import io
import sys
import time
import numpy as np
from optparse import OptionParser


def __iter_C(C, epsilon):
    """Helper for opursuit(...).
    """
    m,n=C.shape
    output = np.zeros(C.shape)
    for i in range(0,n):
        tmp = C[:,i]
        norm_tmp = np.linalg.norm(tmp, ord=2)
        if norm_tmp > epsilon:
            tmp = tmp - tmp*epsilon/norm_tmp;
        else:
            tmp = np.zeros((m,1))
        output[:,i] = tmp.ravel();
    return output


def __iter_L(L, epsilon):
    """Helper for opursuit(...).
    """
    U,S,V = np.linalg.svd(L, full_matrices=False)
    S = np.diag(S)
    for i in range(0,np.amin(S.shape)):
        if S[i,i]>epsilon :
            S[i,i] = S[i,i] - epsilon
        elif S[i,i]<-epsilon:
            S[i,i] = S[i,i] + epsilon
        else:
            S[i,i] = 0
    return np.mat(U)*np.mat(S)*np.mat(V)


def opursuit(M,O=None,gamma=None):
    """Outlier pursuit.

    Paramters
    ---------

    M : numpy matrix, shape (D, N)
        Input data, where D is the dimensionality of the data (e.g.,
        the number of pixels in an image) and N is the number of
        data samples.

    O : numpy matrix, shape (D, N), default: None
        Binary matrix, where a 1 at (i,j) signifies that the value
        has been observed, 0 signifies that this value is missing.
        In case O is 'None', all values are assumed to be existent.

    gamma : float

    Returns
    -------

    L : numpy matrix, shape (D, N)
        Low-rank approximation of M

    C : numpy matrix, shape (D, N)
        Sparse part of the deconvolution.

    term : float
        Termination value.

    iter : int
        Number of iterations it took to converge.
    """

    # sanity check - if no O is given, assume all entries are visible
    if gamma is None:
        raise Exception("\Gamma not given")
    # what is observable
    if O is None:
        O = np.ones(M.shape)

    # misc. variables
    (eta, delta, k) = (0.9, 1e-05, 0)
    # (L_{k), C_{k})
    L_cur = np.zeros(M.shape)
    C_cur = np.zeros(M.shape)
    # (L_{k-1), C_{k-1})
    L_pre = np.zeros(M.shape)
    C_pre = np.zeros(M.shape)
    # t_{k-1}, t_{k}
    t_pre = 1
    t_cur = 1
    # \mu_{k}, \bar{\mu}, for k=0
    m_cur = 0.99 * np.linalg.norm(M, ord=2)
    m_bar = delta * m_cur
    # tolerance
    tol = 1e-06 * np.linalg.norm(M, 'fro')

    stopped = False
    while not stopped:
        # note: YL, YC, GL, GC all change with k

        # lno. (3)
        YL = L_cur + (t_pre-1)/t_cur *(L_cur-L_pre)
        YC = C_cur + (t_pre-1)/t_cur *(C_cur-C_pre)

        # helper for eqs. on line (4)
        M_diff = np.multiply((YL+YC-M),O)

        # lno. (4)
        GL = YL - 0.5*M_diff
        L_new = __iter_L(GL, m_cur/2)

        # lno. (4)
        GC = YC - 0.5*M_diff
        C_new = __iter_C(GC, m_cur*gamma/2)

        # lno. (7)
        t_new = (1 + np.sqrt(4*t_cur**2+1))/2
        m_new = np.amax(eta*m_cur, m_bar)

        # check stopping crit.
        S_L = 2*(YL-L_new)+(L_new+C_new-YL-YC)
        S_C = 2*(YC-C_new)+(L_new+C_new-YL-YC)

        fro_part0 = np.linalg.norm(S_L,ord='fro')
        fro_part1 = np.linalg.norm(S_C,ord='fro')
        term_crit = fro_part0**2 + fro_part1**2

        if term_crit <= tol**2:
            stopped = True
        else:
            # L_{k-1} = L_{k}, L_{k} = L_new
            L_pre = L_cur
            L_cur = L_new
            # C_{k-1} = C_{k}, C_{k} = C_new
            C_pre = C_cur
            C_cur = C_new
            # t_{k-1} = t_{k}, t_{k} = t_new
            t_pre = t_cur
            t_cur = t_new
            # \mu_k = \mu_new
            m_cur = m_new
            k = k+1

    L = L_new
    C = C_new
    return (L, C, term_crit, k)


def main(argv=None):
    if argv is None: argv = sys.argv

    parser = OptionParser()
    parser.add_option("-i", "--dat", help="data file")
    parser.add_option("-g", "--gam", help="gamma value", type="float")
    (options, args) = parser.parse_args()

    X = np.genfromtxt(options.dat)
    t0 = time.clock()
    low_rank, sparse, term, n_iter = opursuit(X, None, options.gam)
    t1 = time.clock()
    print ("%d iterations (%.2g [sec]) to convergence (val=%.10f) on X=(%d x %d) " %
           (n_iter, (t1-t0), term, X.shape[0], X.shape[1]))

if __name__ == "__main__":
    main()
