##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################


"""Permutation testing.
"""


__license__ = "Apache License, Version 2.0"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import os
import sys
import json
import math
import random
import numpy as np

from sklearn import svm
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from itertools import permutations
from optparse import OptionParser


def permutation_testing(bas,ncv,npt,trnb,tstb):
    """

    Performance evaluation and permutation testing
    with SVM classifiers that use precomp. kernels.

    Implementation of test 1, as proposed in

    Ojala, M. and G. Garriga, "Permutation Tests for Studying Classifier
    Performance", in JMLR 11 (2010)

    :param bas: Base directory for experiments.
    :param ncv: #Crossvalidations.
    :param npt: #Permutation tests.
    :param trnb: Basename for trn. data.
    :param tstb: Basename for tst. data.
    """

    A_stat_star = []    # Test statistic
    A_stat_true = []    # True labels (for each A_stat_star)
    stat_cv_perm = []   # CV-wise statistic under H_0 (for each permutation)
    stat_cv_true = []   # CV-wise true labels (for each permutation)

    total_nr_sup = 0     # Nr. of support vectors (SV)
    total_nr_trn = 0     # Nr. of training samples

    # Iterate over all cross-validation runs
    random.seed(1234)
    for cv in xrange(1,ncv+1):

        trn_info_fname = os.path.join(bas,
            "cv-%.4d" % cv,
            "trn-Common.json")
        tst_info_fname = os.path.join(bas,
            "cv-%.4d" % cv,
            "tst-Common.json")

        # Read training/testing labels
        fid_trn = open(trn_info_fname).read()
        fid_tst = open(tst_info_fname).read()

        trn_lab = json.loads(fid_trn)["labels"]
        tst_lab = json.loads(fid_tst)["labels"]

        N = len(trn_lab)
        T = len(tst_lab)

        trn_kern_fname = os.path.join(bas,
            "cv-%.4d" % cv,
            "trn-Common.kern.bin")
        tst_kern_fname = os.path.join(bas,
            "cv-%.4d" % cv,
            "tst-Common.kern.bin")

        # Load training/testing kernels
        trn_kern = np.fromfile(trn_kern_fname,
                               dtype="double").reshape(N,N)
        tst_kern = np.fromfile(tst_kern_fname,
                               dtype="double").reshape(T,N)


        # Train/test classifier to get predictions
        clf = svm.SVC(kernel="precomputed",C=1)
        clf.fit(trn_kern,trn_lab)
        lab_hat = clf.predict(tst_kern)

        # Show confusion matrices
        cm = confusion_matrix(lab_hat, tst_lab)
        print "Confusion matrix for %d-th CV run" % cv
        print cm

        # Keep track of nr. of SVs / nr. training samples
        total_nr_sup = total_nr_sup + sum(clf.n_support_)
        total_nr_trn = total_nr_trn + N

        # Store predictions (just add to the lists)
        A_stat_star += list(lab_hat)
        A_stat_true += list(tst_lab)

        # Compute the classifier predictions under random permutations
        B_stat_perm = []
        B_stat_true = []
        for cnt in range(0,int(npt)):
            random.shuffle(trn_lab)

            # Train / Test SVM with permuted labels
            clf = svm.SVC(kernel="precomputed",C=1)
            clf.fit(trn_kern,trn_lab)
            lab_hat = clf.predict(tst_kern)

            # Store predictions / true labels as seperate lists
            B_stat_perm.append(lab_hat)
            B_stat_true.append(tst_lab)

        stat_cv_perm.append(B_stat_perm)
        stat_cv_true.append(B_stat_true)

    print "SV Fraction: %.2f (%d)" % (total_nr_sup / float(total_nr_trn), total_nr_sup)
    return (A_stat_star,A_stat_true,stat_cv_perm,stat_cv_true)


def f1_score(tp,fp,tn,fn):
    """ Computes F1-score, see http://en.wikipedia.org/wiki/F1_score

    :param tp: True positives  (TP)
    :param fp: False positives (FP)
    :param tn: True negatives  (TN)
    :param fn: False negatives (FN)
    :return: F1-score in [0,1]
    """
    return 2*tp/float(2*tp + fn + fp)


def accuracy(tp,fp,tn,fn):
    """ Compute binary classifier accuracy.

    :param tp: True positives  (TP)
    :param fp: False positives (FP)
    :param tn: True negatives  (TN)
    :param fn: False negatives (FN)
    :return: Classifier accuracy in [0,1]
    """
    return (tp+tn)/float(tp+fp+tn+fn)


def split_confusion_matrix(cm):
    """ Split confusion matrix into it's parts.

    :param cm: 2x2 numpy confusion matrix
    :return: (TP,FP,TN,FN), see http://en.wikipedia.org/wiki/Confusion_matrix
    """

    assert cm.shape == (2,2), "Only binary classifiers supported!"

    tp = cm[0,0]
    fp = cm[1,0]
    tn = cm[1,1]
    fn = cm[0,1]
    return (tp,fp,tn,fn)


def compute_pvalues(T_star_pred,
                    T_star_true,
                    T_perm_pred,
                    T_perm_true,
                    n_permutations):
    """ Compute p-values for classifier accuracy and T1-score.

    This function computes p-values for the two binary classifier performance
    measures accuracy and T1-score. The function is designed to return p-values
    for the situation where T* is obtained by averaging over multiple CV runs.
    Consequently, the permutation test needs to be designed to obtain a p-value
    for the average CV accuracy T*. The p-value is estimated as

    p-pvalue = (#{T_perm > T*}+1) / (n_permutations + 1)

    :param T_star_pred: List of label predictions for T*.
    :param T_star_true: List of ground truth labels for T*.
    :param T_perm_pred: List of CV lists of label predictions for T_perm.
    :param T_perm_true: List of CV lists of ground truth labels for T_perm.
    :return: Tuple (x,y) of p-values for accuracy and T1-score.
    """

    # Compute confusion matrix for the test statistic
    cm = confusion_matrix(T_star_true,T_star_pred)
    assert cm.shape == (2,2), 'Only binary classifiers supported!'
    (tp,fp,tn,fn) = split_confusion_matrix(cm)

    T_acc = accuracy(tp,fp,tn,fn) # T* for accuracy
    T_f1s = f1_score(tp,fp,tn,fn) # T* for F1 score

    assert len(T_perm_pred) == len(T_perm_true), 'Size error!'

    global_cm = [] # CMs for each CV run
    for i,pred_list in enumerate(T_perm_pred):
        cv_cm = [] # CMs for this CV run
        for j,pred in enumerate(pred_list):
            truth = T_perm_true[i][j]
            cm = confusion_matrix(truth,pred)
            cv_cm.append(cm)
        global_cm.append(cv_cm)

    n_cv = len(global_cm) # Number of CV runs
    n_pt = n_permutations # Number of permutations

    tp = np.zeros((n_cv,n_pt))
    fp = np.zeros((n_cv,n_pt))
    tn = np.zeros((n_cv,n_pt))
    fn = np.zeros((n_cv,n_pt))
    for i,tmp in enumerate(global_cm):
        cv_cm = global_cm[i]
        tp[i,] = [split_confusion_matrix(c)[0] for c in cv_cm]
        fp[i,] = [split_confusion_matrix(c)[1] for c in cv_cm]
        tn[i,] = [split_confusion_matrix(c)[2] for c in cv_cm]
        fn[i,] = [split_confusion_matrix(c)[3] for c in cv_cm]
    tp_sum = np.sum(tp,axis=0) # Sum of TPs over all CV runs
    fp_sum = np.sum(fp,axis=0) # Sum of FPs over all CV runs
    tn_sum = np.sum(tn,axis=0) # Sum of TNs over all CV runs
    fn_sum = np.sum(fn,axis=0) # Sum of FNs over all CV runs

    # Compute T_perm statistics
    (T_perm_acc, T_perm_f1s) = ([],[])
    for i in xrange(0,n_pt):
        T_perm_acc.append(accuracy(tp_sum[i],fp_sum[i],tn_sum[i],fn_sum[i]))
        T_perm_f1s.append(f1_score(tp_sum[i],fp_sum[i],tn_sum[i],fn_sum[i]))

    # Compute the probability of obtaining values as extreme as T*
    f1s_pvalue = (sum(T_perm_f1s >= T_f1s)+1) / float(n_pt+1)
    acc_pvalue = (sum(T_perm_acc >= T_acc)+1) / float(n_pt+1)

    result = dict()
    result["accuracy"] = {"stat" : T_acc, "pval" : acc_pvalue }
    result["f1_score"] = {"stat" : T_f1s, "pval" : f1s_pvalue }
    return result


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = OptionParser()
    parser.add_option("", "--dir", help="Experiments directory.")
    parser.add_option("", "--ncv", help="Number of CV splits.", type="int")
    parser.add_option("", "--npt", help="Number of perm. tests",type="int")
    parser.add_option("", "--trnb",help="Basename for trn. data.")
    parser.add_option("", "--tstb",help="Basename for tst. data.")
    (options,args) = parser.parse_args()

    stat = permutation_testing(options.dir,
                               options.ncv,
                               options.npt,
                               options.trnb,
                               options.tstb)

    T_star_pred = stat[0]
    T_star_true = stat[1]
    T_perm_pred = stat[2]
    T_perm_true = stat[3]

    res = compute_pvalues(T_star_pred,
                          T_star_true,
                          T_perm_pred,
                          T_perm_true,
                          options.npt)

    for key in res:
        stat = res[key]["stat"]
        pval = res[key]["pval"]
        print "%s (nperm=%d): %.2f (%.4f)" % (key,options.npt,stat,pval)


if __name__ == "__main__":
    sys.exit(main())
