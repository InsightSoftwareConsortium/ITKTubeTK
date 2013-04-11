"""dcl.py

Demonstrate evaluation of a discriminant SVM graph classifier
using N-Fold cross-validation.
"""

__license__ = "Apache License, Version 2.0 (see TubeTK)"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


# Graph handling
import networkx as nx

# Machine learning
from sklearn.metrics import accuracy_score
from sklearn.cross_validation import KFold
from sklearn.cross_validation import ShuffleSplit
from sklearn.linear_model import LogisticRegression
from sklearn.grid_search import GridSearchCV
from sklearn import preprocessing
from sklearn import svm

# Misc.
from optparse import OptionParser
import logging
import numpy as np
import scipy.sparse
import time
import sys
import os

# pyfsa imports
import core.fsa as fsa
import core.utils as utils


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Setup vanilla CLI parsing and add custom arg(s).
    parser = utils.setup_cli_parsing()
    parser.add_option("",
                      "--codewords",
                      help="number of codewords.",
                      default=50,
                      type="int")
    (options, args) = parser.parse_args()

    # Setup logging
    utils.setup_logging(options)
    logger = logging.getLogger()

    # Read graph file list and label file list
    graph_file_list = utils.read_graph_file_list(options)
    if not options.globalLabelFile is None:
        label_file_list = [options.globalLabelFile] * len(graph_file_list)
    else:
        label_file_list = utils.read_label_file_list(options,
                                                     graph_file_list)

    # Read class info and grouping info
    class_info = utils.read_class_info(options)
    group_info = utils.read_group_info(options)

    assert (group_info.shape[0] ==
            len(class_info) ==
            len(graph_file_list) ==
            len(label_file_list))

    # Zip lists together
    data = zip(graph_file_list,
               label_file_list,
               class_info)

    # Run fine-structure analysis
    fsa_res = fsa.run_fsa(data,
                          options.radii,
                          options.recompute,
                          options.writeAs,
                          options.skip,
                          options.omitDegenerate)
    data_mat = fsa_res['data_mat']
    data_idx = fsa_res['data_idx']

    # Create cross-validation folds
    n_graphs = len(class_info)
    cv = ShuffleSplit(n_graphs,
                      n_iter=options.cvRuns,
                      test_size=0.2,
                      random_state=0)

    # Try inplace feature normalization
    if options.normalize:
        logger.info("Running feature normalization ...")
        scaler = preprocessing.StandardScaler(copy=False)
        scaler.fit_transform(fsa_res['data_mat'])

    scores = []
    for cv_id, (trn, tst) in enumerate(cv):

        # Compose training data
        pos = []
        for i in trn:
            tmp = np.where(data_idx==i)[0]
            pos.extend(list(tmp))
        np_pos = np.array(pos)

        # Learn a codebook from training data
        codebook = fsa.learn_codebook(data_mat[np_pos,:],
                                      options.codewords,
                                      options.seed)

        # Compute BoW histograms for training data
        bow_trn_mat = np.zeros((len(trn), options.codewords))
        for cnt, i in enumerate(trn):
            np_pos = np.where(data_idx==i)[0]
            bow_trn_mat[cnt,:] = np.asarray(fsa.bow(data_mat[np_pos,:],
                                                    codebook))

        # Cross-validate (5-fold) SVM classifier and parameters
        param_selection = [{'kernel': ['rbf'],
                            'gamma': np.logspace(-6,2,10),
                            'C': [1, 10, 100, 1000]},
                           {'kernel': ['linear'],
                            'C': [1, 10, 100, 1000]}]
        clf = GridSearchCV(svm.SVC(C=1), param_selection)
        clf.fit(bow_trn_mat, np.asarray(class_info)[trn], cv=5)

        # Compute BoW histograms for testing data
        bow_tst_mat = np.zeros((len(tst), options.codewords))
        for cnt,i in enumerate(tst):
            pos =  np.where(data_idx==i)[0]
            bow_tst_mat[cnt,:] = fsa.bow(data_mat[pos,:], codebook)

        print "yhat : ", clf.predict(bow_tst_mat)
        print "gold : ", np.asarray(class_info)[tst]

        # Score the classifier
        score = clf.score(bow_tst_mat, np.asarray(class_info)[tst])
        scores.append(score)

        logger.info("Score (%.2d): %.2f" % (cv_id,100*score))

    utils.show_summary(scores)


if __name__ == "__main__":
    main()
