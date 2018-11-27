"""mapcl.py

Demonstrate how to evaluate a maximum a-posteriori
graph classifier using N-fold cross-validation.
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

# Fine-structure analysis
import core.fsa as fsa
import core.utils as utils


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Setup vanilla CLI parsing and add custom arg(s).
    parser = utils.setup_cli_parsing()
    parser.add_option("",
                      "--mixComp",
                      help="number of GMM components.",
                      default=3,
                      type="int")
    (options, args) = parser.parse_args()

    # Setup logging
    utils.setup_logging(options)
    logger = logging.getLogger()

    # Read graph file list and label file list
    graph_file_list = utils.read_graph_file_list(options)
    label_file_list = utils.read_label_file_list(options, graph_file_list)

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

    # Create cross-validation folds (20% testing)
    n_graphs = len(class_info)
    cv = ShuffleSplit(n_graphs,
                      n_iter=options.cvRuns,
                      test_size=0.2,
                      random_state=0)

    # Our unique class labels
    label_set = np.unique(class_info)

    if options.normalize:
        logger.info("Running feature normalization ...")
        scaler = preprocessing.StandardScaler(copy=False)
        scaler.fit_transform(fsa_res['data_mat'])

    scores = []
    for cv_id, (trn, tst) in enumerate(cv):

        models = []
        for l in label_set:
            l_idx = np.where(class_info == l)[0]
            l_idx = np.asarray(l_idx).ravel()
            l_trn = np.intersect1d(l_idx, trn)

            pos = []
            for i in l_trn:
                tmp = np.where(fsa_res['data_idx']==i)[0]
                pos.extend(list(tmp))

            np_pos = np.asarray(pos)
            gmm_model = fsa.estimate_gm(data_mat[np_pos,:], options.mixComp)
            models.append(gmm_model)

        predict = []
        for i in tst:
            pos = np.where(data_idx==i)[0]
            map_idx = fsa.pp_gmm(data_mat[pos,:], models, argmax=True)
            predict.append(label_set[map_idx])

        # Score the MAP classifier
        truth = [class_info[i] for i in tst]
        score = accuracy_score(truth, predict)

        print "yhat :", predict
        print "gold :", truth

        logger.info("Score (%.2d): %.2f" % (cv_id, 100*score))
        scores.append(score)

    utils.show_summary(scores)


if __name__ == "__main__":
    main()
