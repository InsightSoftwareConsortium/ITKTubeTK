"""utils.py

Convenience functions to build applications that use pyfsa.
"""


__license__ = "Apache License, Version 2.0 (see TubeTK)"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import os
import sys
import logging
import numpy as np
from optparse import OptionParser


# Define the log levels
LOGGING_LEVELS = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG}


def _radii_callback(option,opt,value,parser):
    radii = []
    for e in value.split(','):
        radii.append(int(e))
    parser.values.radii = radii


def _globalLabelFile_callback(option,opt,value,parser):
    print value
    if value is None:
        parser.values.globalLabelFile = None
    if not os.path.exists(value):
        raise Exception("Global label file %s does not exist!" % value)
    parser.values.globalLabelFile = value


def _classInfo_callback(option,opt,value,parser):
    if not os.path.exists(value):
        raise Exception("Class info file %s does not exist!" % value)
    parser.values.classInfo = value


def _groupInfo_callback(option,opt,value,parser):
    if not os.path.exists(value):
        raise Exception("Group info file %s does not exist!" % value)
    parser.values.groupInfo = value


def _logLevel_callback(option,opt,value,parser):
    if not value in LOGGING_LEVELS:
        raise Exception("Level : %s not supported!" % value)
    parser.values.logLevel = value


def _graphList_callback(option,opt,value,parser):
    if not os.path.exists(value):
        raise Exception("Graph list file %s does not exist!" % value)
    parser.values.graphList = value


def _baseDir_callback(options,opt,value,parser):
    if not os.path.exists(value):
        raise Exception("Directory %s does not exist!" % value)
    parser.values.baseDir = value


def setup_cli_parsing():
    """Setup CLI parser for common CLI arguments.

    Returns
    -------

    parser : OptionParser object
    """

    parser = OptionParser()
    parser.add_option("",
                      "--labelPrefix",
                      help="prefix of the label files.")
    parser.add_option("",
                      "--skip",
                      default=0,
                      type="int",
                      help="Skip N header entries.")
    parser.add_option("",
                      "--omitDegenerate",
                      action="store_true",
                      default=False,
                      help="Omit single vertex subgraphs.")
    parser.add_option("",
                      "--writeAs",
                      default="/tmp/data",
                      help="feature file base name.")
    parser.add_option("",
                      "--seed",
                      type="int",
                      help="Seed for random number generator.")
    parser.add_option("",
                      "--normalize",
                      action="store_true",
                      default=False,
                      help="Enable feature normalization (mean/var).")
    parser.add_option("",
                      "--cvRuns",
                      default=5,
                      type="int",
                      help="number of cross-validations to run.")
    parser.add_option("",
                      "--recompute",
                      action="store_true",
                      default=False,
                      help="Force feature recomputation.")
    parser.add_option("",
                      "--graphList",
                      type="string",
                      action="callback",
                      callback=_graphList_callback,
                      help="input file with graph adj. file names.")
    parser.add_option("",
                      "--baseDir",
                      type="string",
                      action="callback",
                      callback=_baseDir_callback,
                      help="base directory of graph adj. file names.")
    parser.add_option("", "--logLevel",
                      type="string",
                      default="debug",
                      action="callback",
                      callback=_logLevel_callback,
                      help="set logging level.")
    parser.add_option("",
                      "--classInfo",
                      type="string",
                      action='callback',
                      callback=_classInfo_callback,
                      help="class information file.")
    parser.add_option("",
                      "--groupInfo",
                      type="string",
                      action='callback',
                      callback=_groupInfo_callback,
                      help="group information file.")
    parser.add_option("",
                      "--logTo",
                      help="Specify logging file.")
    parser.add_option("",
                      "--globalLabelFile",
                      type="string",
                      action="callback",
                      callback=_globalLabelFile_callback,
                      help="global label file for all graphs.")
    parser.add_option("", "--radii",
                      type="string",
                      action='callback',
                      callback=_radii_callback,
                      help="list of neighborhood radi(i), e.g., 1,2,3")
    return parser


def setup_logging(options):
    """Configure logger.

    Defines the logging format.
    """

    log_format = '%(asctime)s [%(funcName)s] %(levelname)s: %(message)s'
    logging.basicConfig(level=LOGGING_LEVELS.get(options.logLevel,
                                                 logging.NOTSET),
                        filename=options.logTo,
                        format=log_format,
                        datefmt='%Y-%m-%d %H:%M:%S')


def remap_labels(labels):
    """Map labels in {c_1,...,c_C} to {0,...,C-1}.

    Parameters
    ----------

    labels : numpy array, shape (L,)
        Input labels (could also be 'chars')

    Returns
    -------

    labels : numpy array, shape (L,)
        Labels in {0,...C-1} where C is the number of unique
        labels.
    """

    for i,l in enumerate(np.unique(labels)):
        labels[labels == l] = i
    return labels


def read_graph_file_list(options):
    """Create a list of graph file names.

    Parameters
    ----------

    options : object returned by parse_args() of OptionParser
        CLI options.

    Returns
    -------

    graph_files : list, len=N
        List of N graph file names, contained in the file
        that was specified as the value to '--graphList'.
        In case '--labelPrefix' is provided, the value to
        that option is used as a prefix to each graph file
        name.
    """

    base = options.baseDir
    if not base is None:
        return [os.path.join(base, l.strip()) for l in open(options.graphList)]
    return [l.split() for l in open(options.graphList)]


def read_label_file_list(options, graph_file_list):
    """Create a list of vertex label files (one for each graph).

    Paramters
    ---------

    options : object returned by parse_args() of OptionParser
        CLI options.

    graph_file_list : list
        List of N graph files.

    Returns
    -------

    label_files : list len=N
        List of label files (graph file name + label prefix).
    """

    if options.labelPrefix is None:
        raise Exception("No label prefix given!")

    N = len(graph_file_list)
    label_files = [None for i in range(0, N)]
    for i in range(0, N):
        label_files[i] = "%s.%s" % (graph_file_list[i], options.labelPrefix)
    return label_files


def read_class_info(options):
    """Read class label information from file.

    The class label file needs to follow the convention that
    there is only one label per line, e.g.,

    0
    1
    0
    1
    ...

    Parameters
    ----------

    options : object returned by parse_args() of OptionParser
        CLI options.

    Returns
    -------

    class_labels : list of 'int', len=N
        List of labels (one label per graph).
    """

    return [int(l.strip()) for l in open(options.classInfo)]


def read_group_info(options):
    """Read grouping info from file (in the form of attribute indicators).

    The format of the group info file is a binary N x G matrix, where a
    non-zero entry at the j-th entry of the n-th row signifies that this
    the j-th attribute is present.

    Example
    -------

    0 1 0 0
    0 0 0 1
    0 1 0 1

    Parameters
    ----------

    options : object returned by parse_args() of OptionParser
        CLI options.

    Returns
    -------

    attributes : numpy matrix, shape (N, G)
        Binary attribute indicator matrix.
    """

    grp_info = np.asmatrix(np.genfromtxt(options.groupInfo, dtype="int"))
    if grp_info.shape[0] == 1:
        grp_info = grp_info.transpose()
    return grp_info


def show_summary(scores):
    """Print a classification report to stdout.

    Parameters
    ----------

    scores : numpy array, shape (K,)
        Classifier scores.
    """
    logger = logging.getLogger()

    logger.info("Avg(scores) : %.2f" % (np.mean(scores)*100))
    logger.info("Std(scores) : %.2f" % (np.std(scores)*100))
