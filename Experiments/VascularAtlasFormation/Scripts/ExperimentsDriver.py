import sys
import os
import json
import logging
import numpy as np
from itertools import *
from optparse import OptionParser
from sklearn.cross_validation import KFold
from sklearn.cross_validation import ShuffleSplit
import ExpUtils as Utils

LOGGING_LEVELS = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG}


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = OptionParser()

    parser.add_option("", "--stage", help="Processing stage (0 = Run all)" , type="int" )
    parser.add_option("", "--dest", help="Destination base directory", default="/tmp/")
    parser.add_option("", "--data", help="Data file in JSON format (see README for format)")
    parser.add_option("", "--cvruns", help="Number of cross-validation runs (1 == single split)", type="int", default=1)
    parser.add_option("", "--config", help="Config file with relative executable paths")
    parser.add_option("", "--cells", help="Number of CVT cells to use for ATLAS building", type="int", default=1000)
    parser.add_option("", "--logat", help="Log at the specified logging level")
    parser.add_option("", "--kernelType", help="Graph kernel type (0 ... SP, 1 ... WL)", type="int", default=0)
    parser.add_option("", "--wlHeight", help="Subtree height of the WL subtree kernel", type="int", default=1)
    parser.add_option("", "--defLabel", help="Specify default labeling of graph nodes (0 ... Node ID, 1 ... Node degree)", type="int", default=0)
    parser.add_option("", "--phantomType", help="Specify the phantom type that is used (Supported are: SPL, BrainWeb)")
    parser.add_option("", "--globalLabelFile", help="Specify a global label file to use")
    parser.add_option("", "--logto", help="Log to the specified file")
    parser.add_option("", "--phantom", help="Brainweb phantom")

    (options, args) = parser.parse_args()

    # Logger configuration
    logging.basicConfig(level=LOGGING_LEVELS.get(options.logat, logging.NOTSET), filename=options.logto,
        format='%(asctime)s [%(funcName)s] %(levelname)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S')

    if (not os.path.exists(options.dest) or not os.path.isdir(options.dest)):
        print "Error: Destination directory invalid!"
        return -1
    if (not Utils.check_file(options.data)):
        print "Error: Data file not given or invalid!"
        return -1
    if (not Utils.check_file(options.config)):
        print "Error: Config file is missing!"
        return -1
    if (options.phantom is None):
        print "Error: No phantom given!"
        return -1
    if (options.phantomType is None):
        print "Error: No phantom type given!"
        return -1

    config_fid = open(options.config).read()
    config = json.loads(config_fid)

    json_fid = open(options.data).read()
    json_dat = json.loads( json_fid )

    subject_dir_list = [] # Directories with subject data
    subject_lab_list = [] # The group label for each subject, e.g., 'Male', 'Female'
    for e in json_dat["Data"]:
        subject_dir_list.append(e["Source"])
        subject_lab_list.append(e["Group"].rstrip())

    logger = logging.getLogger()

    N = len(json_dat["Data"])
    cv = ShuffleSplit(N, n_iter=options.cvruns, test_size=0.3, random_state=0)
    Utils.ensure_dir(options.dest)

    logger.debug("Destination directory = %s" % options.dest)
    logger.debug("Phantom file = %s", options.phantom)
    logger.debug("#CVT cells = %d", options.cells)
    logger.debug("Cross-validation runs = %d" % options.cvruns)
    logger.debug("#Subjects = %d", N)

    try:
        stage_opt = dict()
        stage_opt["groupLabel"] = subject_lab_list
        stage_opt["subjects"] = subject_dir_list
        stage_opt["dest"] = options.dest
        stage_opt["phantom"] = options.phantom
        stage_opt["phantomType"] = options.phantomType
        stage_opt["cells"] = options.cells
        stage_opt["kernelType"] = options.kernelType
        stage_opt["wlHeight"] = options.wlHeight
        stage_opt["defLabel"] = options.defLabel
        stage_opt["globalLabelFile"] = options.globalLabelFile
        stage_opt["randomSeed"] = 1234 # Random seed for repeatability

        if (options.stage == 1):
            stage_opt["mra_wSkull_glob"] = "*MRA.mha"            # MRA ToF images (includes skull)
            stage_opt["mri_wSkull_glob"] = "*T1-Flash.mha"       # MRI T1 images  (includes skull)
            stage_opt["mri_nSkull_glob"] = "*SkullStripped*.mha" # MRA Tof images (skull-stripped)
            Utils.compute_registrations(config, stage_opt)
        elif (options.stage == 2):
            Utils.transform_tubes_to_phantom(config, stage_opt)
        elif (options.stage == 3):
            Utils.compute_ind_atlas_edm(config, stage_opt)
        elif (options.stage > 3 and options.stage < 23):
            for cv_id,(train,test) in enumerate(cv):
                stage_opt["id"] = cv_id + 1 # Cross-validation ID
                stage_opt["trn"] = train    # Trn indices
                stage_opt["tst"] = test     # Tst indices

                res = {
                    4 : Utils.compute_grp_atlas_sum,
                    5 : Utils.compute_grp_atlas_cvt,
                    6 : Utils.compute_ind_graph_grp,
                    7 : Utils.compute_grp_graph_grp,
                    8 : Utils.compute_glo_atlas_edm,
                    9 : Utils.compute_glo_atlas_cvt,
                    10: Utils.compute_ind_graph_common,
                    11: Utils.compute_grp_graph_common,
                    12: Utils.compute_ind_graph_grp_testing,
                    13: Utils.compute_ind_graph_common_testing,
                    14: Utils.compute_ind_tube_prob_testing,
                    15: Utils.compute_ind_graph_prob_testing,
                    16: Utils.compute_trn_gk,
                    17: Utils.compute_tst_gk,
                    18: Utils.compute_full_gk,
                    19: Utils.trn_classifier,
                    20: Utils.tst_classifier,
                    21: Utils.evaluate_classifier_from_full_gk,
                    22: Utils.compute_distance_signatures
                }[options.stage](config, stage_opt)
        else:
            print "Error: Stage %d not available!" % options.stage
    except Exception as e:
        print e
        return -1


if __name__ == "__main__":
    sys.exit( main() )
