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
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################


"""Experiment utilities.
"""


__license__ = "Apache License, Version 2.0"
__author__  = "Roland Kwitt, Kitware Inc., 2013"
__email__   = "E-Mail: roland.kwitt@kitware.com"
__status__  = "Development"


import os
import sys
import json
import glob
import cPickle
import logging
import random
import subprocess
import numpy as np
from sklearn import svm
from sklearn import cross_validation as cval


def check_file(file_name):
    """ Checks if file exists.

    :param file_name: Name of the file to check.
    :returns: True if file exists, False otherwise.
    """

    if file is None:
        return False
    return os.path.exists(file_name) and os.path.isfile(file_name)


def ensure_dir(directory):
    """ Ensures that a directory exists (if necessary, it creates it).

    :param directory: Path to the directory to check (or create).
    """

    if not os.path.exists(directory):
        print "Creating directory %s" % directory
        os.makedirs(directory)


def check_dir(directory):
    """ Checks existence of a directory.

    :param directory: Path to the directory to check.
    :returns: True, if directory exists, False otherwise.
    """

    return os.path.exists(directory)


def compute_registrations(config, options):
    """ Performs all necessary registrations.

    Runs all image registrations to compute the transform from MRA space to
    phantom space. This is not included in the cross-validation experiments,
    since we only need to compute that once.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger();

    for subject_dir in options["subjects"]:
        subject_id = int(os.path.basename(subject_dir).rsplit('-')[-1])

        mra_wSkull_listing = glob.glob(os.path.join(subject_dir, "Normal%.3d%s" %(subject_id, options["mra_wSkull_glob"])))
        mri_wSkull_listing = glob.glob(os.path.join(subject_dir, "Normal%.3d%s" %(subject_id, options["mri_wSkull_glob"])))
        mri_nSkull_listing = glob.glob(os.path.join(subject_dir, options["mri_nSkull_glob"]))

        sum_len = 0
        sum_len += len( mri_wSkull_listing )
        sum_len += len( mra_wSkull_listing )
        sum_len += len( mri_nSkull_listing )

        if not sum_len == 3:
            raise Exception( "Mismatch in nr. of MRI/MRA files for subject %s!" % subject_id )

        ensure_dir(os.path.join(options["dest"], os.path.basename(subject_dir)))

        rig_transform_file = os.path.join(options["dest"], os.path.basename(subject_dir), "MRA-T1.tfm" )
        aff_transform_file = os.path.join(options["dest"], os.path.basename(subject_dir), "T1-Phantom.tfm")
        rig_image_file = os.path.join(options["dest"], os.path.basename(subject_dir), "mra-t1.mha")
        aff_image_file = os.path.join(options["dest"], os.path.basename(subject_dir), "t1-phantom.mha")

        param_type = "--registration Rigid"
        param_save = "--saveTransform %s" % rig_transform_file
        param_fimg = "%s" % mra_wSkull_listing[0]
        param_mimg = "%s" % mri_wSkull_listing[0]
        param_samp =  "--resampledImage %s" % rig_image_file

        # MRA to T1 registration (RIGID)
        logger.debug("Running rigid registration to compute %s -> %s transform",
            os.path.basename(param_fimg), os.path.basename(param_mimg))

        cmd_reg0 = [ config["Exec"]["RegisterImages"],
            param_type, param_save, param_fimg, param_mimg, param_samp]
        print cmd_reg0
        subprocess.call(cmd_reg0)

        param_type = "--registration Affine"
        param_iter = "--affineMaxIterations %d" % 200
        param_save = "--saveTransform %s" % aff_transform_file
        param_fimg = "%s" % mri_nSkull_listing[0]
        param_mimg = "%s" % options["phantom"]
        param_samp =  "--resampledImage %s" % aff_image_file

        # T1 to Phantom registration (AFFINE)
        logger.debug("Running affine registration to compute %s -> %s transform",
            os.path.basename(param_fimg), os.path.basename(param_mimg))

        cmd_reg1 = [config["Exec"]["RegisterImages"],
            param_type, param_save, param_fimg, param_mimg, param_iter, param_samp]
        print cmd_reg1
        subprocess.call(cmd_reg1)


def transform_tubes_to_phantom(config, options):
    """ Map tubes to phantom space.

    Transform tubes (i.e., vascular networks) into phantom space, using
    the two image registration transforms, i.e., from MRA to T1 and from
    T1 to Phantom.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    for subject_dir in options["subjects"]:
        tubes_original_file = os.path.join(subject_dir, "VascularNetwork.tre" )

        if (not check_file(tubes_original_file)):
            raise Exception("Original tube file %s not existent!" % tubes_original_file)

        # Output tube files (to be generated in subjects destination directory)
        tubes_mra_t1_file = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-mra-t1.tre" )
        tubes_t1_pha_file = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-t1-phantom.tre")

        # Transform files (assumed to exist in destination directory of subject)
        mra_t1_tfm_file = os.path.join(options["dest"], os.path.basename(subject_dir), "MRA-T1.tfm" )
        t1_pha_tfm_file = os.path.join(options["dest"], os.path.basename(subject_dir), "T1-Phantom.tfm" )

        if (not check_file(mra_t1_tfm_file)):
            raise Exception("MRA to T1 transform file %s missing!" % mra_t1_tfm_file)
        if (not check_file(t1_pha_tfm_file)):
            raise Exception("T1 to Phantom transform file %s missing!" % t1_pha_tfm_file)

        # Tubes (in MRA image space) to MRI T1 image space
        logger.debug("MRA to T1 space for %s (using %s transform)"
            % (tubes_original_file, mra_t1_tfm_file))

        t0_cmd = [config["Exec"]["TubeTransform"],
                  tubes_original_file,
                  tubes_mra_t1_file,
                  "--transformFile %s" % mra_t1_tfm_file]
        subprocess.call(t0_cmd)

        # Tubes (in MRI T1 image space) to Phantom image space
        logger.debug("T1 to Phantom space for %s (using %s transform)"
            % (tubes_mra_t1_file, tubes_t1_pha_file))

        t1_cmd = [config["Exec"]["TubeTransform"],
                  tubes_mra_t1_file,
                  tubes_t1_pha_file,
                  "--transformFile %s" %  t1_pha_tfm_file]
        subprocess.call(t1_cmd)


def compute_ind_atlas_edm(config, options):
    """ Compute Euclidean Distance Map (EDM) for all individuals.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    for subject_dir in options["subjects"]:
        den_image_name = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-t1-phantomDD.mha")
        rad_image_name = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-t1-phantomRad.mha")
        tan_image_name = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-t1-phantomTan.mha")

        # Tubes in phantom space (assumed to exist at that stage)
        tubes_in_phantom_space = os.path.join(options["dest"], os.path.basename(subject_dir), "VascularNetwork-t1-phantom.tre" )
        if (not check_file(tubes_in_phantom_space)):
            raise Exception( "Tube file %s missing!" % tubes_in_phantom_space )

        logger.debug("Create individual EDM for %s"
            % os.path.basename(subject_dir))

        cmd = [config["Exec"]["TubesToDensityImage"],
            tubes_in_phantom_space, den_image_name, rad_image_name, tan_image_name,
            "--inputTemplateImage %s" % options["phantom"],
            "--useSquareDistance"]
        subprocess.call(cmd)


def compute_grp_atlas_sum(config, options):
    """ Compute a group-wise Euclidean Distance Map (EDM).

    Take all EDMs (one per individual) and form a group-wise EDM for
    each group the individuals belong to, e.g., Male/Female.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    ensure_dir(cv_dir) # Create if necessary

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:
        grp_dir = os.path.join(cv_dir, grp) # e.g., <cvdir>/Female
        ensure_dir(grp_dir) # Create if necessary

        pop_cnt = 0
        for i in options["trn"]:
            if not options["groupLabel"][i] == grp:
                continue
            pop_cnt += 1

        atlas_edm_doc = os.path.join(grp_dir, "NormalDD.odc")
        fid_atlas_edm_doc = open(atlas_edm_doc, "w")
        fid_atlas_edm_doc.write('NumberOfObjects = %d\n' % pop_cnt)

        for i in options["trn"]:
            if not options["groupLabel"][i] == grp:
                continue

            ind_edm_file = os.path.join(options["dest"], os.path.basename(options["subjects"][i]),"VascularNetwork-t1-phantomDD.mha")
            if not check_file(ind_edm_file):
                raise Exception("Euclidean distance file %s missing!" % ind_edm_file)

            fid_atlas_edm_doc.write('Type = Image\n')
            fid_atlas_edm_doc.write('Name = %s\n' % ind_edm_file)
            fid_atlas_edm_doc.write('EndObject =\n')

        fid_atlas_edm_doc.close()

        atlas_edm_file = os.path.join(grp_dir, "NormalDD-mean.mha")
        atlas_edm_var_file = os.path.join(grp_dir, "NormalDD-var.mha")

        logger.debug("Compute ATLAS EDM file %s (for group %s)"
            % (atlas_edm_file, grp))

        cmd = None
        if options["phantomType"] == "BrainWeb":
            cmd = [config["Exec"]["AtlasBuilderUsingIntensity"],
                atlas_edm_doc, atlas_edm_file, atlas_edm_var_file,
                "--outputSize %d,%d,%d" % (181, 217, 181),
                "--outputSpacing %d,%d,%d" % (1,1,1)]
        elif options["phantomType"] == "SPL":
            cmd = [config["Exec"]["AtlasBuilderUsingIntensity"],
                atlas_edm_doc, atlas_edm_file, atlas_edm_var_file,
                "--outputSize %d,%d,%d" % (256, 256, 256),
                "--outputSpacing %d,%d,%d" % (1,1,1)]
        else:
            raise Exception("Phantom type %s not supported!" % options["phantomType"])

        print cmd
        subprocess.call(cmd)


def compute_grp_atlas_cvt(config, options):
    """ Compute Central Voronoi Tessellation (CVT) on group-wise EDMs.

    Computes a Central Voronoi Tessellation (CVT) from the group-wise
    EDMs. The CVT represents the space tessellation upon which we
    eventually will define the graph node connectivity.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:
        logger.debug("Currently processing group %s" % grp)

        grp_dir = os.path.join(cv_dir, grp)
        if not check_dir(grp_dir):
            raise Exception("Cross-validation group directory %s missing!" % grp_dir)

        atlas_edm_org_file = os.path.join(grp_dir, "NormalDD-mean.mha")  # Has to exist!
        atlas_edm_thr_file = os.path.join(grp_dir, "bNormalDD-mean.mha") # Is created!

        if (not check_file(atlas_edm_org_file)):
            raise Exception("Mean EDM file %s missing!" % atlas_edm_org_file)

        logger.debug("Create ATLAS EDM file %s (through intensity windowing)"
            % atlas_edm_thr_file)

        cmd = [config["Exec"]["ImageMath"],
            atlas_edm_org_file,
            "-i", "1917", "2047", "0",  "100",
            "-W", "1", atlas_edm_thr_file]
        subprocess.call(cmd)

        # Run CVT
        atlas_cvt_map_file = os.path.join(grp_dir, "bNormalDD-mean.cvt" )
        atlas_cvt_img_file = os.path.join(grp_dir, "bNormalDD-mean.cvt.mha" )

        logger.debug("Create ATLAS CVT map/image %s, %s (with CVT cells = %d)"
            % (atlas_cvt_map_file, atlas_cvt_img_file, options["cells"]))

        cmd = [config["Exec"]["ImageMath"], atlas_edm_thr_file,
            "-S", "%d" % options["randomSeed"], # RNG seed for reproducable results
            "-Z", "%d" % options["cells"], "1000", "100000", atlas_cvt_map_file,
            "-W", "1", atlas_cvt_img_file]
        subprocess.call(cmd)


def graph_from_atlas(config, options, grp, atlas_type, data_type):
    """ Compute graph using a populations atlas.

    For each subject in our dataset, we take 1) the tubes (in phantom space) and
    2) a group atlas (e.g., Male/Female/Common) and computes a graph representation
    of the tubes.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    :param grp: Selects the group for which we want to compute individual graphs.
    :param atlas_type: Selects atlas to use for graph generation.
    :param data_type: Selects the type of data to use, i.e., training (trn), testing (tst).
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation base directory %s missing!" % cv_dir)


    # This specifies the ATLAS file to use (i.e., the CVT tessellation)
    atlas_file = os.path.join(cv_dir, atlas_type, "bNormalDD-mean.cvt.mha")
    if not check_file(atlas_file):
        raise Exception("Atlas %s missing!" % atlas_file)


    for i in options[data_type]: # i.e., trn/tst
        if not options["groupLabel"][i] == grp:
            continue

        subject_dir = os.path.join(options["dest"], os.path.basename(options["subjects"][i]))
        tubes_in_phantom_space = os.path.join(subject_dir, "VascularNetwork-t1-phantom.tre" )
        logger.debug("Using tube file: %s" % tubes_in_phantom_space)

        if (not check_file(tubes_in_phantom_space)):
            raise Exception("Transformed tubes (%s) are missing!" % tubes_in_phantom_space)

        # Directory that contains subject-specific graphs (depend on CVT) for each cross-validation
        cv_subject_dir = os.path.join(cv_dir, os.path.basename(options["subjects"][i]))
        ensure_dir(cv_subject_dir)

        logger.debug("Compute graph for %s individual %s (using ATLAS type %s)"
            % (grp, os.path.basename(options["subjects"][i]), atlas_type))

        output_graph_file = os.path.join(cv_subject_dir, "VascularNetwork-t1-phantom-%s.grp" % atlas_type )
        cmd = [config["Exec"]["TubeToTubeGraph"],
            tubes_in_phantom_space, atlas_file, output_graph_file]
        subprocess.call(cmd)


def compute_ind_graph_grp(config, options):
    """ Driver to compute a graph representation of each individual from a group atlas.

    E.g., for each individual in groups Male and Female, we compute

        - one graph for subject A using the Male atlas
        - one graph for subject A using the Female atlas

    Note: Uses only training examples!

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_atlas(config, options, grp, atlas_selector, "trn")


def graph_from_graphs(config, options, grp, atlas_type, create_image=False):
    """ Summarize a set of graphs into one graph.

    Composes a list (stored in an .odc file) with graph names for a given
    group and uses this list to compute a summary graph. The graphs in the
    list are compiled based on their group membership and the atlas that
    was used to create the graphs.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    :param grp: Selects the group for which we want to compute a summary graph.
    :param atlas_type: Selectes the individual graphs based on the atlas that
                       was used to generate them.
    :param create_image: Run TubeGraphToImage for visualization of results.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    pop_cnt = 0
    for i in options["trn"]:
        if not options["groupLabel"][i] == grp:
            continue
        pop_cnt += 1

    graph_doc_file = os.path.join(cv_dir, grp, "NormalGraph-%s.odc" % atlas_type)
    fid = open( graph_doc_file, "w" )
    fid.write( "NumberOfObjects = %d\n" % pop_cnt)

    for i in options["trn"]:
        if not options["groupLabel"][i] == grp:
            continue

        # Assume existence of tubes (in Phantom space) for each subject
        cv_subject_dir = os.path.join(cv_dir, os.path.basename(options["subjects"][i]))
        tube_graph = os.path.join(cv_subject_dir, "VascularNetwork-t1-phantom-%s.grp" % atlas_type)
        if (not check_file(tube_graph)):
            raise Exception("Graph %s missing!" % tube_graph)

        fid.write( "Type = Image\n" )
        fid.write( "Name = %s\n" % tube_graph )
        fid.write( "EndObject = \n" )

    fid.close()

    cmd = [config["Exec"]["MergeTubeGraphs"],
        graph_doc_file,
        os.path.join(cv_dir, grp, "Normal-mean-%s.grp" % atlas_type),
        "%d" % options["cells"]]
    subprocess.call(cmd)

    # Create a summary graph image
    if create_image:
        atlas_cvt_file = os.path.join(cv_dir, atlas_type, "bNormalDD-mean.cvt")

        logger.debug("Create summary graph image (for group %s) using ATLAS %s"
            % (grp, atlas_cvt_file))

        if not check_file(atlas_cvt_file):
            raise Exception("Atlas CVT file %s missing!" % atlas_cvt_file)

        cmd = [config["Exec"]["TubeGraphToImage"],
            atlas_cvt_file,
            os.path.join(cv_dir, grp, "Normal-mean-%s.grp" % atlas_type),
            os.path.join(cv_dir, grp, "Normal-mean-%s.cvt" % atlas_type)]
        subprocess.call(cmd)


def compute_grp_graph_grp(config, options):
    """ Driver to compute one representative graph for each group (for different atlases).

    E.g., for Female, Male group, compute

        - one Male graph using Female atlas
        - one Male graph using Male atlas
        - one Female graph using Male atlas
        - one Female graph using Female atlas

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_graphs(config, options, grp, atlas_selector, False)


def compute_glo_atlas_edm(config, options):
    """ Compute a global (i.e., Common) EDM from the group-wise EDMs.

    Computes a (glo)bal Euclidean Distance Map for all groups,
    based on the group-specific EDms. Each distance map per group is
    equally weighted.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    common_dir = os.path.join(cv_dir, "Common")
    ensure_dir(common_dir)

    weight = 1.0/len(group_names)  # Equally weigh group DD maps

    dd_list = []
    for grp in group_names:
        grp_dd_file =  os.path.join(cv_dir, grp, "bNormalDD-mean.mha")
        if not check_file(grp_dd_file):
            raise Exception( "Group EDM file %s missing!" % grp_dd_file)
        grp_dd_file_weighted = os.path.join(cv_dir, grp, "bNormalDD-mean-weighted.mha")
        cmd = [config["Exec"]["ImageMath"],
            grp_dd_file,
            "-a", str(weight/2.0), str(weight/2.0),
            grp_dd_file, "-W", "1",
            grp_dd_file_weighted]
        subprocess.call(cmd)
        dd_list.append(grp_dd_file_weighted)

    dd_common_file = os.path.join(common_dir, "bNormalDD-mean.mha" )
    for i in range(len(dd_list)):
        (infile_a, infile_b, out_file) = ("", "", "")
        if i == 0:
            continue
        if i == 1:
            infile_a = dd_list[i-1]
            infile_b = dd_list[i]
            out_file = dd_common_file
        else:
            infile_a = dd_common_file
            infile_b = dd_list[i]
            out_file = dd_common_file
        cmd = [config["Exec"]["ImageMath"],
            infile_a,
            "-a", "1.0", "1.0",
            infile_b, "-W", "1",
            out_file]
        subprocess.call(cmd)


def compute_glo_atlas_cvt(config, options):
    """ Compute a global (i.e., Common) CVT from the global (i.e., Common) EDM.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    glo_dd_cvt_map_file = os.path.join(cv_dir, "Common", "bNormalDD-mean.cvt" )
    glo_dd_cvt_img_file = os.path.join(cv_dir, "Common", "bNormalDD-mean.cvt.mha" )

    common_dir = os.path.join(cv_dir, "Common")
    ensure_dir(common_dir)

    glo_dd_edm_file = os.path.join(common_dir, "bNormalDD-mean.mha")
    cmd = [config["Exec"]["ImageMath"], glo_dd_edm_file,
        "-S", "%d" % options["randomSeed"],
        "-Z", "%d" % options["cells"], "1000", "100000", glo_dd_cvt_map_file,
        "-W", "1", glo_dd_cvt_img_file]
    subprocess.call(cmd)


def compute_ind_graph_common(config, options):
    """ Driver to compute a graph representation of an individual using the global atlas.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    common_dir = os.path.join(cv_dir, "Common")
    if not check_dir(common_dir):
        raise Exception("Directory %s missing!" % common_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    # Runs graph computation using different ATLASes
    for grp in group_names:
        logger.info("Compute graphs (using Common atlas) for group %s!" % grp)
        graph_from_atlas(config, options, grp, "Common", "trn")


def compute_grp_graph_common(config, options):
    """ Driver to compute a summary graph for each group, based on the graphs
    constructed from the global (i.e., Common) atlas.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    common_dir = os.path.join(cv_dir, "Common")
    if not check_dir(common_dir):
        raise Exception("Directory %s missing!" % common_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:
        graph_from_graphs(config, options, grp, "Common", create_image=True)


def compute_ind_graph_grp_testing(config, options):
    """ Driver to compute a graph representation of an individual using the
    group-wise atlases (FOR TESTING DATA!).

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_atlas(config, options, grp, atlas_selector, "tst")


def compute_ind_graph_common_testing(config, options):
    """ Driver to compute a graph representation of an individual using the
    global (i.e., Common) atlas (FOR TESTING DATA!).

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:  # Run on group, e.g., Male, Female
        graph_from_atlas(config, options, grp, "Common", "tst")


def compute_tube_prob(config, options, atlas_type, data_type):
    """ Compute an individuals tube probability, using a group-specific EDM.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    :param atlas_type: Selects the EDM type to use.
    :param data_type: Selects the type of data to use, i.e., training (trn) or testing (trn).
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    atlas_file = os.path.join(cv_dir, atlas_type, "bNormalDD-mean.mha")
    if not check_file(atlas_file):
        raise Exception("Atlas file %s missing!" % atlas_file)

    for i in options[data_type]:

        logger.debug("Running probability computation for subject %s graphs"
            %  os.path.basename(options["subjects"][i]))

        prob_file = os.path.join(cv_dir, os.path.basename(options["subjects"][i]), "on-%s-atlas.treProb.txt" % atlas_type)
        tube_file = os.path.join(options["dest"], os.path.basename(options["subjects"][i]), "VascularNetwork-t1-phantom.tre")

        cmd = [config["Exec"]["ComputeTubeProbability"],
            tube_file, atlas_file, prob_file]
        print cmd

        subprocess.call(cmd)


def compute_ind_tube_prob_testing(config, options):
    """ Driver to compute the tube probability of individuals, using different
    atlases, i.e., group-wise and Common (FOR TESTING DATA!).

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    compute_tube_prob(config, options, "Common", "tst")

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:
        compute_tube_prob(config, options, grp, "tst")


def compute_graph_probability(config, options, grp, atlas_type, data_type):
    """ Compute a pseudo probability of a graph, using the summary graph of
    a group.

    Takes a 1) group-specific summary graph, specified by grp, 2) an atlas type,
    specified by atlas_type (the atlas that was used to generate the individual
    graphs) and computes a pseudo-probability of an individual's graph under that
    summary graph.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    :param grp: Selects the group for which we like to compute individual graph probabilities.
    :param atlas_type: Selects the atlas type that was used to create the individual graphs.
    :param data_type: Selects the type of data to use, i.e., training (trn) or testing (trn).
    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    # This is the basename for the summary graph sub-filetypes .cnt, .mat, .brc
    summary_graph_file = os.path.join(cv_dir, grp, "Normal-mean-%s.grp" % atlas_type)

    for i in options[data_type]:
        in_graph_file = os.path.join(
            cv_dir, os.path.basename(options["subjects"][i]),
            "VascularNetwork-t1-phantom-%s.grp" % atlas_type)

        if not check_file("%s.mat" % in_graph_file):
            raise Exception("Sub-file %s.mat missing!" % in_graph_file)
        if not check_file("%s.rot" % in_graph_file):
            raise Exception("Sub-file %s.rot missing!" % in_graph_file)
        if not check_file("%s.brc" % in_graph_file):
            raise Exception("Sub-file %s.brc missing!" % in_graph_file)

        out_graph_file = os.path.join(cv_dir,
            os.path.basename(options["subjects"][i]),
            "VascularNetwork-t1-phantom-%s-%s.grp" % (atlas_type, grp))

        cmd = [config["Exec"]["ComputeTubeGraphProbability"],
            in_graph_file, summary_graph_file, out_graph_file]
        print cmd
        subprocess.call(cmd)


def compute_ind_graph_prob_testing(config, options):
    """ Driver to compute graph probabilities using different atlases (FOR TESTING DATA!).

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    for grp in group_names:
        compute_graph_probability(config, options, grp, grp, "tst")
        compute_graph_probability(config, options, grp, "Common", "tst")


def create_label_map(options):
    """ Computes a map of textual labels to numeric labels.

    :param options: Configuration options for controlling the experiment flow.
    :returns: A dictionary with textual labels as keys and numeric labels as values.
    """

    group_names = set()
    for grp in options["groupLabel"]:
        group_names.add(grp)

    label_map = dict()
    for label, grp in enumerate(group_names):
        label_map[grp] = label + 1  # Let labels start @ 1
    return label_map


def compose_gk_list(config, options, atlas_type, data_type, graph_list_file):
    """ Composes a list of graphs for graph kernel computation and write the
    list to a JSON file.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    :param atlas_type: Selects which graphs to use, based on the atlas that was used
                       to generate them.
    :param data_type: Selects whether we use training (trn) or testing (tst) data.
    :param graph_list_file: Name of the output JSON file.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    label_map = create_label_map(options)

    graph_list = []
    label_list = []
    for i in options[data_type]:
        label = label_map[options["groupLabel"][i]]
        subject_graph_file = os.path.join(
            cv_dir,
            os.path.basename(options["subjects"][i]),
            "VascularNetwork-t1-phantom-%s.grp.mat" % atlas_type)
        if not check_file(subject_graph_file):
            raise Exception("Subject graph file %s missing!" % subject_graph_file)
        graph_list.append(subject_graph_file)
        label_list.append(label)

    json_dict = { "nGraphs" : len(graph_list), "graphList" : graph_list, "labels" : label_list}
    with open(graph_list_file, "w") as outfile:
        json.dump(json_dict, outfile)


def compute_trn_gk(config, options):
    """Driver for computing the training data graph kernel.

    Creates a list of training graphs from the cross-validation split
    and compute a kernel (Gram) matrix using a graph kernel.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    trn_common_kern_file = os.path.join(cv_dir, "trn-Common.kern")
    trn_common_json_file = os.path.join(cv_dir, "trn-Common.json")
    compose_gk_list(config, options, "Common", "trn", trn_common_json_file)

    cmd = [config["Exec"]["TubeGraphKernel"],
        trn_common_json_file,
        trn_common_json_file,
        trn_common_kern_file,
        "--defaultLabelType %d" % options["defaultLabelType"],
        "--graphKernelType %d" % options["graphKernelType"],
        "--subtreeHeight %d" % options["subtreeHeight"]]

    if not options["globalLabelFile"] is None:
        full_path_to_label_file = os.path.join(cv_dir, "Common", options["globalLabelFile"])
        cmd.append("--globalLabelFile %s" % full_path_to_label_file)

    print cmd
    subprocess.call(cmd)


def compute_tst_gk(config, options):
    """ Driver for computing the testing data graph kernel.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    trn_common_json_file = os.path.join(cv_dir, "trn-Common.json")
    if not check_file(trn_common_json_file):
        raise Exception("JSON file %s for training data missing!" % trn_common_json_file)

    tst_common_json_file = os.path.join(cv_dir, "tst-Common.json")
    tst_common_kern_file = os.path.join(cv_dir, "tst-Common.kern")
    compose_gk_list(config, options, "Common", "tst", tst_common_json_file)

    cmd = [config["Exec"]["TubeGraphKernel"],
        tst_common_json_file,
        trn_common_json_file,
        tst_common_kern_file,
        "--defaultLabelType %d" % options["defaultLabelType"],
        "--graphKernelType %d" % options["graphKernelType"],
        "--subtreeHeight %d" % options["subtreeHeight"]]

    if not options["globalLabelFile"] is None:
        full_path_to_label_file = os.path.join(cv_dir, "Common", options["globalLabelFile"])
        cmd.append("--globalLabelFile %s" % full_path_to_label_file)

    print cmd
    subprocess.call(cmd)


def compute_full_gk(config, options):
    """ Compute graph kernel from all (i.e., training + testing) samples.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    trn_common_json_file = os.path.join(cv_dir, "trn-Common.json")
    tst_common_json_file = os.path.join(cv_dir, "tst-Common.json")

    compose_gk_list(config, options, "Common", "trn", trn_common_json_file)
    compose_gk_list(config, options, "Common", "tst", tst_common_json_file)

    if not check_file(trn_common_json_file) or not check_file(tst_common_json_file):
        raise Exception("JSON files for training/testing missing!")

    # Read dictionaries
    trn_fid = open(trn_common_json_file).read()
    tst_fid = open(tst_common_json_file).read()

    trn_data = json.loads(trn_fid)
    tst_data = json.loads(tst_fid)

    full_data = dict()
    full_data['labels'] = trn_data['labels'] + tst_data['labels']
    full_data['graphList'] = trn_data['graphList'] + tst_data['graphList']
    full_data['nGraphs'] = trn_data['nGraphs'] + tst_data['nGraphs']

    # Write dictionary
    graph_list_file = os.path.join(cv_dir, "full-Common.json" )
    with open(graph_list_file, "w") as outfile:
        json.dump(full_data, outfile)

    full_common_kern_file = os.path.join(cv_dir, "full-Common.kern")
    cmd = [config["Exec"]["TubeGraphKernel"],
        graph_list_file,
        graph_list_file,
        full_common_kern_file,
        "--defaultLabelType %d" % options["defaultLabelType"],
        "--graphKernelType %d" % options["graphKernelType"],
        "--subtreeHeight %d" % options["subtreeHeight"]]
    subprocess.call(cmd)


def normalize_kernel(kernel_data):
    """ Normalize kernel matrix.

    Use K_ij / sqrt( K_ii * K_jj) for normalization

    :param kernel_data: A NxN matrix of kernel values.
    :returns: A normalized kernel matrix.
    """

    diag_mat = np.matrix(kernel_data.diagonal())
    norm_mat = 1./(diag_mat.T * diag_mat)
    norm_mat[np.isinf(norm_mat)] = 0.0

    return np.array(np.multiply(kernel_data, np.sqrt(norm_mat)))


def crossvalidate_cost_param(clf, kernel, labels, n_cv_runs, cost_values, trn_fraction):
    """ Cross-validate cost factor C of SVM using precomputed kernel.

    :param clf: An existing SVM classifier.
    :param kernel: A NxN matrix of training kernel values.
    :param labels: A list of N training labels (one for each row in the kernel).
    :param n_cv_runs: Number of cross-validation runs.
    :param cost_values: An array of cost values to use for cross-validation.
    :param trn_fraction: Fraction of the training data to use for cross-validation.
    :returns: The index into cost_values that gave the best cross-validation result on
              the training portion.
    """

    logger = logging.getLogger()

    n_samples = kernel.shape[0] # Assume NxN kernel matrix
    logger.info("Crossvalidate cost on %d samples" % n_samples)

    n_trn_samples = int(trn_fraction * n_samples)
    logger.info("Training portion = %d samples" % n_trn_samples)

    full_set = set(xrange(n_samples))
    map_scores = list()

    for C in cost_values:
        score_list = list()
        for cv_run in range(1, n_cv_runs):
            trn_set = set(random.sample(xrange(n_samples), n_trn_samples))
            tst_lst = list(full_set - trn_set)
            trn_lst = list(trn_set)

            # Make numpy matrix
            kernel_mat = np.matrix(kernel)

            # Extract training kernel entries (i.e., all train rows, all train cols)
            trn_kernel_mat = np.delete(kernel_mat, tst_lst, axis=0)
            trn_kernel_mat = np.delete(trn_kernel_mat, tst_lst, axis=1)
            trn_lab = np.delete(labels, tst_lst)

            # Extract testing kernel entries (i.e., all test rows, all training cols)
            tst_kernel_mat = np.delete(kernel_mat, trn_lst, axis=0)
            tst_kernel_mat = np.delete(tst_kernel_mat, tst_lst, axis=1)
            tst_lab = np.delete(labels, trn_lst)

            # Next, train classifier with current C value
            clf = svm.SVC(kernel="precomputed", C=1, verbose=False)
            clf.fit(trn_kernel_mat, trn_lab)
            label_hat = clf.predict(tst_kernel_mat)
            score_hat = clf.score(tst_kernel_mat, tst_lab)
            score_list.append(score_hat)

        print "mAP (%.2f) at C=%.2f" % (np.mean(score_list)*100,C)
        map_scores.append(np.mean(score_list)*100)

    return map_scores.index(max(map_scores))


def compute_distance_signatures(config, options):
    """ Computes SGD (summary-graph-distance) features.

    1) For all individuals (trn/tst), compile a list file for SG
       distance computation (i.e., a list of graphs)
    2) Run the kernel on all distance files
    3) Extract relevant kernel entries and compute the induced
       distance, i.e., d(x,y) = || phi(x) - phi(y) ||, where
       phi() is the feature-mapping.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)


    std_groupings = set()
    all_groupings = set()
    for group in options["groupLabel"]:
        std_groupings.add(group)
        all_groupings.add(group)
    all_groupings.add("Common")


    logger.info("Creating SGD list!")

    label_list = []
    graph_list = []

    # Add the summary graphs to the list
    for grp_atlas in std_groupings:
        for atlas in all_groupings:
            # Only with the 'Common' atlas, we stay in the same space
            if not (atlas == "Common"):
                continue
            summary_graph_file = os.path.join(cv_dir, grp_atlas,
                "Normal-mean-%s.grp.mat" % atlas)

            if not check_file(summary_graph_file):
                raise Exception("Summary graph file %s missing!" % summary_graph_file)

            graph_list.append(summary_graph_file)

    # Create summary-graph list JSON file
    json_dict = {"nGraphs" : len(graph_list), "graphList" : graph_list, "labels" : [-1]*len(graph_list)}
    sgd_file = os.path.join(cv_dir, "sgd-list.txt")
    with open(sgd_file, "w") as outfile:
        json.dump(json_dict, outfile)


    # Compute the k(y,y) entries ... (used in kernel-induced distance)
    k2_file = os.path.join(cv_dir, "k2-sgd.kern")
    cmd = [config["Exec"]["TubeGraphKernel"],
            sgd_file,
            sgd_file,
            k2_file,
            "--defaultLabelType %d" % options["defaultLabelType"],
            "--graphKernelType %d" % options["graphKernelType"],
            "--subtreeHeight %d" % options["subtreeHeight"]]
    subprocess.call(cmd)


    logger.info("Creating subject-specific graph lists for SGD")

    label_map = create_label_map(options)
    for cnt,subject in enumerate(options["subjects"]):
        label_list = []
        graph_list = []
        for atlas_used_for_graph in all_groupings:

            if not (atlas_used_for_graph == "Common"):
                continue

            # Add the subject's graph file (for one atlas type)
            subject_graph_file = os.path.join(cv_dir,os.path.basename(subject),
                "VascularNetwork-t1-phantom-%s.grp.mat" % atlas_used_for_graph)

            if not check_file(subject_graph_file):
                raise Exception("File %s missing!" % subject_graph_file)

            graph_list.append(subject_graph_file)
            label_list.append(label_map[options["groupLabel"][cnt]])

        # Dictionary entries ...
        json_dict = {
            "nGraphs"   : len(graph_list),
            "graphList" : graph_list,
            "labels"    : label_list }

        # Write subject graph list to subject's directory
        subject_graph_list_file = os.path.join(cv_dir, os.path.basename(subject),
            "sgd-compare-graph-list.txt")
        with open(subject_graph_list_file, "w") as outfile:
            json.dump(json_dict, outfile)

        # We can now compute the kernel matrices
        k0_file = os.path.join(cv_dir, os.path.basename(subject), "k0-sgd.kern")
        k1_file = os.path.join(cv_dir, os.path.basename(subject), "k1-sgd.kern")

        # First, the k(x,x) entries ... (used in kernel-induced distance)
        cmd = [config["Exec"]["TubeGraphKernel"],
            subject_graph_list_file,
            subject_graph_list_file,
            k0_file,
            "--defaultLabelType %d" % options["defaultLabelType"],
            "--graphKernelType %d" % options["graphKernelType"],
            "--subtreeHeight %d" % options["subtreeHeight"]]
        subprocess.call(cmd)

        # Second, the k(x,y) entries ... (used in kernel-induced distance)
        cmd = [config["Exec"]["TubeGraphKernel"],
            subject_graph_list_file,
            sgd_file,
            k1_file,
            "--defaultLabelType %d" % options["defaultLabelType"],
            "--graphKernelType %d" % options["graphKernelType"],
            "--subtreeHeight %d" % options["subtreeHeight"]]

        subprocess.call(cmd)


def compute_glo_label_map(config, options):
    """ Computes a global label map from CVT cell ID to a discrete label.

    The created mapping files can later be used as input to the graph
    kernel(s) to provide better (more discriminative) labelings of nodes.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    atlas_cvt_file = os.path.join(cv_dir, "Common", "bNormalDD-mean.cvt.mha")
    if not check_file(atlas_cvt_file):
        raise Exception("CVT file %s for group Common not found!" % atlas_cvt_file)

    # Newly created CVT files, one (w)ith background cell counting, the other has (n)o background cell counting
    relabeled_wBack_atlas_cvt_img_file = os.path.join(cv_dir, "Common", "relabeled-wBack-bNormalDD-mean.cvt.mha")
    relabeled_nBack_atlas_cvt_img_file = os.path.join(cv_dir, "Common", "relabeled-nBack-bNormalDD-mean.cvt.mha")

    # CVT cell ID to new label mapping, one (w)ith background cell counting, the other has (n)o background cell counting
    relabeled_wBack_atlas_cvt_map_file = os.path.join(cv_dir, "Common", "relabeled-wBack-bNormalDD-mean.cvt.map")
    relabeled_nBack_atlas_cvt_map_file = os.path.join(cv_dir, "Common", "relabeled-nBack-bNormalDD-mean.cvt.map")

    if options["segmentationImage"] is None:
        raise Exception("No segmentation file given!")

    # Run with background cell labeling
    cmd = [config["Exec"]["TransferLabelsToRegions"],
            atlas_cvt_file,
            options["segmentationImage"],
            relabeled_wBack_atlas_cvt_img_file,
            relabeled_wBack_atlas_cvt_map_file,
            "--omitRegions %d" % 0]
    subprocess.call(cmd)

    # Rum without background cell labeling
    cmd = [config["Exec"]["TransferLabelsToRegions"],
            atlas_cvt_file,
            options["segmentationImage"],
            relabeled_nBack_atlas_cvt_img_file,
            relabeled_nBack_atlas_cvt_map_file]
    subprocess.call(cmd)


def evaluate_classifier_from_full_gk(config, options):
    """ Train and evaluate a SVM classifier, assuming that we have access to all
    data, i.e., training and testing.

    Train a SVM using a full (training + testing) kernel. During training,
    only those kernel values that correspond to training sample pairs are
    considered, obviously.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    N = len(options["trn"]) # M training samples
    M = len(options["tst"]) # N testing samples

    # Create mapping of class to label
    label_map = create_label_map(options)

    # This is the file created by TubeGraphKernel
    kernel_data_file = os.path.join(cv_dir, "full-Common.kern.bin")
    if not check_file(kernel_data_file):
        raise Exception("Training kernel data file %s missing!" % kernel_data_file)

    # Reads kernel from binary file and normalizes it
    kernel_mat = np.fromfile(kernel_data_file, dtype="double").reshape(N+M,N+M);

    # Create training/testing/all labels
    tst_labels = []
    trn_labels = []
    all_labels = []
    for i in options["trn"]:
        label = label_map[options["groupLabel"][i]]
        trn_labels.append(label)
    for i in options["tst"]:
        label = label_map[options["groupLabel"][i]]
        tst_labels.append(label)
    all_labels = trn_labels + tst_labels


    # The set of cost factors, we use for CV
    cost_factors = [0.1, 1, 5, 10, 20, 50];

    # Normalize the kernel data
    normalized_kernel_mat = normalize_kernel(kernel_mat)

    # Extract the training portion of the kernel
    list_of_tst_samples = range(N, N+M)

    assert len(list_of_tst_samples) == M, 'Mismatch in list size of testing samples!'

    trn_kernel_mat = np.matrix(normalized_kernel_mat)
    pruned_trn_kernel_mat = np.delete(trn_kernel_mat, list_of_tst_samples, axis=0)
    pruned_trn_kernel_mat = np.delete(pruned_trn_kernel_mat, list_of_tst_samples, axis=1)

    # Extract the testing portion of the kernel
    list_of_trn_samples = range(0,N)
    assert len(list_of_trn_samples) == N, 'Mismatch in list size of training samples!'

    tst_kernel_mat = np.matrix(normalized_kernel_mat)
    pruned_tst_kernel_mat = np.delete(tst_kernel_mat, list_of_trn_samples, axis=0)
    pruned_tst_kernel_mat = np.delete(pruned_tst_kernel_mat, list_of_tst_samples, axis=1)

    # Create a SVM classifier and cross-validate the cost factor
    clf = svm.SVC(kernel="precomputed", C=1, verbose=True)
    optimal_cost_idx = crossvalidate_cost_param(
        clf, pruned_trn_kernel_mat, trn_labels, 20, cost_factors, 0.5)

    logger.info("Optimal cost = %.2f", cost_factors[optimal_cost_idx])

    clf = svm.SVC(kernel="precomputed", C=cost_factors[optimal_cost_idx], verbose=True)
    score = clf.fit(pruned_trn_kernel_mat, trn_labels).score(pruned_trn_kernel_mat, trn_labels)
    logger.info("Training finished (Score = %.2f)!" % score)

    y = clf.predict(pruned_tst_kernel_mat)
    score = clf.score(pruned_tst_kernel_mat, tst_labels)
    logger.info("Final score %.2f", 100 * score)


def trn_classifier(config, options):
    """ Train a support vector machine classifer.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    N = len(options["trn"])
    label_map = create_label_map(options)

    logger.debug("%d subjects in training data, %d label types!"
        % (N, len(label_map)))

    kernel_data_file = os.path.join(cv_dir, "trn-Common.kern.bin")
    if not check_file(kernel_data_file):
        raise Exception("Training kernel data file %s missing!" % kernel_data_file)

    # Reads kernel from binary file and normalizes it
    kernel_data = np.fromfile(kernel_data_file, dtype="double").reshape(N,N);

    labels = []
    for i in options["trn"]:
        label = label_map[options["groupLabel"][i]]
        labels.append(label)

    logger.debug("Training C-SVM with precomputed kernel ...")
    clf = svm.SVC(kernel="precomputed", C=200, verbose=True)

    score = clf.fit(kernel_data, labels).score(kernel_data, labels)
    logger.debug("Training finished (Score = %.2f)!" % score)
    print float(sum(clf.n_support_))/N

    clf_out_file = os.path.join(cv_dir, "svm-Common.clf")
    with open(clf_out_file, 'wb') as fid:
        cPickle.dump(clf, fid)
    logger.debug("Wrote classifer to %s!" % clf_out_file)

    # Write SVM information file
    svm_info_file = os.path.join(cv_dir, "svm-Info.json")
    svm_info = dict()
    svm_info["num_sv"] = float(sum(clf.n_support_))/N
    svm_info["trn_score"] = score
    with open(svm_info_file, 'wb') as fp:
        json.dump(svm_info, fp)
    logger.debug("Wrote classifier info to %s!" % svm_info_file)


def tst_classifier(config, options):
    """ Test the support vector machine classifier.

    :param config: Configuration options for the executables.
    :param options: Configuration options for controlling the experiment flow.
    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    N = len(options["trn"])
    M = len(options["tst"])
    label_map = create_label_map(options)

    logger.debug("Testing on %d samples (%d used in training)" % (M,N))

    labels = []
    for i in options["tst"]:
        label = label_map[options["groupLabel"][i]]
        labels.append(label)

    kernel_data_file = os.path.join(cv_dir, "tst-Common.kern.bin")
    if not check_file(kernel_data_file):
        raise Exception("Testing kernel data file %s missing!" % kernel_data_file)
    kernel_data = np.fromfile(kernel_data_file, dtype="double").reshape(M,N);

    # Load classifier from disk
    clf_file = os.path.join(cv_dir, "svm-Common.clf")
    if not check_file(clf_file):
        raise Exception("Classifier %s missing!" % clf_file)
    with open(clf_file, "rb") as fid:
        clf = cPickle.load(fid)

    y = clf.predict(kernel_data)

    score = clf.score(kernel_data, labels)
    logger.debug("Testing score = %.2f (%d/%d)" % (score, int(score*M), M))

    # Try to update the SVM information file
    svm_info = dict()
    svm_info_file = os.path.join(cv_dir, "svm-Info.json")
    if check_file(svm_info_file):
        fp = open(svm_info_file)
        svm_info = json.load(fp)
        svm_info["tst_score"] = score
        svm_info["truth"] = list(labels)
        svm_info["guess"] = list(y)
        fp.close()

        with open(svm_info_file, 'wb') as fp:
            json.dump(svm_info, fp)
