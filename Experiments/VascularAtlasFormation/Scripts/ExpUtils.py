import os
import sys
import json
import glob
import cPickle
import logging
import subprocess
import numpy as np
from sklearn import svm


def check_file(file):
    """ Checks if file exists """
    if file is None:
        return False
    return os.path.exists(file) and os.path.isfile(file)


def ensure_dir(d):
    """ Checks if directory exists; if 'False', it creates it """
    if not os.path.exists(d):
        print "Creating directory %s" % d
        os.makedirs(d)


def check_dir(d):
    """ Checks existence of a directory """
    return os.path.exists(d)


def compute_registrations(config, options):
    """

    Runs image registrations to compute transform from MRA to Phantom.
    No need for cross-validation here, since the registrations do not
    depend on the training/testing splits.

    """

    logger = logging.getLogger();

    for dir in options["subjects"]:
        mra_wSkull_listing = glob.glob(os.path.join(dir, options["mra_wSkull_glob"]))
        mri_wSkull_listing = glob.glob(os.path.join(dir, options["mri_wSkull_glob"]))
        mri_nSkull_listing = glob.glob(os.path.join(dir, options["mri_nSkull_glob"]))

        sum_len = 0
        sum_len += len( mri_wSkull_listing )
        sum_len += len( mra_wSkull_listing )
        sum_len += len( mri_nSkull_listing )

        if not sum_len == 3:
            raise Exception( "Mismatch in nr. of MRI/MRA files!" )

        ensure_dir(os.path.join(options["dest"], os.path.basename(dir)))

        rig_transform_file = os.path.join(options["dest"], os.path.basename(dir), "MRA-T1.tfm" )
        aff_transform_file = os.path.join(options["dest"], os.path.basename(dir), "T1-Phantom.tfm")
        rig_image_file = os.path.join(options["dest"], os.path.basename(dir), "mra-t1.mha")
        aff_image_file = os.path.join(options["dest"], os.path.basename(dir), "t1-phantom.mha")

        param_type = "--registration Rigid"
        param_save = "--saveTransform %s" % rig_transform_file
        param_fimg = "%s" % mra_wSkull_listing[0]
        param_mimg = "%s" % mri_wSkull_listing[0]
        param_samp =  "--resampledImage %s" % rig_image_file

        # MRA to T1 registration
        logger.debug("Running rigid registration to compute %s -> %s transform",
            os.path.basename(param_fimg), os.path.basename(param_mimg))

        cmd_reg0 = [ config["Exec"]["RegisterImages"],
            param_type, param_save, param_fimg, param_mimg, param_samp]
        print cmd_reg0
        subprocess.call(cmd_reg0)

        param_type = "--registration Affine"
        param_save = "--saveTransform %s" % aff_transform_file
        param_fimg = "%s" % mri_nSkull_listing[0]
        param_mimg = "%s" % options["phantom"]
        param_samp =  "--resampledImage %s" % aff_image_file

        # T1 to Phantom registration
        logger.debug("Running rigid registration to compute %s -> %s transform",
            os.path.basename(param_fimg), os.path.basename(param_mimg))

        cmd_reg1 = [config["Exec"]["RegisterImages"],
            param_type, param_save, param_fimg, param_mimg, param_samp]
        subprocess.call(cmd_reg1)


def transform_tubes_to_phantom(config, options):
    """

    Transform tubes (i.e., vascular networks) into Phantom space, using
    the two image registration transforms, i.e., from MRA to T1 and from
    T1 to Phantom.

    """

    logger = logging.getLogger()

    for dir in options["subjects"]:
        tubes_original_file = os.path.join(dir, "color-new.tre" )

        if (not check_file(tubes_original_file)):
            raise Exception("Original tube file %s not existent!" % tubes_original_file)

        # Output tube files (to be generated in subjects destination directory)
        tubes_mra_t1_file = os.path.join(options["dest"], os.path.basename(dir), "color-new-mra-t1.tre" )
        tubes_t1_pha_file = os.path.join(options["dest"], os.path.basename(dir), "color-new-t1-phantom.tre")

        # Transform files (assumed to exist in destination directory of subject)
        mra_t1_tfm_file = os.path.join(options["dest"], os.path.basename(dir), "MRA-T1.tfm" )
        t1_pha_tfm_file = os.path.join(options["dest"], os.path.basename(dir), "T1-Phantom.tfm" )

        if (not check_file(mra_t1_tfm_file)):
            raise Exception("MRA to T1 transform file %s missing!" % mra_t1_tfm_file)
        if (not check_file(t1_pha_tfm_file)):
            raise Exception("T1 to Phantom transform file %s missing!" % t1_pha_tfm_file)

        # Tubes (in MRA image space) to MRI T1 image space
        logger.debug("MRA to T1 space for %s (using %s transform)"
            % (tubes_original_file, mra_t1_tfm_file))

        t0_cmd = [config["Exec"]["tubeTransform"], tubes_original_file, mra_t1_tfm_file, tubes_mra_t1_file]
        subprocess.call(t0_cmd)

        # Tubes (in MRI T1 image space) to Phantom image space
        logger.debug("T1 to Phantom space for %s (using %s transform)"
            % (tubes_mra_t1_file, tubes_t1_pha_file))

        t1_cmd = [config["Exec"]["tubeTransform"], tubes_mra_t1_file, t1_pha_tfm_file, tubes_t1_pha_file]
        subprocess.call(t1_cmd)


def compute_ind_atlas_edm(config, options):
    """

    Computes (ind)ividual (e)euclidean (d)istance (m)ap -- EDM -- (using Daniellson's algorithm)

    """
    logger = logging.getLogger()

    for dir in options["subjects"]:
        den_image_name = os.path.join(options["dest"], os.path.basename(dir), "color-new-t1-phantomDD.mha")
        rad_image_name = os.path.join(options["dest"], os.path.basename(dir), "color-new-t1-phantomRad.mha")
        tan_image_name = os.path.join(options["dest"], os.path.basename(dir), "color-new-t1-phantomTan.mha")

        # Tubes in phantom space (assumed to exist at that stage)
        tubes_in_phantom_space = os.path.join(options["dest"], os.path.basename(dir), "color-new-t1-phantom.tre" )
        if (not check_file(tubes_in_phantom_space)):
            raise Exception( "Tube file %s missing!" % tubes_in_phantom_space )

        logger.debug("Create individual EDM for %s"
            % os.path.basename(dir))

        cmd = [config["Exec"]["tubeDensityImageRadiusBuilder"],
            tubes_in_phantom_space, den_image_name, rad_image_name, tan_image_name,
            "--inputTemplateImage %s" % options["phantom"],
            "--useSquareDistance"]
        subprocess.call(cmd)


def compute_grp_atlas_sum(config, options):
    """

    Compute group-wise EDMs (i.e., the ATLAS EDM) from individual EDMs

    """
    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    ensure_dir(cv_dir) # Create if necessary

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for grp in group_names:
        grp_dir = os.path.join(cv_dir, grp) # e.g., <cvdir>/Female
        ensure_dir(grp_dir) # Create if necessary

        pop_cnt = 0
        for i in options["trn"]:
            if not options["grouplabel"][i] == grp:
                continue
            pop_cnt += 1

        atlas_edm_doc = os.path.join(grp_dir, "NormalDD.odc")
        fid_atlas_edm_doc = open(atlas_edm_doc, "w")
        fid_atlas_edm_doc.write('NumberOfObjects = %d\n' % pop_cnt)

        for i in options["trn"]:
            if not options["grouplabel"][i] == grp:
                continue

            ind_edm_file = os.path.join(options["dest"], os.path.basename(options["subjects"][i]),"color-new-t1-phantomDD.mha")
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

        cmd = [config["Exec"]["tubeAtlasSummation"],
            atlas_edm_doc, atlas_edm_file, atlas_edm_var_file,
            "--outputSize %d,%d,%d" % (181, 217, 181),
            "--outputSpacing %d,%d,%d" % (1,1,1)]
        subprocess.call(cmd)


def compute_grp_atlas_cvt(config, options):
    """

    Computes a Central Voronoi Tesselation (CVT) from the Euclidean
    distance map per group. This will represent the group-specific
    atlas.

    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["grouplabel"]:
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
            "-S", "1234", # RNG seed for reproducable results
            "-Z", "%d" % options["cells"], "1000", "10000", atlas_cvt_map_file,
            "-W", "1", atlas_cvt_img_file]
        subprocess.call(cmd)


def graph_from_atlas(config, options, grp, atlas_type, data_type):
    """ Graph from ATLAS Computation

    For each subject in our dataset, the routine takes 1) the tubes (in phantom)
    space and 2) a group atlas and computes a graph representation of the tubes.

    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation base directory %s missing!" % cv_dir)

    atlas_file = os.path.join(cv_dir, atlas_type, "bNormalDD-mean.cvt.mha")
    if not check_file(atlas_file):
        raise Exception("Atlas %s missing!" % atlas_file)

    for i in options[data_type]:
        if not options["grouplabel"][i] == grp:
            continue

        subject_dir = os.path.join(options["dest"], os.path.basename(options["subjects"][i]))
        tubes_in_phantom_space = os.path.join(subject_dir, "color-new-t1-phantom.tre" )

        if (not check_file(tubes_in_phantom_space)):
            raise Exception("Transformed tubes (%s) are missing!" % tubes_in_phantom_space)

        # Directory that contains subject-specific graphs (depend on CVT) for each cross-validation
        cv_subject_dir = os.path.join(cv_dir, os.path.basename(options["subjects"][i]))
        ensure_dir(cv_subject_dir)

        logger.debug("Compute graph for %s individual %s (using ATLAS type %s)"
            % (grp, os.path.basename(options["subjects"][i]), atlas_type))

        output_graph_file = os.path.join(cv_subject_dir, "color-new-t1-phantom-%s.grp" % atlas_type )
        cmd = [config["Exec"]["tubeToGraph"],
            tubes_in_phantom_space, atlas_file, output_graph_file]
        subprocess.call(cmd)


def compute_ind_graph_grp(config, options):
    """

    Driver for computing one graph for each subject based on a given
    group atlas.

    E.g., for subject A and Male, Female groups, we compute

        - one graph for subject A using the Male atlas
        - one graph for subject A using the Female atlas

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_atlas(config, options, grp, atlas_selector, "trn")


def graph_from_graphs(config, options, grp, atlas_type, create_image=False):
    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    pop_cnt = 0
    for i in options["trn"]:
        if not options["grouplabel"][i] == grp:
            continue
        pop_cnt += 1

    graph_doc_file = os.path.join(cv_dir, grp, "NormalGraph-%s.odc" % atlas_type)
    fid = open( graph_doc_file, "w" )
    fid.write( "NumberOfObjects = %d\n" % pop_cnt)

    for i in options["trn"]:
        if not options["grouplabel"][i] == grp:
            continue

        # Assume existence of tubes (in Phantom space) for each subject
        cv_subject_dir = os.path.join(cv_dir, os.path.basename(options["subjects"][i]))
        tube_graph = os.path.join(cv_subject_dir, "color-new-t1-phantom-%s.grp" % atlas_type)
        if (not check_file(tube_graph)):
            raise Exception("Graph %s missing!" % tube_graph)

        fid.write( "Type = Image\n" )
        fid.write( "Name = %s\n" % tube_graph )
        fid.write( "EndObject = \n" )

    fid.close()

    cmd = [config["Exec"]["tubeGraphsToGraph"],
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

        cmd = [config["Exec"]["tubeGraphToImage"],
            atlas_cvt_file,
            os.path.join(cv_dir, grp, "Normal-mean-%s.grp" % atlas_type),
            os.path.join(cv_dir, grp, "Normal-mean-%s.cvt" % atlas_type)]
        subprocess.call(cmd)


def compute_grp_graph_grp(config, options):
    """

    Driver routine for computing one representative (summary) graph per
    group from different atlases.

    E.g., for Female, Male group, compute

        - one Male graph using Female atlas
        - one Male graph using Male atlas
        - one Female graph using Male atlas
        - one Female graph using Female atlas

    """

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_graphs(config, options, grp, atlas_selector, False)


def compute_glo_atlas_edm(config, options):
    """

    Computes a (glo)bal Euclidean distance map for all groups,
    based on the group-specific Euclidean distance maps. Each
    distance map per group is equally weighted.

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["grouplabel"]:
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
    """

    Computes a (glo)bal CVT from the gobal Euclidean distance map.
    This represents the so called 'Common' atlas.

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    glo_dd_cvt_map_file = os.path.join(cv_dir, "Common", "bNormalDD-mean.cvt" )
    glo_dd_cvt_img_file = os.path.join(cv_dir, "Common", "bNormalDD-mean.cvt.mha" )

    common_dir = os.path.join(cv_dir, "Common")
    ensure_dir(common_dir)

    glo_dd_edm_file = os.path.join(common_dir, "bNormalDD-mean.mha")
    cmd = [config["Exec"]["ImageMath"],
        glo_dd_edm_file,
        "-Z", "%d" % options["cells"],
        "1000", "10000",
        glo_dd_cvt_map_file,
        "-W", "1", glo_dd_cvt_img_file]
    subprocess.call(cmd)


def compute_ind_graph_common(config, options):
    """

    Compute one graph for each individual using the 'Common' atlas.

    """

    logger = logging.getLogger()

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    common_dir = os.path.join(cv_dir, "Common")
    if not check_dir(common_dir):
        raise Exception("Directory %s missing!" % common_dir)

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    # Runs graph computation using different ATLASes
    for grp in group_names:
        graph_from_atlas(config, options, grp, "Common", "trn")


def compute_grp_graph_common(config, options):
    """

    Compute one graph for each group from the individual graphs that
    were computed using the 'Common' atlas.

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    common_dir = os.path.join(cv_dir, "Common")
    if not check_dir(common_dir):
        raise Exception("Directory %s missing!" % common_dir)

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for grp in group_names:
        graph_from_graphs(config, options, grp, "Common", create_image=True)


def compute_ind_graph_grp_testing(config, options):
    """

    Computes a graph representation for all individuals in the
    groups, based on the corresponding group ATLAS(es).


    Uses only testing data!

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for atlas_selector in group_names:
        for grp in group_names:
            graph_from_atlas(config, options, grp, atlas_selector, "tst")


def compute_ind_graph_common_testing(config, options):
    """

    Computes a graph representation for all individuals in the
    groups, using a common ATLAS.

    Uses only testing data!

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for grp in group_names:  # Run on group, e.g., Male, Female
        graph_from_atlas(config, options, grp, "Common", "tst")


def compute_tube_prob(config, options, atlas_type, data_type):
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
        tube_file = os.path.join(options["dest"], os.path.basename(options["subjects"][i]), "color-new-t1-phantom.tre")

        cmd = [config["Exec"]["tubeDensityProbability"],
            tube_file, atlas_file, prob_file]
        subprocess.call(cmd)


def compute_ind_tube_prob_testing(config, options):
    """ Compute Pseudo-Probability of a Testing Tube File

    Driver to compute a pseudo-probability for each testing tube file under
    different atlas EDM files, e.g., Female, Male, Common.

    """

    compute_tube_prob(config, options, "Common", "tst")

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for grp in group_names:
        compute_tube_prob(config, options, grp, "tst")


def compute_graph_probability(config, options, grp, atlas_type, data_type):
    """ Compute Pseudo-Probability of a Graph

    Takes a 1) summary group graph, indicated by 'grp', 2) an atlas type, specified
    by 'atlas_type' (the atlas that was used to generate the individual graphs) and
    computes a pseudo-probability of an individual's graph under that summary graph.

    """

    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    # This is the basename for the summary graph sub-filetypes .cnt, .mat, .brc
    summary_graph_file = os.path.join(cv_dir, grp, "Normal-mean-%s.grp" % atlas_type)

    for i in options[data_type]:
        in_graph_file = os.path.join(
            cv_dir, os.path.basename(options["subjects"][i]),
            "color-new-t1-phantom-%s.grp" % atlas_type)

        if not check_file("%s.mat" % in_graph_file):
            raise Exception("Sub-file %s.mat missing!" % in_graph_file)
        if not check_file("%s.rot" % in_graph_file):
            raise Exception("Sub-file %s.rot missing!" % in_graph_file)
        if not check_file("%s.brc" % in_graph_file):
            raise Exception("Sub-file %s.brc missing!" % in_graph_file)

        out_graph_file = os.path.join(cv_dir,
            os.path.basename(options["subjects"][i]),
            "color-new-t1-phantom-%s-%s.grp" % (atlas_type, grp))

        cmd = [config["Exec"]["tubeGraphProbability"],
            in_graph_file, summary_graph_file, out_graph_file]
        print cmd
        subprocess.call(cmd)


def compute_ind_graph_prob_testing(config, options):
    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    for grp in group_names:
        compute_graph_probability(config, options, grp, grp, "tst")
        compute_graph_probability(config, options, grp, "Common", "tst")



def create_label_map(options):
    """ Returns a map of textual labels to numeric values """

    group_names = set()
    for grp in options["grouplabel"]:
        group_names.add(grp)

    label_map = dict()
    for label, grp in enumerate(group_names):
        label_map[grp] = label + 1  # Let labels start @ 1
    return label_map


def compose_gk_list(config, options, atlas_type, data_type, graph_list_file):
    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("CV directory %s missing!" % cv_dir)

    label_map = create_label_map(options)

    graph_list = []
    label_list = []
    for i in options[data_type]:
        label = label_map[options["grouplabel"][i]]
        subject_graph_file = os.path.join(
            cv_dir,
            os.path.basename(options["subjects"][i]),
            "color-new-t1-phantom-%s.grp.mat" % atlas_type)
        if not check_file(subject_graph_file):
            raise Exception("Subject graph file %s missing!" % subject_graph_file)
        graph_list.append(subject_graph_file)
        label_list.append(label)

    json_dict = { "nGraphs" : len(graph_list), "graphList" : graph_list, "labels" : label_list}
    with open(graph_list_file, "w") as outfile:
        json.dump(json_dict, outfile)


def compute_trn_gk(config, options):
    """Driver for training data graph kernel

    Creates a list of training graphs from the cross-validation split
    and compute a kernel matrix using a graph-kernel.

    """
    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    trn_common_kern_file = os.path.join(cv_dir, "trn-Common.kern")
    trn_common_json_file = os.path.join(cv_dir, "trn-Common.json")
    compose_gk_list(config, options, "Common", "trn", trn_common_json_file)

    cmd = [config["Exec"]["tubeGraphKernel"],
        trn_common_json_file, trn_common_json_file, trn_common_kern_file]
    subprocess.call(cmd)


def compute_tst_gk(config, options):
    """ Driver for testing data graph kernel

    see compute_trn_gk for details!

    """
    cv_dir = os.path.join(options["dest"], "cv-%.4d" % options["id"])
    if not check_dir(cv_dir):
        raise Exception("Cross-validation directory %s missing!" % cv_dir)

    trn_common_json_file = os.path.join(cv_dir, "trn-Common.json")
    if not check_file(trn_common_json_file):
        raise Exception("JSON file %s for training data missing!" % trn_common_json_file)

    tst_common_json_file = os.path.join(cv_dir, "tst-Common.json")
    tst_common_kern_file = os.path.join(cv_dir, "tst-Common.kern")
    compose_gk_list(config, options, "Common", "tst", tst_common_json_file)

    cmd = [config["Exec"]["tubeGraphKernel"],
        tst_common_json_file, trn_common_json_file, tst_common_kern_file]
    subprocess.call(cmd)


def trn_classifier(config, options):
    """ Train a support vector machine classifer """

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

    # Reads kernel from binary file
    kernel_data = np.fromfile(kernel_data_file, dtype="double").reshape(N,N);

    labels = []
    for i in options["trn"]:
        label = label_map[options["grouplabel"][i]]
        labels.append(label)

    logger.debug("Training C-SVM with precomputed kernel ...")
    clf = svm.SVC(kernel="precomputed")
    score = clf.fit(kernel_data, labels).score(kernel_data, labels)
    logger.debug("Training finished (Score = %.2f)!" % score)

    clf_out_file = os.path.join(cv_dir, "svm-Common.clf")
    with open(clf_out_file, 'wb') as fid:
        cPickle.dump(clf, fid)
    logger.debug("Wrote classifer to %s!" % clf_out_file)


def tst_classifier(config, options):
    """ Test the support vector machine classifier """

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
        label = label_map[options["grouplabel"][i]]
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
