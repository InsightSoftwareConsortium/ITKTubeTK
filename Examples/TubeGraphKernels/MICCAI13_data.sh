#!/bin/bash

#------------------------------------------------------------------------------
# Data generation pipeline
#
# Author: Roland Kwitt, Kitware Inc., 2013
# E-Mail: roland.kwitt@kitware.com, rkwitt@gmx.at
#------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# READ THIS BEFORE YOU BEGIN
#------------------------------------------------------------------------------
#
#    It is highly recommended that you run stages 1-3 (registration, mapping of
#    vessels into phantom space and density image computation) before running
#    this script, since these steps need to be executed only once. Consequently
#    runtime can be considerably reduced.
#
#    THIS SCRIPT IS CONFIGURED TO RUN FROM STAGE 4!
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# ADJUST THIS
#-------------------------------------------------------------------------------
CVT_CELLS=( 512 )
GLOBAL_LABEL_FILES=( "relabeled-wBack-bNormalDD-mean.cvt.map"
                     "relabeled-nBack-bNormalDD-mean.cvt.map" )
SUBTREE_HEIGHTS=( 2 4 6 )
EXPERIMENTS_BASE_DIR="/tmp/Output/"
#EXPERIMENTS_BASE_DIR="/Users/rkwitt/Data/MICCAI2013"
DRIVER_EXEC="/opt/local/bin/python2.7 expdrive.py"
#DATA_CONFIG_FILE="/Users/rkwitt/Data/Data-Config.json"
DATA_CONFIG_FILE="/Volumes/Skunk/Data/Data-Config.json"
#GENERAL_CONFIG_FILE="/Users/rkwitt/Data/General-Config.json"
GENERAL_CONFIG_FILE="/Volumes/Skunk/Data/General-Config.json"
#PHANTOM="/Users/rkwitt/Data/MICCAI2013/spl-phantom-2012-SkullStripped.mha"
PHANTOM="/Volumes/Skunk/Data/MICCAI2013/spl-phantom-2012-SkullStripped.mha"
#SEGMENTATION_IMAGE="/Users/rkwitt/Data/MICCAI2013/spl-phantom-2012-labels.nrrd"
SEGMENTATION_IMAGE="/Volumes/Skunk/Data/MICCAI2013/spl-phantom-2012-labels.nrrd"
PHANTOM_TYPE="SPL"
N_CV_RUNS=5
#------------------------------------------------------------------------------

for cells in ${CVT_CELLS[@]}; do
    printf -v dir_num "%.4d" $cells
    out_dir_name="${EXPERIMENTS_BASE_DIR}/CVT${dir_num}"

    case $1 in
        #-----------------------------------------------------------------------
        # Stage 1: ATLAS construction
        #-----------------------------------------------------------------------
        #
        # Run throuh the various ATLAS construction stages; This will create
        # N_CV_RUNS directories cv-xxxx under the directory that you specified
        # in EXPERIMENTS_BASE_DIR/CVT{$CVT_CELLS}, e.g.,
        #
        # ...
        # /tmp/data/CVT1024/cv-0000
        # /tmp/data/CVT1024/cv-0001
        # ...
        #-----------------------------------------------------------------------
        1)
            for stage in $(seq 4 15); do
                cmd="${DRIVER_EXEC} \
                    --data ${DATA_CONFIG_FILE} \
                    --config ${GENERAL_CONFIG_FILE} \
                    --phantom ${PHANTOM} \
                    --cvruns ${N_CV_RUNS} \
                    --phantomType ${PHANTOM_TYPE} \
                    --stage ${stage} \
                    --cells ${cells} \
                    --dest ${EXPERIMENTS_BASE_DIR}"
                echo $cmd
                $cmd
            done
            mkdir -p ${out_dir_name} # Create new directory

            # Move cv-xxxx results to the newly created directory
            cmd="mv ${EXPERIMENTS_BASE_DIR}/cv* ${out_dir_name}"
            $cmd
            ;;

        #-----------------------------------------------------------------------
        # Kernel computation
        #-----------------------------------------------------------------------
        # We first move the cross-validation directories of stage 1 back
        # to the base EXPERIMENTS_BASE_DIR directory and then run a WL kernel
        # in different configurations on the data to create the training and
        # testing Gram matrices for each CV run.
        #
        # The kernel results are then moved into subdirectories that
        # identify the kernel height + labeling strategy, e.g.,
        #
        # ...
        # /tmp/data/CVT1024/height_02_label_01/cv-0000
        # /tmp/data/CVT1024/height_02_label_00/cv-0000
        # /tmp/data/CVT1024/relabeled-wBack-bNormalDD-mean-height_02/
        # ...
        #
        # For the WL-subtree kernel we vary:
        #   - Subtree height
        #   - Node labeling strategy (degree, ID, global label files)
        #-----------------------------------------------------------------------
        2)
            # Temorarily move cv-xxxx back to base experiments directory
            cmd="mv ${EXPERIMENTS_BASE_DIR}/CVT${dir_num}/cv-* \
                 ${EXPERIMENTS_BASE_DIR}"
            $cmd

            # Then, compute kernels on graphs with DEFAULT node labelings
            for height in ${SUBTREE_HEIGHTS[@]}; do

                # Label defaults for nodes, (0,1) = (degree, ID)
                for node_label in 0 1; do

                    for stage in $(seq 17 18); do
                        cmd="${DRIVER_EXEC} \
                            --data ${DATA_CONFIG_FILE} \
                            --config ${GENERAL_CONFIG_FILE} \
                            --phantom ${PHANTOM} \
                            --cvruns ${N_CV_RUNS} \
                            --phantomType ${PHANTOM_TYPE} \
                            --stage ${stage} \
                            --cells ${cells} \
                            --dest ${EXPERIMENTS_BASE_DIR} \
                            --graphKernelType 1 \
                            --defaultLabelType ${node_label} \
                            --subtreeHeight ${height}"
                        echo $cmd
                        $cmd
                    done # END for stage ...

                    # Store the kernel results in ATLAS + CV specific folders
                    printf -v subdir_id \
                        "height_%.2d_label_%.2d" ${height} ${node_label}
                    subdir=${out_dir_name}/${subdir_id}

                    for cv in $(seq 1 ${N_CV_RUNS}); do
                        printf -v cv_dir "%.4d" ${cv}
                        mkd_cmd="mkdir -p ${subdir}/cv-${cv_dir}"
                        mv0_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-${cv_dir}/trn* \
                            ${subdir}/cv-${cv_dir}/"
                        mv1_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-${cv_dir}/tst* \
                            ${subdir}/cv-${cv_dir}/"

                        $mkd_cmd
                        $mv0_cmd
                        $mv1_cmd
                    done # END for cv ...

                done # END for node_label ...

                # KERNELS on graphs with PRECOMPUTED (GLOBAL) node labelings
                for label_file in ${GLOBAL_LABEL_FILES[@]}; do

                    global_id=`basename $label_file .cvt.map`
                    printf -v suffix "height_%.2d" ${height}
                    subdir="${out_dir_name}/${global_id}_${suffix}"

                    for stage in $(seq 16 18); do
                        cmd="${DRIVER_EXEC} \
                            --data ${DATA_CONFIG_FILE} \
                            --config ${GENERAL_CONFIG_FILE} \
                            --phantom ${PHANTOM} \
                            --cvruns ${N_CV_RUNS} \
                            --phantomType ${PHANTOM_TYPE} \
                            --stage ${stage} \
                            --cells ${cells} \
                            --dest ${EXPERIMENTS_BASE_DIR} \
                            --graphKernelType 1 \
                            --subtreeHeight ${height} \
                            --segmentationImage ${SEGMENTATION_IMAGE} \
                            --globalLabelFile ${label_file}"
                        $cmd
                    done # END for stage ...

                    for cv in $(seq 1 ${N_CV_RUNS}); do
                        printf -v cv_dir "%.4d" ${cv}
                        mkd_cmd="mkdir -p ${subdir}/cv-${cv_dir}"
                        mv0_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-${cv_dir}/trn* \
                            ${subdir}/cv-${cv_dir}/"
                        mv1_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-${cv_dir}/tst* \
                            ${subdir}/cv-${cv_dir}/"
                        $mkd_cmd
                        $mv0_cmd
                        $mv1_cmd
                    done # END for cv ...

                done # END for label_file

            done # END for height

            # Move the cv-xxxx data back to the CVT specific directory
            cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-* ${out_dir_name}"
            $cmd
            ;;

        #-----------------------------------------------------------------------
        # Run SVM training and testing based on the results of stage 2
        #-----------------------------------------------------------------------
        #
        # All computations are done in the base experiments directory, i.e.,
        # EXPERIMENTS_BASE_DIR.
        #
        # So, first we move the training/testing kernel data for each kernel
        # type/kernel height/labeling to the base experiments directory; then
        # we run SVM training and testing and eventually move the results back
        # into the corresponding directories.
        #-----------------------------------------------------------------------
        3)
            for height in ${SUBTREE_HEIGHTS[@]}; do
                for node_label in 0 1; do

                    printf -v subdir_id "height_%.2d_label_%.2d" ${height} ${node_label}
                    mv_cmd="mv ${EXPERIMENTS_BASE_DIR}/CVT${dir_num}/${subdir_id}/cv-* \
                            ${EXPERIMENTS_BASE_DIR}"
                    echo $mv_cmd
                    ${mv_cmd}

                    # Call driver to execute training and testing
                    for stage in $(seq 20 21); do
                        cmd="${DRIVER_EXEC} \
                            --data ${DATA_CONFIG_FILE} \
                            --config ${GENERAL_CONFIG_FILE} \
                            --cvruns ${N_CV_RUNS} \
                            --phantomType ${PHANTOM_TYPE} \
                            --stage ${stage} \
                            --phantom ${PHANTOM} \
                            --dest ${EXPERIMENTS_BASE_DIR}"
                        echo $cmd
                        $cmd
                    done # END for stage ...

                    # Move everthing back to the param-specific directory
                    mv_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-* \
                            ${EXPERIMENTS_BASE_DIR}/CVT${dir_num}/${subdir_id}"
                    echo ${mv_cmd}
                    ${mv_cmd}

                done # END for node_label

                # now the global labels
                for label_file in ${GLOBAL_LABEL_FILES[@]}; do

                    global_id=`basename $label_file .cvt.map`
                    printf -v suffix "height_%.2d" ${height}
                    subdir="${global_id}_${suffix}"

                    mv_cmd="mv ${EXPERIMENTS_BASE_DIR}/CVT${dir_num}/${subdir}/cv-* \
                        ${EXPERIMENTS_BASE_DIR}"
                    $mv_cmd

                    for stage in $(seq 20 21); do
                        cmd="${DRIVER_EXEC} \
                            --data ${DATA_CONFIG_FILE} \
                            --config ${GENERAL_CONFIG_FILE} \
                            --cvruns ${N_CV_RUNS} \
                            --phantomType ${PHANTOM_TYPE} \
                            --stage ${stage} \
                            --phantom ${PHANTOM} \
                            --dest ${EXPERIMENTS_BASE_DIR}"
                        echo $cmd
                        $cmd
                    done # END for stage ...

                    mv_cmd="mv ${EXPERIMENTS_BASE_DIR}/cv-* \
                        ${EXPERIMENTS_BASE_DIR}/CVT${dir_num}/${subdir}"
                    echo $mv_cmd
                    $mv_cmd
                done # END for label_file

            done # END for height
            ;;
    esac

done # END for cells ...
