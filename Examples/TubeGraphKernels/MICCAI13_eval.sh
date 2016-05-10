#!/bin/bash

#-------------------------------------------------------------------------------
# This is the driver script to compute confusion matrices and P-VALUES for
# ACCURACY and the F1-SCORE of a binary (gender) classifier! The number of
# permutations is set the 1000.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# ADJUST THIS
#-------------------------------------------------------------------------------
# Adjust the Python path, the BASE and the LIST variable to reflect your
# experiment configuration!
#-------------------------------------------------------------------------------

PYTHON="/opt/local/bin/python2.7"
BASE="/Volumes/Skunk/Data/MICCAI2013"
LIST=("CVT0256/height_02_label_00"
      "CVT0256/height_04_label_00"
      "CVT0256/height_06_label_00"
      "CVT0256/height_02_label_01"
      "CVT0256/height_04_label_01"
      "CVT0256/height_06_label_01"

      "CVT0512/height_02_label_00"
      "CVT0512/height_04_label_00"
      "CVT0512/height_06_label_00"
      "CVT0512/height_02_label_01"
      "CVT0512/height_04_label_01"
      "CVT0512/height_06_label_01"

      "CVT1024/height_02_label_00"
      "CVT1024/height_04_label_00"
      "CVT1024/height_06_label_00"
      "CVT1024/height_02_label_01"
      "CVT1024/height_04_label_01"
      "CVT1024/height_06_label_01"

      "CVT2048/height_02_label_00"
      "CVT2048/height_04_label_00"
      "CVT2048/height_06_label_00"
      "CVT2048/height_02_label_01"
      "CVT2048/height_04_label_01"
      "CVT2048/height_06_label_01"

      "CVT0256/relabeled-wBack-bNormalDD-mean_height_02"
      "CVT0256/relabeled-wBack-bNormalDD-mean_height_04"
      "CVT0256/relabeled-wBack-bNormalDD-mean_height_06"
      "CVT0256/relabeled-nBack-bNormalDD-mean_height_02"
      "CVT0256/relabeled-nBack-bNormalDD-mean_height_04"
      "CVT0256/relabeled-nBack-bNormalDD-mean_height_06"

      "CVT0512/relabeled-wBack-bNormalDD-mean_height_02"
      "CVT0512/relabeled-wBack-bNormalDD-mean_height_04"
      "CVT0512/relabeled-wBack-bNormalDD-mean_height_06"
      "CVT0512/relabeled-nBack-bNormalDD-mean_height_02"
      "CVT0512/relabeled-nBack-bNormalDD-mean_height_04"
      "CVT0512/relabeled-nBack-bNormalDD-mean_height_06"

      "CVT1024/relabeled-wBack-bNormalDD-mean_height_02"
      "CVT1024/relabeled-wBack-bNormalDD-mean_height_04"
      "CVT1024/relabeled-wBack-bNormalDD-mean_height_06"
      "CVT1024/relabeled-nBack-bNormalDD-mean_height_02"
      "CVT1024/relabeled-nBack-bNormalDD-mean_height_04"
      "CVT1024/relabeled-nBack-bNormalDD-mean_height_06"

      "CVT2048/relabeled-wBack-bNormalDD-mean_height_02"
      "CVT2048/relabeled-wBack-bNormalDD-mean_height_04"
      "CVT2048/relabeled-wBack-bNormalDD-mean_height_06"
      "CVT2048/relabeled-nBack-bNormalDD-mean_height_02"
      "CVT2048/relabeled-nBack-bNormalDD-mean_height_04"
      "CVT2048/relabeled-nBack-bNormalDD-mean_height_06")

      #-------------------------------------------------------------------------------
# DO NOT MODIFY THIS!
#-------------------------------------------------------------------------------

for experiment in ${LIST[@]}; do
    cmd="/opt/local/bin/python2.7 permtest.py \
        --dir ${BASE}/${experiment}
        --ncv 5 \
        --npt 1000 \
        --trnb trn-Common\
        --tstb tst-Common"
    echo ${experiment}
    $cmd
done
