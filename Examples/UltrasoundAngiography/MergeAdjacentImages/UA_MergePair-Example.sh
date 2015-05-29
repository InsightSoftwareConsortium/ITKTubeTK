#!/bin/bash

time ./UA_MergeAdjacentImages.sh ${1}_bmode_left_2*.mhd ${1}_bmode_right_2*mhd ${1}_bmode_merged.mha
time ./UA_MergeAdjacentImages.sh ${1}_df_left_2*.mhd ${1}_df_right_2*mhd ${1}_bmode_merged.mha.tfm ${1}_df_merged.mha
