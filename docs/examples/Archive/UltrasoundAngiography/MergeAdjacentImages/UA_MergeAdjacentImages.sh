#!/bin/bash

# Merge the two given adjacent images into the output image
#
if [ "$#" -eq 0 ]; then
  echo "Usage:"
  echo "  UA_MergeAdjacentImages.sh <leftImage> <rightImage> <outputImage>"
  echo "or"
  echo "  UA_MergeAdjacentImages.sh <leftImage> <rightImage> <tfmFile> <outputImage>"
  exit
fi

if [ "$#" -eq 3 ]; then
  echo "Initial merge estimate..."
  mkdir ${3}_Temp
  ImageMath ${1} -b 0.2 -r 3 -W 2 ${3}_Temp/left.mha > ${3}.log.txt
  ImageMath ${2} -b 0.2 -r 3 -W 2 ${3}_Temp/right.mha >> ${3}.log.txt
  MergeAdjacentImages -a -L ${TubeTKExperimentsPath}/Projects/UltrasoundAngiography/MergeAdjacentImages/UA_MergeInitialTransform.tfm -S ${3}_Temp/finalTransform.tfm -r 0.01 -M -s 0.03 ${3}_Temp/left.mha ${3}_Temp/right.mha ${3}_Temp/temp.mha >> ${3}.log.txt
  echo "Precise merge..."
  MergeAdjacentImages -L ${3}_Temp/finalTransform.tfm -S ${3}.tfm -s 0.01 -i 50 -M ${1} ${2} ${3} >> ${3}.log.txt
  exit
else
  echo "Merging using pre-computed transform..."
  MergeAdjacentImages -L ${3} -o 0 -r 0 -i 1 -M ${1} ${2} ${4} >> ${4}.log.txt
  exit
fi
