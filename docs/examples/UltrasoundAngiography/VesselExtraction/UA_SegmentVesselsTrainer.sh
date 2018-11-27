#!/bin/bash

if [ ! "$#" -eq 3 ]; then
  echo "TubeTKSegmentVesselsTrainer <seedThreshold> <image> <vessels.tubex.tre>"
  exit
fi

if [ ! -e $2 ]; then
  echo "Cannot find the image file = " $2
  exit
fi

if [ ! -e $3 ]; then
  echo "Cannot find the vessel file = " $3
  exit
fi

#mkdir $3_$1_train
#
#echo "01/14: Normalizing input image"
#ImageMath $2 -i 0 255 0 1 -w $3_$1_train/inputFull.mha -r 2 -w $3_$1_train/input.mha >> $3_$1_train/trainLog.txt
#
#echo "02/14: Converting TubeX tubes to TubeTK tubes"
#ConvertTRE $3 $3_$1_train/tubetk.so.tre > $3_$1_train/trainLog.txt
#
echo "03/14: Entering subdirectory: $3_$1_train"
cd $3_$1_train
#
#echo "04/14: Rendering tubes into mask"
#ConvertTubesToImage -r ../$2 tubetk.so.tre tubetk.so.tre.mha >> trainLog.txt
#
#echo "05/14: Creating vessel training mask"
#ImageMath tubetk.so.tre.mha -M 0 1 1 0 -W 2 erode1.mha >> trainLog.txt
#ImageMath tubetk.so.tre.mha -M 1 1 1 0 -W 2 dilate1.mha >> trainLog.txt
#ImageMath tubetk.so.tre.mha -M 1 4 1 0 -W 2 dilate4.mha >> trainLog.txt
#ImageMath dilate4.mha -a 1 -1 dilate1.mha -a 127 255 erode1.mha -W 2 vessMaskFull.mha >> trainLog.txt
#
#echo "06/14: Creating half-sized images for speed"
#ImageMath vessMaskFull.mha -r 2 -w vessMask.mha >> trainLog.txt
#
#echo "07/14: Enhancing vessels: using scales 0.2,0.4"
#EnhanceTubesUsingDiscriminantAnalysis --saveDiscriminantInfo vess.mrs --tubeScales 0.2,0.4 --labelmap vessMask.mha --useIntensityOnly --outputSeedScaleImage vessEnh-scales.mha input.mha vessEnh.mha >> trainLog.txt
#
#echo "08/14: Computing seeds"
#ImageMath vessEnh.mha -t 0 ${1} 0 1 -W 2 vessEnh-vessMask$1.mha >> trainLog.txt
#SegmentBinaryImageSkeleton vessEnh-vessMask$1.mha vessEnh-seedMask${1}.mha >> trainLog.txt
#
#echo "09/14: Computing SegmentTubes parameters"
#ComputeSegmentTubesParameters -v 1 vessEnh-seedMask$1.mha vessEnh-scales.mha input.mha vess.mtp >> trainLog.txt

#echo "11/14: Computing full-size images"
#ResampleImage --matchImage ../$2 --interpolator NearestNeighbor vessEnh-seedMask${1}.mha vessEnh-seedMask${1}Full_fat.mha >> trainLog.txt
#SegmentBinaryImageSkeleton vessEnh-seedMask${1}Full_fat.mha vessEnh-seedMask${1}Full.mha >> trainLog.txt
#ResampleImage --matchImage ../$2 vessEnh-scales.mha vessEnh-scalesFull.mha >> trainLog.txt

echo "12/14: Segmenting tubes"
#SegmentTubes -o ../$3_$1_train.so.tre.mha -P vess.mtp -S vessEnh-scalesFull.mha -M vessEnh-seedMask${1}Full.mha -T 10 inputFull.mha ../$3_$1_train.so.tre >> trainLog.txt
SegmentTubes -o ../$3_$1_train.so.tre.mha -P vess.mtp -S vessEnh-scales.mha -M vessEnh-seedMask${1}.mha -T 10 input.mha ../$3_$1_train.so.tre >> trainLog.txt

echo "13/14: Converting results from .so.tre format to TubeX.tre format"
ConvertTRE -r ../$3_$1_train.so.tre ../$3_$1_train.tx.tre >> trainLog.txt

echo "14/14: Returning to data directory .."
cd ..

echo "Done."
