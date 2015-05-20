#!/bin/bash

if [ ! "$#" -eq 4 ]; then
  echo "UA_SegmentVessels.sh <trainDir> <seedThreshold> <inputImage> <outputVesselsBaseName>"
  exit
fi

if [ ! -e $3 ]; then
  echo "Cannot find the image file = $3"
  exit
fi

if [ ! -e ${TubeTKExperimentsPath}/Projects/UltrasoundAngiography/VesselExtraction/SegmentTubesParameters/$1/vess.mtp ]; then
  echo "Cannot find the training directory ${TubeTKExperimentsPath}/Projects/UltrasoundAngiography/VesselExtraction/SegmentTubesParameters/$1"
  exit
fi

rm -rf $4_tubetk
mkdir $4_tubetk

cd $4_tubetk

echo "1/6: Normalizing input image"
ImageMath ../$3 -i 0 255 0 1 -w inputFull.mha -r 2 -w input.mha  > trainLog.txt

echo "2/6: Enhancing vessels: using scales 0.2,0.4"
EnhanceTubesUsingDiscriminantAnalysis --loadDiscriminantInfo ${TubeTKExperimentsPath}/Projects/UltrasoundAngiography/VesselExtraction/SegmentTubesParameters/$1/vess.mrs --outputSeedScaleImage vessEnh-scales.mha --useIntensityOnly input.mha vessEnh.mha >> trainLog.txt

echo "3/6: Computing seeds"
ImageMath vessEnh.mha -t 0 $2 0 1 -W 2 vessEnh-vessMask$2.mha >> trainLog.txt
SegmentBinaryImageSkeleton vessEnh-vessMask$2.mha vessEnh-seedMask$2.mha >> trainLog.txt

echo "4/6: Resampling images"
ResampleImage --matchImage ../$3 --interpolator NearestNeighbor vessEnh-seedMask${2}.mha vessEnh-seedMask${2}Full_fat.mha >> trainLog.txt
SegmentBinaryImageSkeleton vessEnh-seedMask${2}Full_fat.mha vessEnh-seedMask${2}Full.mha >> trainLog.txt
ResampleImage --matchImage ../$3 vessEnh-scales.mha vessEnh-scalesFull.mha >> trainLog.txt

echo "5/6: Segmenting tubes"
SegmentTubes -o ../$4.so.tre.mha -P ${TubeTKExperimentsPath}/Projects/UltrasoundAngiography/VesselExtraction/SegmentTubesParameters/$1/vess.mtp -S vessEnh-scalesFull.mha -M vessEnh-seedMask${2}Full.mha -T 10 inputFull.mha ../$4.so.tre >> trainLog.txt

cd ..

echo "6/6: Converting results to $4"
ConvertTRE -r $4.so.tre $4.tx.tre >> $4_tubetk/trainLog.txt

echo "Done."
