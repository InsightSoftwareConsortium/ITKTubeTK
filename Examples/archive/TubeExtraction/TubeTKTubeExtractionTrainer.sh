#!/bin/bash

export bufferRadius=2

if [ "$#" -eq 0 ]; then
  echo "TubeTKTubeExtractionTrainer.sh <inputImage> <vesselEnhancedImage> <vesselEnhancedScalesImage> <vesselThreshold> <trainingResultsBasename>"
  echo "  <inputImage> = image used to learn the parameters"
  echo "  <vesselEnhancedImage> = output image from TubeTKTubeEnhance.sh script"
  echo "  <vesselEnhancedScalesImage> = output scales image from TubeTKTubeEnhance.sh script"
  echo "  <vesselThreshold> = value to threshold the <vesselEnhancedImage> to define seeds"
  echo "  <trainingResultsBasename> = basename for files created by training"
  echo "  Note that other parameters are set within this script file:"
  echo "        bufferRadius = $bufferRadius"
  echo "            region around vessel use to define background"
  echo ""
  exit
fi

rm -rf $5.mtp.log.txt

mkdir $5.mtp_Temp

echo "1/3: Create seed mask"
ImageMath $2 -t 0 $4 0 1 -W 2 $5.mtp_Temp/vessMask.mha >> $5.mtp.log.txt

SegmentBinaryImageSkeleton $5.mtp_Temp/vessMask.mha $5.mtp_Temp/vessSeedsMask.mha >> $5.mtp.log.txt

echo "2/3: Creating vessel training mask"
ImageMath $5.mtp_Temp/vessMask.mha -M 1 $bufferRadius 1 0 -W 2 $5.mtp_Temp/dilate1.mha -M 1 $bufferRadius 1 0 -a 1 -1 $5.mtp_Temp/dilate1.mha -a 127 255 $5.mtp_Temp/vessSeedsMask.mha -W 2 $5.mtp_Temp/vessSeedsBkgMask.mha >> $5.mtp.log.txt

echo "3/3: Computing parameters"
ComputeSegmentTubesParameters $5.mtp_Temp/vessSeedsBkgMask.mha $3 $1 $5.mtp >> $5.mtp.log.txt
