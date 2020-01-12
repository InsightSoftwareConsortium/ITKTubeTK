#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "TubeTKTubeEnhancementApplier <inputImage> <trainingBasename> <outputImage>"
  echo "  <inputImage> = input image (must already be normalized 0 to 1)"
  echo "  <trainingBasename> = basename for trained parameters"
  echo "  <outputImage> = enhanced image"
  echo ""
  exit
fi

rm $3.trainLog.txt

echo "1/1: Enhancing vessels"
EnhanceTubesUsingDiscriminantAnalysis --loadDiscriminantInfo $2.mrs --outputSeedScaleImage $3.scales.mha --useIntensityOnly $1 $3 >> $3.log.txt
