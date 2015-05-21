#!/bin/bash

export inputImageIntensityMin=0
export inputImageIntensityMax=255

export bufferRadius=2

export scales="1,2,4"

if [ "$#" -eq 0 ]; then
  echo "TubeTKTubeEnhancementTrainer <inputImage> <vesselMask> <trainingResultsBasename>"
  echo "  <inputImage> = image used to learn the parameters of the enhancement filter"
  echo "  <vesselMask> = Voxel value 1 = within tube, 0 = unknown"
  echo "  <trainingResultsBasename> = basename for a set of parameter files generates"
  echo "  Note that other parameters are set within this script file:"
  echo "        inputImageIntensityMin = $inputImageIntensityMin"
  echo "        inputImageIntensityMax = $inputImageIntensityMax"
  echo "        bufferRadius = $bufferRadius"
  echo "        scales = $scales"
  echo ""
  exit
fi

rm $3.trainLog.txt

echo "1/3: Normalizing input image"
echo "   Mapping from $inputImageIntensityMin - $inputImageIntensityMax to 0 - 1"
ImageMath $1 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $3.input.mha > $3.trainLog.txt

echo "2/3: Creating vessel training mask"
echo "   Using bufferRadius = $bufferRadius"
ImageMath $2 -M 1 $bufferRadius 1 0 -W 2 $3.dilate1.mha >> $3.trainLog.txt
ImageMath $3.dilate1.mha -M 1 $bufferRadius 1 0 -W 2 $3.dilate2.mha >> $3.trainLog.txt
ImageMath $3.dilate2.mha -a 1 -1 $3.dilate1.mha -a 127 255 $2 -W 2 $3.tubeMask.mha >> $3.trainLog.txt

echo "3/3: Enhancing vessels"
echo "   Using scales = $scales"
EnhanceTubesUsingDiscriminantAnalysis --saveDiscriminantInfo $3.mrs --outputSeedScaleImage $3.inputEnhScales.mha --tubeScales $scales --labelmap $3.tubeMask.mha --useIntensityOnly $3.input.mha $3.inputEnh.mha >> $3.trainLog.txt
