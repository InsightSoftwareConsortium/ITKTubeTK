#!/bin/bash

export inputImageIntensityMin=0
export inputImageIntensityMax=550

export bufferRadius=2

export scales="0.5,1,2,4"

if [ "$#" -eq 0 ]; then
  echo "TubeTKTubeEnhancementTrainer <inputImage1> <inputImage2> <inputImage3> <vesselMask> <trainingResultsBasename>"
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

mkdir $5_Temp
rm $5.log.txt

echo "1/3: Normalizing input image"
echo "   Mapping from $inputImageIntensityMin - $inputImageIntensityMax to 0 - 1"
ImageMath $1 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm1.mha >> $5.log.txt
ImageMath $2 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm2.mha >> $5.log.txt
ImageMath $3 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm3.mha >> $5.log.txt

ImageMath $5_Temp/norm3.mha -a 1 -1 $5_Temp/norm1.mha -w $5_Temp/sub31.mha >> $5.log.txt
ImageMath $5_Temp/norm3.mha -a 1 -1 $5_Temp/norm2.mha -w $5_Temp/sub32.mha >> $5.log.txt
ImageMath $5_Temp/norm2.mha -a 1 -1 $5_Temp/norm1.mha -w $5_Temp/sub21.mha >> $5.log.txt

echo "2/3: Creating vessel training mask"
echo "   Using bufferRadius = $bufferRadius"
ImageMath $4 -M 1 $bufferRadius 1 0 -W 2 $5_Temp/vess.dilate1.mha >> $5.log.txt
ImageMath $5_Temp/vess.dilate1.mha -M 1 $bufferRadius 1 0 -W 2 $5_Temp/vess.dilate2.mha >> $5.log.txt
ImageMath $5_Temp/vess.dilate2.mha -a 1 -1 $5_Temp/vess.dilate1.mha -a 127 255 $4 -W 2 $5_Temp/vess.tubeMask.mha >> $5.log.txt

echo "3/3: Enhancing vessels"
echo "   Using scales = $scales"
EnhanceTubesUsingDiscriminantAnalysis --saveDiscriminantInfo $5.mrs --outputSeedScaleImage $5.enhScales.mha --tubeScales $scales --labelmap $5_Temp/vess.tubeMask.mha $5_Temp/sub31.mha,$5_Temp/sub32.mha,$5_Temp/sub21.mha,$5_Temp/norm3.mha $5.enh.mha --useIntensityOnly >> $5.log.txt
