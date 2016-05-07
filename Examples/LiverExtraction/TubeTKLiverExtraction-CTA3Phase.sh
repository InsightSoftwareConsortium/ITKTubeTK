#!/bin/bash

export boundarySmoothness=3
export imageSmoothing=0.5
export histogramSmoothing=10
export holeFillIterations=10
export inputImageIntensityMin=0
export inputImageIntensityMax=550

if [ "$#" -eq 0 ]; then
  echo "TubeTKLiverExtraction-CTA3Phase <inputImage_Phase1> <inputImage_Phase2> <inputImage_Phase3> <inputLiverTrainMask> <outputLiverMask>"
  echo "  <inputImage_PhaseX> : Liver CT scan with 3 contract uptake phases"
  echo "  <inputLiverTrainMask> : Voxel value 1 = outside of liver, 2 = inside of liver"
  echo "     - this image is usually created by painting a few voxels in the image"
  echo "  <outputLiverMask> : will have 1 inside liver and 0 elsewhere"
  echo "  Note that other parameters are set within this script file:"
  echo "        boundarySmoothness = $boundarySmoothness"
  echo "          - e.g., 12 for 0.8mm isotropic CT"
  echo "          - e.g., 10 for 1mm isotropic CT"
  echo "          - e.g., 5 for 2mm isotropic CT"
  echo "        imageSmoothing = $imageSmoothing"
  echo "          - e.g., 0 for standard energy CT"
  echo "          - e.g., 2 for low-dose CT"
  exit
fi

rm $3.log.txt

mkdir $5_Temp

echo "1/2: Normalizing images"
ImageMath $1 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm1.mha >> $5.log.txt
ImageMath $2 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm2.mha >> $5.log.txt
ImageMath $3 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w $5_Temp/norm3.mha >> $5.log.txt

ImageMath $5_Temp/norm3.mha -a 1 -1 $5_Temp/norm1.mha -w $5_Temp/sub31.mha >> $5.log.txt
ImageMath $5_Temp/norm3.mha -a 1 -1 $5_Temp/norm2.mha -w $5_Temp/sub32.mha >> $5.log.txt
ImageMath $5_Temp/norm2.mha -a 1 -1 $5_Temp/norm1.mha -w $5_Temp/sub21.mha >> $5.log.txt

echo "2/2: Creating probability image"
SegmentConnectedComponentsUsingParzenPDFs --objectId 2,1 $5_Temp/norm3.mha --inputVolume2 $5_Temp/sub21.mha --inputVolume3 $5_Temp/sub31.mha --inputVolume4 $5_Temp/sub32.mha $4 $5_Temp/liverMask.mha --holeFillIterations $holeFillIterations --erodeRadius $boundarySmoothness --probImageSmoothingStdDev $imageSmoothing --histogramSmoothingStdDev $histogramSmoothing

ImageMath $5_Temp/liverMask.mha -t 2 2 1 0 -W 0 $5
