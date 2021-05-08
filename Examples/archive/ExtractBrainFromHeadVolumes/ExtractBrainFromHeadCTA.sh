#!/bin/bash

export boundarySmoothness=10
export imageSmoothing=0

if [ "$#" -eq 0 ]; then
  echo "TubeTKBrainExtractionFromCT <inputImage> <inputSeedMask> <outputBrainMask>"
  echo "  <inputImage> : image used to learn the parameters of the enhancement filter"
  echo "  <inputSeedMask> : Voxel value 1 = outside of brain, 2 = inside of brain"
  echo "     - this image is usually created by painting a few voxels in the image"
  echo "  <outputBrainMask> : will have 1 insight brain and 0 elsewhere"
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

rm $3.trainLog.txt

echo "1/2: Extracting brain"
SegmentConnectedComponentsUsingParzenPDFs --objectId 2,1 $1 $2 $3 --holeFillIterations 10 --erodeRadius $boundarySmoothness --probSmoothingStdDev $imageSmoothing

echo "2/2: Extracting brain"
ImageMath $1 -i $inputImageIntensityMin $inputImageIntensityMax 0 1 -w input.mha > $3.trainLog.txt
