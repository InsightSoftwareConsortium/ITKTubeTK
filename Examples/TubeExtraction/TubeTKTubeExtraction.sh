#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "TubeTKTubeExtracter.sh <inputImage> <inputImageScaleImage> <seedImage> <trained.mtp> <outputBasename>"
  echo "  <inputImage> = input image"
  echo "  <inputImageScaleImage> = input image's corresponding ScaleImage from TubeEnhancer"
  echo "  <trained.mtp> = files created by training using TubeTKTubeExtractionTrainer"
  echo "  <outputBasename> = basename for saving extracted tubes"
  echo ""
  exit
fi

rm -rf $5
mkdir $5

rm $5/log.txt

echo "1/1: Segmenting tubes"
SegmentTubes -P $4 -M $3 -S $2 -o $5.outVess.mha $1 $5.outVess.tre >> $5/log.txt
