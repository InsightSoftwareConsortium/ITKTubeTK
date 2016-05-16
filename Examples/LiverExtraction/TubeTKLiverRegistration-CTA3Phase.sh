#!/bin/bash

if [ "$#" -eq 0 ]; then
  echo "TubeTKLiverRegistration-CTA3Phase <inputImage_Phase1> <inputImage_Phase2> <inputImage_Phase3>"
  echo "  <inputImage_PhaseX> : Liver CT scan with 3 contrast uptake phases"
  exit
fi

rm $3.log.txt

RegisterImages $1 $2 --registration PipelineAffine --sampleFromOverlap --rigidSamplingRatio 0.05 --rigidMaxIterations 300 --affineSamplingRatio 0.05 --affineMaxIterations 300 --resampledImage $2.reg.mha --saveTransform $2.reg.tfm --initialization None

RegisterImages $1 $3 --registration PipelineAffine --sampleFromOverlap --rigidSamplingRatio 0.05 --rigidMaxIterations 300 --affineSamplingRatio 0.05 --affineMaxIterations 300 --resampledImage $3.reg.mha --saveTransform $3.reg.tfm --initialization None
