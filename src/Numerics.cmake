##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################

set( TubeTK_Numerics_H_Files
  Numerics/itktubeBasisFeatureVectorGenerator.h
  Numerics/itktubeBlurImageFunction.h
  Numerics/itktubeComputeImageSimilarityMetrics.h
  Numerics/itktubeComputeImageStatistics.h
  Numerics/itktubeFeatureVectorGenerator.h
  Numerics/itktubeImageRegionMomentsCalculator.h
  Numerics/itktubeJointHistogramImageFunction.h
  Numerics/itktubeNJetFeatureVectorGenerator.h
  Numerics/itktubeNJetImageFunction.h
  Numerics/itktubeRecordOptimizationParameterProgressionCommand.h
  Numerics/itktubeRidgeFFTFeatureVectorGenerator.h
  Numerics/itktubeVectorImageToListGenerator.h
  Numerics/itktubeVotingResampleImageFunction.h
  Numerics/tubeBrentOptimizer1D.h
  Numerics/tubeGoldenMeanOptimizer1D.h
  Numerics/tubeMatrixMath.h
  Numerics/tubeOptimizer1D.h
  Numerics/tubeOptimizerND.h
  Numerics/tubeParabolicFitOptimizer1D.h
  Numerics/tubeSpline1D.h
  Numerics/tubeSplineApproximation1D.h
  Numerics/tubeSplineND.h
  Numerics/tubeUserFunction.h )

set( TubeTK_Numerics_HXX_Files
  Numerics/itktubeBasisFeatureVectorGenerator.hxx
  Numerics/itktubeBlurImageFunction.hxx
  Numerics/itktubeComputeImageSimilarityMetrics.hxx
  Numerics/itktubeComputeImageStatistics.hxx
  Numerics/itktubeFeatureVectorGenerator.hxx
  Numerics/itktubeImageRegionMomentsCalculator.hxx
  Numerics/itktubeJointHistogramImageFunction.hxx
  Numerics/itktubeNJetFeatureVectorGenerator.hxx
  Numerics/itktubeNJetImageFunction.hxx
  Numerics/itktubeRecordOptimizationParameterProgressionCommand.hxx
  Numerics/itktubeRidgeFFTFeatureVectorGenerator.hxx
  Numerics/itktubeVectorImageToListGenerator.hxx
  Numerics/itktubeVotingResampleImageFunction.hxx
  Numerics/tubeMatrixMath.hxx )

set( TubeTK_Numerics_CXX_Files
  Numerics/tubeBrentOptimizer1D.cxx
  Numerics/tubeGoldenMeanOptimizer1D.cxx
  Numerics/tubeOptimizer1D.cxx
  Numerics/tubeOptimizerND.cxx
  Numerics/tubeParabolicFitOptimizer1D.cxx
  Numerics/tubeSpline1D.cxx
  Numerics/tubeSplineApproximation1D.cxx
  Numerics/tubeSplineND.cxx )

list( APPEND TubeTK_SRCS
  ${TubeTK_Numerics_H_Files}
  ${TubeTK_Numerics_HXX_Files}
  ${TubeTK_Numerics_CXX_Files} )
