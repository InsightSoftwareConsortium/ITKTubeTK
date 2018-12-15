##############################################################################
#
# Library:   TubeTKLib
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

####
# Track local source files
####
set( TubeTKLib_Base_Numerics_H_Files
  itktubeBasisFeatureVectorGenerator.h
  itktubeBlurImageFunction.h
  itktubeComputeImageSimilarityMetrics.h
  itktubeComputeImageStatistics.h
  itktubeFeatureVectorGenerator.h
  itktubeImageRegionMomentsCalculator.h
  itktubeJointHistogramImageFunction.h
  itktubeNJetFeatureVectorGenerator.h
  itktubeNJetImageFunction.h
  itktubeRecordOptimizationParameterProgressionCommand.h
  itktubeRidgeFFTFeatureVectorGenerator.h
  itktubeVectorImageToListGenerator.h
  itktubeVotingResampleImageFunction.h
  tubeBrentOptimizer1D.h
  tubeGoldenMeanOptimizer1D.h
  tubeMatrixMath.h
  tubeOptimizer1D.h
  tubeOptimizerND.h
  tubeParabolicFitOptimizer1D.h
  tubeSpline1D.h
  tubeSplineApproximation1D.h
  tubeSplineND.h
  tubeTubeMath.h
  tubeUserFunction.h )

set( TubeTKLib_Base_Numerics_HXX_Files
  itktubeBasisFeatureVectorGenerator.hxx
  itktubeBlurImageFunction.hxx
  itktubeComputeImageSimilarityMetrics.hxx
  itktubeComputeImageStatistics.hxx
  itktubeFeatureVectorGenerator.hxx
  itktubeImageRegionMomentsCalculator.hxx
  itktubeJointHistogramImageFunction.hxx
  itktubeNJetFeatureVectorGenerator.hxx
  itktubeNJetImageFunction.hxx
  itktubeRecordOptimizationParameterProgressionCommand.hxx
  itktubeRidgeFFTFeatureVectorGenerator.hxx
  itktubeVectorImageToListGenerator.hxx
  itktubeVotingResampleImageFunction.hxx
  tubeMatrixMath.hxx
  tubeTubeMath.hxx )

set( TubeTKLib_Base_Numerics_SRCS
  tubeBrentOptimizer1D.cxx
  tubeGoldenMeanOptimizer1D.cxx
  tubeOptimizer1D.cxx
  tubeOptimizerND.cxx
  tubeParabolicFitOptimizer1D.cxx
  tubeSpline1D.cxx
  tubeSplineApproximation1D.cxx
  tubeSplineND.cxx )

####
# Create project Library
####
add_library( TubeTKLibNumerics STATIC
  ${TubeTKLib_Base_Numerics_H_Files}
  ${TubeTKLib_Base_Numerics_HXX_Files}
  ${TubeTKLib_Base_Numerics_SRCS} )

if( TubeTK_USE_ARRAYFIRE )
  target_link_libraries( TubeTKLibNumerics PUBLIC
    ${ArrayFire_LIBRARIES} )
endif()

target_link_libraries( TubeTKLibNumerics PUBLIC
  TubeTKLibCommon
  ${ITK_LIBRARIES} )

target_include_directories( TubeTKLibNumerics PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

####
# Do testing
####
if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_Numerics_H_Files}
  ${TubeTKLib_Base_Numerics_HXX_Files}
  DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS TubeTKLibNumerics )
