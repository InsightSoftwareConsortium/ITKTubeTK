##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
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

project( TubeTKNumerics )

set( TubeTK_Base_Numerics_H_Files
  itktubeBasisFeatureVectorGenerator.h
  itktubeBlurImageFunction.h
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

set( TubeTK_Base_Numerics_HXX_Files
  itktubeBasisFeatureVectorGenerator.hxx
  itktubeBlurImageFunction.hxx
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

set( TubeTK_Base_Numerics_SRCS
  tubeBrentOptimizer1D.cxx
  tubeGoldenMeanOptimizer1D.cxx
  tubeOptimizer1D.cxx
  tubeOptimizerND.cxx
  tubeParabolicFitOptimizer1D.cxx
  tubeSpline1D.cxx
  tubeSplineApproximation1D.cxx
  tubeSplineND.cxx )

add_library( TubeTKNumerics STATIC
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_HXX_Files}
  ${TubeTK_Base_Numerics_SRCS} )
if( TubeTK_USE_GPU_ARRAYFIRE )
  target_link_libraries( TubeTKNumerics ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_GPU_ARRAYFIRE )
target_link_libraries( TubeTKNumerics PUBLIC TubeTKCommon )
target_include_directories( TubeTKNumerics PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_HXX_Files}
  DESTINATION include )

TubeTKMacroInstallPlugins( TubeTKNumerics )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKNumerics )
