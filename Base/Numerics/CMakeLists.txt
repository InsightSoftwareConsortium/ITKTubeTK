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

project( TubeNumerics )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_Numerics_H_Files
  itkImageRegionMomentsCalculator.h
  itkJointHistogramImageFunction.h
  itkTubeBlurImageFunction.h
  itkTubeLDAGenerator.h
  itkTubeNJetImageFunction.h
  itkTubeNJetLDAGenerator.h
  itkTubeNJetLDAGenerator2F.h
  itkVectorImageToListGenerator.h
  itkVotingResampleImageFunction.h
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
  tubeUserFunction.h
  )

set( TubeTK_Base_Numerics_TXX_Files
  itkImageRegionMomentsCalculator.txx
  itkJointHistogramImageFunction.txx
  itkTubeLDAGenerator.txx
  itkTubeNJetImageFunction.txx
  itkTubeNJetLDAGenerator.txx
  itkVectorImageToListGenerator.txx
  itkVotingResampleImageFunction.txx
  tubeMatrixMath.txx
  tubeTubeMath.txx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  )

set( TubeTK_Base_Numerics_SRCS
  itkTubeNJetLDAGenerator2F.cxx
  tubeBrentOptimizer1D.cxx
  tubeGoldenMeanOptimizer1D.cxx
  tubeOptimizer1D.cxx
  tubeOptimizerND.cxx
  tubeParabolicFitOptimizer1D.cxx
  tubeSpline1D.cxx
  tubeSplineApproximation1D.cxx
  tubeSplineND.cxx
  )

add_library( TubeNumerics STATIC
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_TXX_Files}
  ${TubeTK_Base_Numerics_SRCS}
  )

target_link_libraries( TubeNumerics ${ITK_LIBRARIES} )

slicer4_set_plugins_output_path( TubeNumerics )
set( TARGETS TubeNumerics )
slicer4_install_plugins( ${TARGETS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_TXX_Files}
  DESTINATION include
  )
