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

project( TubeFiltering )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_Numerics_H_Files
  itkImageRegionMomentsCalculator.h
  itkJointHistogramImageFunction.h
  itkTubeBlurImageFunction.h
  itkTubeLDAGenerator.h
  itkTubeNJetImageFunction.h
  itkTubeNJetLDAGenerator2F.h
  itkTubeNJetLDAGenerator.h
  itkVectorImageToListGenerator.h
  itkVotingResampleImageFunction.h
  tubeMatrixMath.h
  tubeOptBrent1D.h
  tubeOptGoldenMean1D.h
  tubeOptimizer1D.h
  tubeOptimizerND.h
  tubeOptParabolicFit1D.h
  tubeSpline1D.h
  tubeSplineApproximation1D.h
  tubeSplineND.h
  tubeTubeMath.h
  tubeUserFunc.h
  )

set( TubeTK_Base_Numerics_TXX_Files
  itkImageRegionMomentsCalculator.txx
  itkJointHistogramImageFunction.txx
  itkTubeLDAGenerator.txx
  itkTubeNJetImageFunction.txx
  itkTubeNJetLDAGenerator.txx
  itkVectorImageToListGenerator.txx
  itkVotingResampleImageFunction.txx
  tubeTubeMath.txx
  tubeMatrixMath.txx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  )

set( TubeTK_Base_Numerics_SRCS
  itkTubeNJetLDAGenerator2F.cxx
  tubeOptBrent1D.cxx
  tubeOptGoldenMean1D.cxx
  tubeOptParabolicFit1D.cxx
  tubeOptimizer1D.cxx
  tubeOptimizerND.cxx
  tubeSpline1D.cxx
  tubeSplineApproximation1D.cxx
  tubeSplineND.cxx
  )

add_library( TubeNumerics STATIC ${TubeTK_Base_Numerics_SRCS} )
target_link_libraries( TubeNumerics ${ITK_LIBRARIES} )

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
#if( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  #set_target_properties(TubeNumerics PROPERTIES COMPILE_FLAGS "-fPIC")
#endif( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )

if( NOT TubeTK_USE_ITKV4 )
  slicer3_set_plugins_output_path( TubeNumerics )
  set( TARGETS TubeNumerics )
  slicer3_install_plugins( ${TARGETS} )
endif( NOT TubeTK_USE_ITKV4 )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install(FILES
  ${TubeTK_Base_Numerics_H_Files}
  ${TubeTK_Base_Numerics_TXX_Files}
  DESTINATION include
  )
