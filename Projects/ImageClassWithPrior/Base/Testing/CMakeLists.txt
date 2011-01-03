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

include_regular_expression( "^.*$" )

find_package( Tclsh )
if( TCL_TCLSH )
  add_test( PrintSelf-tubeProjImageClassBase ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl
    ${TubeTK_SOURCE_DIR}/Projects/ImageClassWithPrior/Base )
endif( TCL_TCLSH )

set( BASE_CLASSIFICATION_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeProjImageClassBaseTests
  )

set( BASE_CLASSIFICATION_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeProjImageClassBaseHeaderTest
  )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE
${TubeTK_SOURCE_DIR}/Projects/ImageClassWithPrior/Base/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeProjImageClassBase_SRCS
  tubeProjImageClassBasePrintTest.cxx
  itkFeatureGeneratingImageFunctionTest.cxx
  itkPatchFeatureGeneratingImageFunctionTest.cxx
  itkStandardFeatureGeneratingImageFunctionTest.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  ${TubeTK_SOURCE_DIR}/Projects/ImageClassWithPrior/Base
  )

add_executable( tubeProjImageClassBaseHeaderTest
  tubeProjImageClassBaseHeaderTest.cxx
  )
target_link_libraries( tubeProjImageClassBaseHeaderTest
  TubeNumerics
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_executable( tubeProjImageClassBaseTests
  tubeProjImageClassBaseTests.cxx
  ${tubeProjImageClassBase_SRCS}
  )
target_link_libraries( tubeProjImageClassBaseTests
  TubeNumerics
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_test( tubeProjImageClassBaseHeaderTest
  ${BASE_CLASSIFICATION_HEADER_TEST}
  )

add_test( tubeProjImageClassBasePrintTest
  ${BASE_CLASSIFICATION_TESTS}
  tubeProjImageClassBasePrintTest
  )

add_test( itkFeatureGeneratingImageFunctionTest1
  ${BASE_CLASSIFICATION_TESTS}
  itkFeatureGeneratingImageFunctionTest
  ${TubeTK_SOURCE_DIR}/Data/ES0015_1.mha
  )
add_test( itkStandardFeatureGeneratingImageFunctionTest1
  ${BASE_CLASSIFICATION_TESTS}
  itkStandardFeatureGeneratingImageFunctionTest
  ${TubeTK_SOURCE_DIR}/Data/GDS0015_Large_match.mha
  ${TubeTK_SOURCE_DIR}/Data/ES0015_Large.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_addMeanHist.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_addStdDevHist.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_nomMeanHist.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_nomStdDevHist.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_subMeanHist.mha
  ${TubeTK_SOURCE_DIR}/Data/Modified_GDS0015_Large_subStdDevHist.mha
  )

add_test( itkPatchFeatureGeneratingImageFunctionTest1
  ${BASE_CLASSIFICATION_TESTS}
  itkPatchFeatureGeneratingImageFunctionTest
  ${TubeTK_SOURCE_DIR}/Data/ES0015_1.mha
  ${TubeTK_SOURCE_DIR}/Data/GDS0015_1_match.mha
  )
