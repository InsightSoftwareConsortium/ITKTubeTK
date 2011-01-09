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
  add_test( PrintSelf-BaseSegmentation ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl
    ${TubeTK_SOURCE_DIR}/Base/Segmentation
    )
endif( TCL_TCLSH )

set( BASE_SEGMENTATION_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseSegmentationTests
  )

set( BASE_SEGMENTATION_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseSegmentationHeaderTest
  )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Segmentation/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseSegmentation_SRCS
  tubeBaseSegmentationPrintTest.cxx
  itkTubeCVTImageFilterTest.cxx
  itkTubeLDAGeneratorTest.cxx
  itkTubePDFSegmenterTest.cxx
  itkTubeRadiusExtractorTest.cxx
  itkTubeRadiusExtractorTest2.cxx
  itkTubeRidgeExtractorTest.cxx
  itkTubeRidgeExtractorTest2.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  )

add_executable( tubeBaseSegmentationHeaderTest
  tubeBaseSegmentationHeaderTest.cxx
  )
target_link_libraries( tubeBaseSegmentationHeaderTest
  TubeNumerics
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_executable( tubeBaseSegmentationTests
  tubeBaseSegmentationTests.cxx
  ${tubeBaseSegmentation_SRCS}
  )
target_link_libraries( tubeBaseSegmentationTests
  TubeNumerics
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_test( tubeBaseSegmentationHeaderTest
  ${BASE_SEGMENTATION_HEADER_TEST}
  )

add_test( tubeBaseSegmentationPrintTest
  ${BASE_SEGMENTATION_TESTS}
  tubeBaseSegmentationPrintTest
  )

add_test( itkTubeCVTImageFilterTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeCVTImageFilterTest
    ${INPUT_DATA}/GDS0015_1.mha
    ${TEMP}/itkTubeCVTImageFilterTest.mha
  )

add_test( itkTubeLDAGeneratorTest
  ${BASE_SEGMENTATION_TESTS}
  --compare ${BASELINE}/itkTubeLDAGeneratorTest_lda0.mha
            ${TEMP}/itkTubeLDAGeneratorTest_lda0.mha
  itkTubeLDAGeneratorTest
    ${INPUT_DATA}/ES0015_Large.mha
    ${INPUT_DATA}/GDS0015_Large-TrainingMask.mha
    ${TEMP}/itkTubeLDAGeneratorTest_lda0.mha
    ${TEMP}/itkTubeLDAGeneratorTest_lda1.mha
  )

add_test( itkTubePDFSegmenterTest
  ${BASE_SEGMENTATION_TESTS}
  --compare ${BASELINE}/itkTubePDFSegmenterTest_mask.mha
            ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  itkTubePDFSegmenterTest
    ${INPUT_DATA}/ES0015_Large.mha
    false
    ${INPUT_DATA}/GDS0015_Large-TrainingMask.mha
    ${TEMP}/itkTubePDFSegmenterTest_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  )

add_test( itkTubePDFSegmenterTest2
  ${BASE_SEGMENTATION_TESTS}
  --compare ${BASELINE}/itkTubePDFSegmenterTest2_mask.mha
            ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  itkTubePDFSegmenterTest
    ${INPUT_DATA}/ES0015_Large.mha
    true
    ${INPUT_DATA}/GDS0015_Large-TrainingMask.mha
    ${TEMP}/itkTubePDFSegmenterTest2_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest2_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  )

add_test( itkTubeRidgeExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  --compare ${BASELINE}/itkTubeRidgeExtractorTest.mha
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  itkTubeRidgeExtractorTest
    ${INPUT_DATA}/Branch.n010.sub.mha
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  )
 
add_test( itkTubeRidgeExtractorTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRidgeExtractorTest2
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )
 
add_test( itkTubeRadiusExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )
 
add_test( itkTubeRadiusExtractorTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest2
    ${INPUT_DATA}/Branch.n010.mha
    ${INPUT_DATA}/Branch-truth.tre
  )
