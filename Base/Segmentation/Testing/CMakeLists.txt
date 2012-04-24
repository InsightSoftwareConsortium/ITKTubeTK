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

Include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)


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
  itkTubeOtsuThresholdMaskedImageFilterTest.cxx
  itkTubePDFSegmenterTest.cxx
  itkTubeRadiusExtractorTest.cxx
  itkTubeRadiusExtractorTest2.cxx
  itkTubeRidgeExtractorTest.cxx
  itkTubeRidgeExtractorTest2.cxx
  itkTubeTubeExtractorTest.cxx
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

add_test( itkTubeOtsuThresholdMaskedImageFilterTest
  ${BASE_SEGMENTATION_TESTS}
  --compare ${BASELINE}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
            ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
  itkTubeOtsuThresholdMaskedImageFilterTest
    ${INPUT_DATA}/ES0015_1.mha
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest1.mha
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
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
    ${INPUT_DATA}/Branch.n010.sub.mha
    ${INPUT_DATA}/Branch-truth_Subs.tre
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

add_test( itkTubeTubeExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeTubeExtractorTest
    ${INPUT_DATA}/Branch.n010.sub.mha
    ${INPUT_DATA}/Branch-truth.tre
  )

## MIDAS TESTs
MIDAS_ADD_TEST( itkTubeCVTImageFilterMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeCVTImageFilterTest
    MIDAS{GDS0015_1.mha.md5}
    ${TEMP}/itkTubeCVTImageFilterTest.mha
  )

MIDAS_ADD_TEST( itkTubeOtsuThresholdMaskedImageFilterMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubeOtsuThresholdMaskedImageFilterTest2.mha.md5}
            ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
  itkTubeOtsuThresholdMaskedImageFilterTest
    MIDAS{ES0015_1.mha.md5}
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest1.mha
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
  )

MIDAS_ADD_TEST( itkTubePDFSegmenterMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubePDFSegmenterTest_mask.mha.md5}
            ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  itkTubePDFSegmenterTest
    MIDAS{ES0015_Large.mha.md5}
    false
    MIDAS{GDS0015_Large-TrainingMask.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  )

MIDAS_ADD_TEST( itkTubePDFSegmenterMIDASTest2
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubePDFSegmenterTest2_mask.mha.md5}
            ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  itkTubePDFSegmenterTest
    MIDAS{ES0015_Large.mha.md5}
    true
    MIDAS{GDS0015_Large-TrainingMask.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest2_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest2_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  )

MIDAS_ADD_TEST( itkTubeRidgeExtractorMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubeRidgeExtractorTest.mha.md5}
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  itkTubeRidgeExtractorTest
    MIDAS{Branch.n010.sub.mha.md5}
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  )

MIDAS_ADD_TEST( itkTubeRidgeExtractorMIDASTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRidgeExtractorTest2
    MIDAS{Branch.n010.sub.mha.md5}
    MIDAS{Branch-truth_Subs.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeRadiusExtractorMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest
    MIDAS{Branch.n010.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeRadiusExtractorMIDASTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest2
    MIDAS{Branch.n010.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeTubeExtractorMIDASTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeTubeExtractorTest
    MIDAS{Branch.n010.sub.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )
