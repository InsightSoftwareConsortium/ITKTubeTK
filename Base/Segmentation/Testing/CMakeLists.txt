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

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseSegmentation ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Segmentation
    )
endif()

set( BASE_SEGMENTATION_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseSegmentationTests
  )

set( BASE_SEGMENTATION_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseSegmentationHeaderTest
  )

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
  itkTubeRidgeSeedGeneratorTest.cxx
  itkTubeTubeExtractorTest.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  )

add_executable( tubeBaseSegmentationHeaderTest
  tubeBaseSegmentationHeaderTest.cxx
  )
target_link_libraries( tubeBaseSegmentationHeaderTest
  TubeNumerics
  ${ITK_LIBRARIES}
  )

add_executable( tubeBaseSegmentationTests
  tubeBaseSegmentationTests.cxx
  ${tubeBaseSegmentation_SRCS}
  )
target_link_libraries( tubeBaseSegmentationTests
  TubeNumerics
  ${ITK_LIBRARIES}
  )

add_test( tubeBaseSegmentationHeaderTest
  ${BASE_SEGMENTATION_HEADER_TEST}
  )

add_test( tubeBaseSegmentationPrintTest
  ${BASE_SEGMENTATION_TESTS}
  tubeBaseSegmentationPrintTest
  )

MIDAS_ADD_TEST( itkTubeCVTImageFilterTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeCVTImageFilterTest
    MIDAS{GDS0015_1.mha.md5}
    ${TEMP}/itkTubeCVTImageFilterTest.mha
  )

MIDAS_ADD_TEST( itkTubeOtsuThresholdMaskedImageFilterTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubeOtsuThresholdMaskedImageFilterTest2.mha.md5}
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
  itkTubeOtsuThresholdMaskedImageFilterTest
    MIDAS{ES0015_1.mha.md5}
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest1.mha
    ${TEMP}/itkTubeOtsuThresholdMaskedImageFilterTest2.mha
  )

MIDAS_ADD_TEST( itkTubePDFSegmenterTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubePDFSegmenterTest_mask.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  itkTubePDFSegmenterTest
    MIDAS{ES0015_Large.mha.md5}
    MIDAS{ES0015_Large.mha.md5}
    true
    5.0
    MIDAS{GDS0015_Large-TrainingMask.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest_pdf0.mha
    ${TEMP}/itkTubePDFSegmenterTest_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest_pdf1.mha
    ${TEMP}/itkTubePDFSegmenterTest_mask.mha
  )

MIDAS_ADD_TEST( itkTubePDFSegmenterTest2
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubePDFSegmenterTest2_mask.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  itkTubePDFSegmenterTest
    MIDAS{im0001.crop.mha.md5}
    MIDAS{im0001.crop.contrast.mha.md5}
    false
    0.4
    MIDAS{im0001.vk.mask.crop.mha.md5}
    ${TEMP}/itkTubePDFSegmenterTest2_prob0.mha
    ${TEMP}/itkTubePDFSegmenterTest2_pdf0.mha
    ${TEMP}/itkTubePDFSegmenterTest2_prob1.mha
    ${TEMP}/itkTubePDFSegmenterTest2_pdf1.mha
    ${TEMP}/itkTubePDFSegmenterTest2_mask.mha
  )

MIDAS_ADD_TEST( itkTubeRidgeExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubeRidgeExtractorTest.mha.md5}
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  itkTubeRidgeExtractorTest
    MIDAS{Branch.n010.sub.mha.md5}
    ${TEMP}/itkTubeRidgeExtractorTest.mha
  )

MIDAS_ADD_TEST( itkTubeRidgeExtractorTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRidgeExtractorTest2
    MIDAS{Branch.n010.sub.mha.md5}
    MIDAS{Branch-truth_Subs.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeRadiusExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest
    MIDAS{Branch.n010.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeRadiusExtractorTest2
  ${BASE_SEGMENTATION_TESTS}
  itkTubeRadiusExtractorTest2
    MIDAS{Branch.n010.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeTubeExtractorTest
  ${BASE_SEGMENTATION_TESTS}
  itkTubeTubeExtractorTest
    MIDAS{Branch.n010.sub.mha.md5}
    MIDAS{Branch-truth.tre.md5}
  )

MIDAS_ADD_TEST( itkTubeRidgeSeedGeneratorTest
  ${BASE_SEGMENTATION_TESTS}
  --compare MIDAS{itkTubeRidgeSeedGeneratorTest.lda00.mha.md5}
    ${TEMP}/itkTubeRidgeSeedGeneratorTest.lda00.mha
  itkTubeRidgeSeedGeneratorTest
    MIDAS{im0001.crop.contrast.mha.md5}
    MIDAS{im0001.vk.mask.crop.mha.md5}
    ${TEMP}/itkTubeRidgeSeedGeneratorTest
  )
