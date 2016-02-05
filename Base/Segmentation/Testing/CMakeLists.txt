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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseSegmentation ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Segmentation )
endif( PYTHON_EXECUTABLE )

set( BASE_SEGMENTATION_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseSegmentationTests> )

set( BASE_SEGMENTATION_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseSegmentationHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseSegmentation_SRCS
  tubeBaseSegmentationPrintTest.cxx
  itktubeCVTImageFilterTest.cxx
  itktubePDFSegmenterTest.cxx
  itktubeRadiusExtractor2Test.cxx
  itktubeRadiusExtractor2Test2.cxx
  itktubeRidgeExtractorTest.cxx
  itktubeRidgeExtractorTest2.cxx
  itktubeRidgeSeedFilterTest.cxx
  itktubeTubeExtractorTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation )

add_executable( tubeBaseSegmentationHeaderTest
  tubeBaseSegmentationHeaderTest.cxx )
target_link_libraries( tubeBaseSegmentationHeaderTest
  TubeTKSegmentation TubeTKFiltering TubeTKNumerics TubeTKIO TubeTKCommon
  ${ITK_LIBRARIES} )

SEMMacroBuildCLI(
  NAME tubeBaseSegmentationTests
  ADDITIONAL_SRCS
    ${tubeBaseSegmentation_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    TubeCLI TubeTKSegmentation TubeTKFiltering TubeTKNumerics TubeTKIO TubeTKCommon
    ${ITK_LIBRARIES} ITKIOMeta ITKIOSpatialObjects
  INSTALL_RUNTIME_DESTINATION bin
  EXECUTABLE_ONLY
  )

add_test( NAME tubeBaseSegmentationHeaderTest
  COMMAND ${BASE_SEGMENTATION_HEADER_TEST} )

add_test( NAME tubeBaseSegmentationPrintTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
  tubeBaseSegmentationPrintTest )

Midas3FunctionAddTest( NAME itktubeCVTImageFilterTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeCVTImageFilterTest
      MIDAS{GDS0015_1.mha.md5}
      ${TEMP}/itktubeCVTImageFilterTest.mha )

Midas3FunctionAddTest( NAME itktubePDFSegmenterTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compare MIDAS{itktubePDFSegmenterTest_mask.mha.md5}
      ${TEMP}/itktubePDFSegmenterTest_mask.mha
    itktubePDFSegmenterTest
      MIDAS{ES0015_Large.mha.md5}
      MIDAS{ES0015_Large.mha.md5}
      true
      5.0
      MIDAS{GDS0015_Large-TrainingMask.mha.md5}
      ${TEMP}/itktubePDFSegmenterTest_prob0.mha
      ${TEMP}/itktubePDFSegmenterTest_pdf0.mha
      ${TEMP}/itktubePDFSegmenterTest_prob1.mha
      ${TEMP}/itktubePDFSegmenterTest_pdf1.mha
      ${TEMP}/itktubePDFSegmenterTest_mask.mha
      ${TEMP}/itktubePDFSegmenterTest_labeledFeatureSpace.mha )

Midas3FunctionAddTest( NAME itktubePDFSegmenterTest2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compare MIDAS{itktubePDFSegmenterTest2_mask.mha.md5}
      ${TEMP}/itktubePDFSegmenterTest2_mask.mha
    itktubePDFSegmenterTest
      MIDAS{im0001.crop.mha.md5}
      MIDAS{im0001.crop.contrast.mha.md5}
      false
      0.4
      MIDAS{im0001.vk.mask.crop.mha.md5}
      ${TEMP}/itktubePDFSegmenterTest2_prob0.mha
      ${TEMP}/itktubePDFSegmenterTest2_pdf0.mha
      ${TEMP}/itktubePDFSegmenterTest2_prob1.mha
      ${TEMP}/itktubePDFSegmenterTest2_pdf1.mha
      ${TEMP}/itktubePDFSegmenterTest2_mask.mha
      ${TEMP}/itktubePDFSegmenterTest2_labeledFeatureSpace.mha )

Midas3FunctionAddTest( NAME itktubeRidgeExtractorTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compare MIDAS{itktubeRidgeExtractorTest.mha.md5}
      ${TEMP}/itktubeRidgeExtractorTest.mha
    itktubeRidgeExtractorTest
      MIDAS{Branch.n010.sub.mha.md5}
      ${TEMP}/itktubeRidgeExtractorTest.mha )

Midas3FunctionAddTest( NAME itktubeRidgeExtractorTest2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRidgeExtractorTest2
      MIDAS{Branch.n010.sub.mha.md5}
      MIDAS{Branch-truth_Subs.tre.md5} )

Midas3FunctionAddTest( NAME itktubeRadiusExtractor2Test
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRadiusExtractor2Test
      MIDAS{Branch.n010.mha.md5}
      MIDAS{Branch-truth.tre.md5} )

Midas3FunctionAddTest( NAME itktubeRadiusExtractor2Test2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRadiusExtractor2Test2
      MIDAS{Branch.n010.mha.md5}
      MIDAS{Branch-truth.tre.md5} )

Midas3FunctionAddTest( NAME itktubeTubeExtractorTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeTubeExtractorTest
      MIDAS{Branch.n010.sub.mha.md5}
      MIDAS{Branch-truth.tre.md5} )

Midas3FunctionAddTest( NAME itktubeRidgeSeedFilterTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compareNumberOfPixelsTolerance 100
    --compare MIDAS{itktubeRidgeSeedFilterTest_Output.mha.md5}
      ${TEMP}/itktubeRidgeSeedFilterTest_Output.mha
    itktubeRidgeSeedFilterTest
      MIDAS{im0001.crop.contrast.mha.md5}
      MIDAS{im0001.vk.maskRidge.crop.mha.md5}
      255 127
      ${TEMP}/itktubeRidgeSeedFilterTest_Feature0Image.mha
      ${TEMP}/itktubeRidgeSeedFilterTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterTest_MaxScale.mha )
