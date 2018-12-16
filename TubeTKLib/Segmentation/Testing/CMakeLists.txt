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

include_regular_expression( "^.*$" )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-TubeTKLibSegmentation ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/TubeTKLib/Segmentation )
endif()

set( BASE_SEGMENTATION_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseSegmentationTests> )

set( BASE_SEGMENTATION_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseSegmentationHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseSegmentationTest_SRCS
  tubeBaseSegmentationTests.cxx
  tubeBaseSegmentationPrintTest.cxx
  itktubePDFSegmenterParzenTest.cxx
  itktubeRadiusExtractor2Test.cxx
  itktubeRadiusExtractor2Test2.cxx
  itktubeRidgeExtractorTest.cxx
  itktubeRidgeExtractorTest2.cxx
  itktubeRidgeSeedFilterTest.cxx
  itktubeTubeExtractorTest.cxx )

if( TubeTK_USE_LIBSVM )
  list( APPEND tubeBaseSegmentationTest_SRCS 
    itktubePDFSegmenterSVMTest.cxx )
endif()

if( TubeTK_USE_RANDOMFOREST )
  list( APPEND tubeBaseSegmentationTest_SRCS 
    itktubePDFSegmenterRandomForestTest.cxx )
endif()

include_directories(
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Common
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Numerics
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Filtering
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Segmentation )

add_executable( tubeBaseSegmentationHeaderTest
  tubeBaseSegmentationHeaderTest.cxx )
target_link_libraries( tubeBaseSegmentationHeaderTest
  TubeTKLibSegmentation TubeTKLibFiltering TubeTKLibNumerics TubeTKLibIO
  TubeTKLibCommon
  ${ITK_LIBRARIES} )
if( TubeTK_USE_LIBSVM )
  target_link_libraries( tubeBaseSegmentationHeaderTest
    ${LIBSVM_LIBRARIES} )
endif()

set( tubeBaseSegmentationTests_External_Libraries )
if( TubeTK_USE_LIBSVM )
  list( APPEND tubeBaseSegmentationTests_External_Libraries
    ${LIBSVM_LIBRARIES} )
endif()

add_executable( tubeBaseSegmentationTests 
  ${tubeBaseSegmentationTest_SRCS} )
target_link_libraries( tubeBaseSegmentationTests 
  TubeTKLibSegmentation TubeTKLibFiltering TubeTKLibNumerics
  TubeTKLibIO TubeTKLibCommon ${ITK_LIBRARIES} ITKIOMeta
  ITKIOSpatialObjects
  ${tubeBaseSegmentationTests_External_Libraries} )

add_test( NAME tubeBaseSegmentationHeaderTest
  COMMAND ${BASE_SEGMENTATION_HEADER_TEST} )

add_test( NAME tubeBaseSegmentationPrintTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
  tubeBaseSegmentationPrintTest )

ExternalData_Add_Test( TubeTKData
  NAME itktubePDFSegmenterParzenTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compareNumberOfPixelsTolerance 50
    --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterParzenTest_mask.mha}
      ${TEMP}/itktubePDFSegmenterParzenTest_mask.mha
    itktubePDFSegmenterParzenTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      true
      5.0
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
      ${TEMP}/itktubePDFSegmenterParzenTest_prob0.mha
      ${TEMP}/itktubePDFSegmenterParzenTest_pdf0.mha
      ${TEMP}/itktubePDFSegmenterParzenTest_prob1.mha
      ${TEMP}/itktubePDFSegmenterParzenTest_pdf1.mha
      ${TEMP}/itktubePDFSegmenterParzenTest_mask.mha
      ${TEMP}/itktubePDFSegmenterParzenTest_labeledFeatureSpace.mha )

ExternalData_Add_Test( TubeTKData
  NAME itktubePDFSegmenterParzenTest2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compareNumberOfPixelsTolerance 50
    --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterParzenTest2_mask.mha}
      ${TEMP}/itktubePDFSegmenterParzenTest2_mask.mha
    itktubePDFSegmenterParzenTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      false
      0.4
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
      ${TEMP}/itktubePDFSegmenterParzenTest2_prob0.mha
      ${TEMP}/itktubePDFSegmenterParzenTest2_pdf0.mha
      ${TEMP}/itktubePDFSegmenterParzenTest2_prob1.mha
      ${TEMP}/itktubePDFSegmenterParzenTest2_pdf1.mha
      ${TEMP}/itktubePDFSegmenterParzenTest2_mask.mha
      ${TEMP}/itktubePDFSegmenterParzenTest2_labeledFeatureSpace.mha )

if( TubeTK_USE_LIBSVM )

  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterSVMTest
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterSVMTest_mask.mha}
        ${TEMP}/itktubePDFSegmenterSVMTest_mask.mha
      itktubePDFSegmenterSVMTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        true
        5.0
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${TEMP}/itktubePDFSegmenterSVMTest_prob0.mha
        ${TEMP}/itktubePDFSegmenterSVMTest_prob1.mha
        ${TEMP}/itktubePDFSegmenterSVMTest_mask.mha )

  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterSVMTest2
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterSVMTest2_mask.mha}
        ${TEMP}/itktubePDFSegmenterSVMTest2_mask.mha
      itktubePDFSegmenterSVMTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        false
        0.4
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
        ${TEMP}/itktubePDFSegmenterSVMTest2_prob0.mha
        ${TEMP}/itktubePDFSegmenterSVMTest2_prob1.mha
        ${TEMP}/itktubePDFSegmenterSVMTest2_mask.mha )

endif()

if( TubeTK_USE_RANDOMFOREST )

  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterRandomForestTest
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterRandomForestTest_mask.mha}
        ${TEMP}/itktubePDFSegmenterRandomForestTest_mask.mha
      itktubePDFSegmenterRandomForestTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        true
        5.0
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${TEMP}/itktubePDFSegmenterRandomForestTest_prob0.mha
        ${TEMP}/itktubePDFSegmenterRandomForestTest_prob1.mha
        ${TEMP}/itktubePDFSegmenterRandomForestTest_mask.mha )

  ExternalData_Add_Test( TubeTKData
    NAME itktubePDFSegmenterRandomForestTest2
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterRandomForestTest2_mask.mha}
        ${TEMP}/itktubePDFSegmenterRandomForestTest2_mask.mha
      itktubePDFSegmenterRandomForestTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        false
        0.4
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
        ${TEMP}/itktubePDFSegmenterRandomForestTest2_prob0.mha
        ${TEMP}/itktubePDFSegmenterRandomForestTest2_prob1.mha
        ${TEMP}/itktubePDFSegmenterRandomForestTest2_mask.mha )

endif()

ExternalData_Add_Test( TubeTKData
  NAME itktubeRidgeExtractorTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeExtractorTest.mha}
      ${TEMP}/itktubeRidgeExtractorTest.mha
    itktubeRidgeExtractorTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      ${TEMP}/itktubeRidgeExtractorTest.mha )

ExternalData_Add_Test( TubeTKData
  NAME itktubeRidgeExtractorTest2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRidgeExtractorTest2
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth_Subs.tre} )

ExternalData_Add_Test( TubeTKData
  NAME itktubeRadiusExtractor2Test
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRadiusExtractor2Test
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

ExternalData_Add_Test( TubeTKData
  NAME itktubeRadiusExtractor2Test2
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeRadiusExtractor2Test2
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

ExternalData_Add_Test( TubeTKData
  NAME itktubeTubeExtractorTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    itktubeTubeExtractorTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

ExternalData_Add_Test( TubeTKData
  NAME itktubeRidgeSeedFilterParzenTest
  COMMAND ${BASE_SEGMENTATION_TESTS}
    --compareNumberOfPixelsTolerance 100
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterParzenTest_Output.mha}
      ${TEMP}/itktubeRidgeSeedFilterParzenTest_Output.mha
    itktubeRidgeSeedFilterTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
      255 127 0
      ${TEMP}/itktubeRidgeSeedFilterParzenTest_Feature0Image.mha
      ${TEMP}/itktubeRidgeSeedFilterParzenTest_Output.mha
      ${TEMP}/itktubeRidgeSeedFilterParzenTest_MaxScale.mha )

if( TubeTK_USE_LIBSVM )

  ExternalData_Add_Test( TubeTKData
    NAME itktubeRidgeSeedFilterSVMTest
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 700
      --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterSVMTest_Output.mha}
        ${TEMP}/itktubeRidgeSeedFilterSVMTest_Output.mha
      itktubeRidgeSeedFilterTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
        255 127 1
        ${TEMP}/itktubeRidgeSeedFilterSVMTest_Feature0Image.mha
        ${TEMP}/itktubeRidgeSeedFilterSVMTest_Output.mha
        ${TEMP}/itktubeRidgeSeedFilterSVMTest_MaxScale.mha )

endif()

if( TubeTK_USE_RANDOMFOREST )

  ExternalData_Add_Test( TubeTKData
    NAME itktubeRidgeSeedFilterRandomForestTest
    COMMAND ${BASE_SEGMENTATION_TESTS}
      --compareNumberOfPixelsTolerance 700
      --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterRandomForestTest_Output.mha}
        ${TEMP}/itktubeRidgeSeedFilterRandomForestTest_Output.mha
      itktubeRidgeSeedFilterTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
        255 127 2
        ${TEMP}/itktubeRidgeSeedFilterRandomForestTest_Feature0Image.mha
        ${TEMP}/itktubeRidgeSeedFilterRandomForestTest_Output.mha
        ${TEMP}/itktubeRidgeSeedFilterRandomForestTest_MaxScale.mha )

endif( TubeTK_USE_RANDOMFOREST )
