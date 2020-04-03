##############################################################################
#
# Library:   TubeTK
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

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  itk_add_test( PrintSelf-TubeTKSegmentation ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/src/Segmentation )
endif()

CreateTestDriver( tubeSegmentationHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeSegmentationHeaderTest.cxx" )
itk_add_test( NAME tubeSegmentationHeaderTest
  COMMAND tubeSegmentationHeaderTestDriver
    tubeSegmentationHeaderTest )

set( tubeSegmentationTest_SRCS
  tubeSegmentationPrintTest.cxx
  itktubePDFSegmenterParzenTest.cxx
  itktubeRadiusExtractor2Test.cxx
  itktubeRadiusExtractor2Test2.cxx
  itktubeRidgeExtractorTest.cxx
  itktubeRidgeExtractorTest2.cxx
  itktubeRidgeSeedFilterTest.cxx
  itktubeTubeExtractorTest.cxx )

#if( TubeTK_USE_LIBSVM )
  #list( APPEND tubeSegmentationTest_SRCS 
    #itktubePDFSegmenterSVMTest.cxx )
#endif()

#if( TubeTK_USE_RANDOMFOREST )
  #list( APPEND tubeSegmentationTest_SRCS 
    #itktubePDFSegmenterRandomForestTest.cxx )
#endif()

CreateTestDriver( tubeSegmentation
  "${TubeTK-Test_LIBRARIES}"
  "${tubeSegmentationTest_SRCS}" )

itk_add_test( NAME tubeSegmentationPrintTest
  COMMAND tubeSegmentationTestDriver
  tubeSegmentationPrintTest )

itk_add_test(
  NAME itktubePDFSegmenterParzenTest
  COMMAND tubeSegmentationTestDriver
    --compareNumberOfPixelsTolerance 50
    --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterParzenTest_mask.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_mask.mha
    itktubePDFSegmenterParzenTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
      true
      5.0
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_prob0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_pdf0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_prob1.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_pdf1.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_mask.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest_labeledFeatureSpace.mha )

itk_add_test(
  NAME itktubePDFSegmenterParzenTest2
  COMMAND tubeSegmentationTestDriver
    --compareNumberOfPixelsTolerance 50
    --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterParzenTest2_mask.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_mask.mha
    itktubePDFSegmenterParzenTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      false
      0.4
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_prob0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_pdf0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_prob1.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_pdf1.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_mask.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterParzenTest2_labeledFeatureSpace.mha )

if( TubeTK_USE_LIBSVM )

  itk_add_test(
    NAME itktubePDFSegmenterSVMTest
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterSVMTest_mask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest_mask.mha
      itktubePDFSegmenterSVMTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        true
        5.0
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest_prob0.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest_prob1.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest_mask.mha )

  itk_add_test(
    NAME itktubePDFSegmenterSVMTest2
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterSVMTest2_mask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest2_mask.mha
      itktubePDFSegmenterSVMTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        false
        0.4
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest2_prob0.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest2_prob1.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterSVMTest2_mask.mha )

endif()

if( TubeTK_USE_RANDOMFOREST )

  itk_add_test(
    NAME itktubePDFSegmenterRandomForestTest
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterRandomForestTest_mask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest_mask.mha
      itktubePDFSegmenterRandomForestTest
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        true
        5.0
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest_prob0.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest_prob1.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest_mask.mha )

  itk_add_test(
    NAME itktubePDFSegmenterRandomForestTest2
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 300
      --compare DATA{${TubeTK_DATA_ROOT}/itktubePDFSegmenterRandomForestTest2_mask.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest2_mask.mha
      itktubePDFSegmenterRandomForestTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        false
        0.4
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest2_prob0.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest2_prob1.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubePDFSegmenterRandomForestTest2_mask.mha )

endif()

itk_add_test(
  NAME itktubeRidgeExtractorTest
  COMMAND tubeSegmentationTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeExtractorTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeExtractorTest.mha
    itktubeRidgeExtractorTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeExtractorTest.mha )

itk_add_test(
  NAME itktubeRidgeExtractorTest2
  COMMAND tubeSegmentationTestDriver
    itktubeRidgeExtractorTest2
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth_Subs.tre} )

itk_add_test(
  NAME itktubeRadiusExtractor2Test
  COMMAND tubeSegmentationTestDriver
    itktubeRadiusExtractor2Test
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

itk_add_test(
  NAME itktubeRadiusExtractor2Test2
  COMMAND tubeSegmentationTestDriver
    itktubeRadiusExtractor2Test2
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

itk_add_test(
  NAME itktubeTubeExtractorTest
  COMMAND tubeSegmentationTestDriver
    itktubeTubeExtractorTest
      DATA{${TubeTK_DATA_ROOT}/Branch.n010.sub.mha}
      DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre} )

itk_add_test(
  NAME itktubeRidgeSeedFilterParzenTest
  COMMAND tubeSegmentationTestDriver
    --compareNumberOfPixelsTolerance 100
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterParzenTest_Output.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterParzenTest_Output.mha
    itktubeRidgeSeedFilterTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
      255 127 0
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterParzenTest_Feature0Image.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterParzenTest_Output.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterParzenTest_MaxScale.mha )

if( TubeTK_USE_LIBSVM )

  itk_add_test(
    NAME itktubeRidgeSeedFilterSVMTest
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 700
      --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterSVMTest_Output.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterSVMTest_Output.mha
      itktubeRidgeSeedFilterTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
        255 127 1
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterSVMTest_Feature0Image.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterSVMTest_Output.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterSVMTest_MaxScale.mha )

endif()

if( TubeTK_USE_RANDOMFOREST )

  itk_add_test(
    NAME itktubeRidgeSeedFilterRandomForestTest
    COMMAND tubeSegmentationTestDriver
      --compareNumberOfPixelsTolerance 700
      --compare DATA{${TubeTK_DATA_ROOT}/itktubeRidgeSeedFilterRandomForestTest_Output.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterRandomForestTest_Output.mha
      itktubeRidgeSeedFilterTest
        DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
        DATA{${TubeTK_DATA_ROOT}/im0001.vk.maskRidge.crop.mha}
        255 127 2
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterRandomForestTest_Feature0Image.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterRandomForestTest_Output.mha
        ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeSeedFilterRandomForestTest_MaxScale.mha )

endif( TubeTK_USE_RANDOMFOREST )
