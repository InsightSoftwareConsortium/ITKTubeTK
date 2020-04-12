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

if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-Numerics ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/TubeTKLib/Numerics )
endif()

set( tubeNumericsTests_SRCS
  tubeNumericsPrintTest.cxx
  itktubeBlurImageFunctionTest.cxx
  itktubeImageRegionMomentsCalculatorTest.cxx
  itktubeJointHistogramImageFunctionTest.cxx
  itktubeNJetBasisFeatureVectorGeneratorTest.cxx
  itktubeNJetFeatureVectorGeneratorTest.cxx
  itktubeNJetImageFunctionTest.cxx
  itktubeSingleValuedCostFunctionImageSourceTest.cxx
  itktubeRecordOptimizationParameterProgressionCommandTest.cxx
  itktubeRidgeBasisFeatureVectorGeneratorTest.cxx
  itktubeRidgeFFTFeatureVectorGeneratorTest.cxx
  itktubeVotingResampleImageFunctionTest.cxx
  tubeBrentOptimizer1DTest.cxx
  tubeBrentOptimizerNDTest.cxx
  tubeGoldenMeanOptimizer1DTest.cxx
  tubeMatrixMathTest.cxx
  tubeParabolicFitOptimizer1DTest.cxx
  tubeSplineApproximation1DTest.cxx
  tubeSplineNDTest.cxx
  tubeUserFunctionTest.cxx )

CreateTestDriver( tubeNumericsHeader
  "${TubeTK-Test_LIBRARIES}"
  "tubeNumericsHeaderTest.cxx" )
if( TubeTK_USE_ARRAYFIRE )
  target_include_directories( tubeNumericsHeaderTestDriver PUBLIC
    ${ArrayFire_INCLUDE_DIRS} )
  target_link_libraries( tubeNumericsHeaderTestDriver LINK_PUBLIC
    ${ArrayFire_LIBRARIES} )
endif()
itk_add_test( NAME tubeNumericsHeaderTest
  COMMAND tubeNumericsHeaderTestDriver
    tubeNumericsHeaderTest )

CreateTestDriver( tubeNumerics
  "${TubeTK-Test_LIBRARIES}"
  "${tubeNumericsTests_SRCS}" )
if( TubeTK_USE_ARRAYFIRE )
  target_include_directories( tubeNumericsTestDriver PUBLIC
    ${ArrayFire_INCLUDE_DIRS} )
  target_link_libraries( tubeNumericsTestDriver LINK_PUBLIC
    ${ArrayFire_LIBRARIES} )
endif()

itk_add_test( NAME tubeNumericsPrintTest
  COMMAND tubeNumericsTestDriver
    tubeNumericsPrintTest )

itk_add_test( NAME tubeMatrixMathTest
  COMMAND tubeNumericsTestDriver
    tubeMatrixMathTest )

itk_add_test( NAME tubeBrentOptimizer1DTest
  COMMAND tubeNumericsTestDriver
    tubeBrentOptimizer1DTest )

# Minimization - initialized at middle of interval
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 0.01 -0.5 -1 )

# impossible maximization - not unimodal distribution
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest1
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 0 350 100 -1 )
set_tests_properties( tubeGoldenMeanOptimizer1DTest1
  PROPERTIES WILL_FAIL true )

# minimization - initialized at xMax
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest2
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 1 350 100 -1 )

# minimization - initialized at xMin
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest3
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 -3 -0.5 -1 )

# maximization - initialized at middle of interval
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest4
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 100 0 1 )

# maximization - initialized at xMax
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest5
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 300 0 1 )

# maximization - initialized at xMin
itk_add_test( NAME tubeGoldenMeanOptimizer1DTest6
  COMMAND tubeNumericsTestDriver
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 -300 0 1 )

itk_add_test( NAME tubeParabolicFitOptimizer1DTest
  COMMAND tubeNumericsTestDriver
    tubeParabolicFitOptimizer1DTest )

itk_add_test( NAME tubeBrentOptimizerNDTest
  COMMAND tubeNumericsTestDriver
    tubeBrentOptimizerNDTest )

itk_add_test( NAME tubeUserFunctionTest
  COMMAND tubeNumericsTestDriver
    tubeUserFunctionTest )

itk_add_test(
  NAME itktubeJointHistogramImageFunctionTest01
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itkJointHistogramImageFunctionTest01.mha}
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest01.mha
    --compareNumberOfPixelsTolerance 7
    itktubeJointHistogramImageFunctionTest
      DATA{${TubeTK_DATA_ROOT}/GDS0015_1_match_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_1_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_1_match_Subs.mask.mha}
      0
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest01.mha
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest01Mean.mha
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest01StdDev.mha )

itk_add_test(
  NAME itktubeJointHistogramImageFunctionTest02
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itkJointHistogramImageFunctionTest02.mha}
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest02.mha
    --compareNumberOfPixelsTolerance 7
    itktubeJointHistogramImageFunctionTest
      DATA{${TubeTK_DATA_ROOT}/GDS0015_1_match_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/ES0015_1_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_1_match_Subs.mask.mha}
      1
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest02.mha
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest02Mean.mha
      ${ITK_TEST_OUTPUT_DIR}/itkJointHistogramImageFunctionTest02StdDev.mha )

itk_add_test(
  NAME itktubeImageRegionMomentsCalculatorTest
  COMMAND tubeNumericsTestDriver
    itktubeImageRegionMomentsCalculatorTest
      DATA{${TubeTK_DATA_ROOT}/scoring-test.png} )

foreach( testNum RANGE 0 38 )
  itk_add_test(
    NAME itktubeNJetImageFunctionTest${testNum}
    COMMAND tubeNumericsTestDriver
      --compare DATA{${TubeTK_DATA_ROOT}/itktubeNJetImageFunctionTest${testNum}.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubeNJetImageFunctionTest${testNum}.mha
        --compareIntensityTolerance 0.000000000001
      itktubeNJetImageFunctionTest
        ${testNum}
        DATA{${TubeTK_DATA_ROOT}/im0001.crop2.mha}
        ${ITK_TEST_OUTPUT_DIR}/itktubeNJetImageFunctionTest${testNum}.mha )
endforeach( testNum )

itk_add_test(
  NAME itktubeVotingResampleImageFunctionTest0
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeVotingResampleImageFunctionTest0.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest0.png
    itktubeVotingResampleImageFunctionTest
      0
      DATA{${TubeTK_DATA_ROOT}/greyscale01.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest0.png )

itk_add_test(
  NAME itktubeVotingResampleImageFunctionTest1
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeVotingResampleImageFunctionTest1.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest1.png
    itktubeVotingResampleImageFunctionTest
      1
      DATA{${TubeTK_DATA_ROOT}/greyscale01.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest1.png )

itk_add_test(
  NAME itktubeVotingResampleImageFunctionTest2
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeVotingResampleImageFunctionTest2.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest2.png
    itktubeVotingResampleImageFunctionTest
      2
      DATA{${TubeTK_DATA_ROOT}/greyscale01.png}
      ${ITK_TEST_OUTPUT_DIR}/itktubeVotingResampleImageFunctionTest2.png )

itk_add_test(
  NAME tubeSplineApproximation1DTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/tubeSplineApprox1DTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/tubeSplineApprox1DTest.mha
    tubeSplineApproximation1DTest
      ${ITK_TEST_OUTPUT_DIR}/tubeSplineApprox1DTest.mha )

itk_add_test(
  NAME tubeSplineNDTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/tubeSplineNDTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/tubeSplineNDTest.mha
    tubeSplineNDTest
      ${ITK_TEST_OUTPUT_DIR}/tubeSplineNDTest.mha )

itk_add_test( NAME itktubeRecordOptimizationParameterProgressionCommandTest
  COMMAND tubeNumericsTestDriver
    itktubeRecordOptimizationParameterProgressionCommandTest
      ${ITK_TEST_OUTPUT_DIR}/itkRecordOptimizationParameterProgressionCommandTest.h5 )

itk_add_test(
  NAME itktubeBlurImageFunctionTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeBlurImageFunctionTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeBlurImageFunctionTest.mha
    itktubeBlurImageFunctionTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeBlurImageFunctionTest.mha )

itk_add_test(
  NAME itktubeNJetFeatureVectorGeneratorTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeNJetFeatureVectorGeneratorTest_f0.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetFeatureVectorGeneratorTest_f0.mha
    itktubeNJetFeatureVectorGeneratorTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask_Subs.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetFeatureVectorGeneratorTest_f0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetFeatureVectorGeneratorTest_f1.mha )

itk_add_test(
  NAME itktubeNJetBasisFeatureVectorGeneratorTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeNJetBasisFeatureVectorGeneratorTest_basis.0.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetBasisFeatureVectorGeneratorTest_basis.0.mha
    itktubeNJetBasisFeatureVectorGeneratorTest
      DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Subs.mha}
      DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask_Subs.mha}
      255 0
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetBasisFeatureVectorGeneratorTest_basis
      ${ITK_TEST_OUTPUT_DIR}/itktubeNJetBasisFeatureVectorGeneratorTest_feature )

itk_add_test(
  NAME itktubeSingleValuedCostFunctionImageSourceTest
  COMMAND tubeNumericsTestDriver
    --compare DATA{${TubeTK_DATA_ROOT}/itktubeSingleValuedCostFunctionImageSourceTest.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeSingleValuedCostFunctionImageSourceTest.mha
    itktubeSingleValuedCostFunctionImageSourceTest
      ${ITK_TEST_OUTPUT_DIR}/itktubeSingleValuedCostFunctionImageSourceTest.mha )

itk_add_test(
  NAME itktubeRidgeFFTFeatureVectorGeneratorTest
  COMMAND tubeNumericsTestDriver
    itktubeRidgeFFTFeatureVectorGeneratorTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFeatureVectorGeneratorTest_f0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeFFTFeatureVectorGeneratorTest_f1.mha )

itk_add_test(
  NAME itktubeRidgeBasisFeatureVectorGeneratorTest
  COMMAND tubeNumericsTestDriver
    itktubeRidgeBasisFeatureVectorGeneratorTest
      DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
      DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
      255 127
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeBasisFeatureVectorGeneratorTest_lda0.mha
      ${ITK_TEST_OUTPUT_DIR}/itktubeRidgeBasisFeatureVectorGeneratorTest_lda1.mha )
