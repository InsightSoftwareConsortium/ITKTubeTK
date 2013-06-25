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

include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseNumerics ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Numerics )
endif()

set( BASE_NUMERICS_TESTS
      ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseNumericsTests )

set( BASE_NUMERICS_HEADER_TEST
     ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseNumericsHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseNumerics_SRCS
  tubeBaseNumericsPrintTest.cxx
  itktubeBlurImageFunctionTest.cxx
  itktubeImageRegionMomentsCalculatorTest.cxx
  itktubeJointHistogramImageFunctionTest.cxx
  itktubeNJetBasisFeatureVectorGeneratorTest.cxx
  itktubeNJetFeatureVectorGeneratorTest.cxx
  itktubeNJetImageFunctionTest.cxx
  itktubeRecordOptimizationParameterProgressionCommandTest.cxx
  itktubeRidgeBasisFeatureVectorGeneratorTest.cxx
  itktubeRidgeFeatureVectorGeneratorTest.cxx
  itktubeVotingResampleImageFunctionTest.cxx
  tubeBrentOptimizer1DTest.cxx
  tubeBrentOptimizerNDTest.cxx
  tubeGoldenMeanOptimizer1DTest.cxx
  tubeMatrixMathTest.cxx
  tubeParabolicFitOptimizer1DTest.cxx
  tubeSplineApproximation1DTest.cxx
  tubeSplineNDTest.cxx
  tubeUserFunctionTest.cxx )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  ${TubeTK_SOURCE_DIR}/Base/Numerics )

add_executable( tubeBaseNumericsHeaderTest
  tubeBaseNumericsHeaderTest.cxx )
target_link_libraries( tubeBaseNumericsHeaderTest
  TubeNumerics
  ${ITK_LIBRARIES} )

add_executable( tubeBaseNumericsTests
  tubeBaseNumericsTests.cxx
  ${tubeBaseNumerics_SRCS} )
target_link_libraries( tubeBaseNumericsTests
  TubeNumerics
  ${ITK_LIBRARIES} )

add_test( tubeBaseNumericsHeaderTest
  ${BASE_NUMERICS_HEADER_TEST} )

add_test( tubeBaseNumericsPrintTest
  ${BASE_NUMERICS_TESTS}
    tubeBaseNumericsPrintTest )

add_test( tubeMatrixMathTest
  ${BASE_NUMERICS_TESTS}
    tubeMatrixMathTest )

add_test( tubeBrentOptimizer1DTest
  ${BASE_NUMERICS_TESTS}
    tubeBrentOptimizer1DTest )

# Minimization - initialized at middle of interval
add_test( tubeGoldenMeanOptimizer1DTest
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 0.01 -0.5 -1 )

# impossible maximization - not unimodal distribution
add_test( tubeGoldenMeanOptimizer1DTest1
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 0 350 100 -1 )
set_tests_properties( tubeGoldenMeanOptimizer1DTest1
  PROPERTIES WILL_FAIL true )

# minimization - initialized at xMax
add_test( tubeGoldenMeanOptimizer1DTest2
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 1 350 100 -1 )

# minimization - initialized at xMin
add_test( tubeGoldenMeanOptimizer1DTest3
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 -3 -0.5 -1 )

# maximization - initialized at middle of interval
add_test( tubeGoldenMeanOptimizer1DTest4
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 100 0 1 )

# maximization - initialized at xMax
add_test( tubeGoldenMeanOptimizer1DTest5
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 300 0 1 )

# maximization - initialized at xMin
add_test( tubeGoldenMeanOptimizer1DTest6
  ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 -300 0 1 )

add_test( tubeParabolicFitOptimizer1DTest
  ${BASE_NUMERICS_TESTS}
    tubeParabolicFitOptimizer1DTest )

add_test( tubeBrentOptimizerNDTest
  ${BASE_NUMERICS_TESTS}
    tubeBrentOptimizerNDTest )

add_test( tubeUserFunctionTest
  ${BASE_NUMERICS_TESTS}
    tubeUserFunctionTest )

MIDAS_ADD_TEST( itktubeJointHistogramImageFunctionTest01
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest01.mha.md5}
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
    itktubeJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      0
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01StdDev.mha )

MIDAS_ADD_TEST( itktubeJointHistogramImageFunctionTest02
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest02.mha.md5}
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
    itktubeJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      1
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02StdDev.mha )

MIDAS_ADD_TEST( itktubeImageRegionMomentsCalculatorTest
  ${BASE_NUMERICS_TESTS}
    itktubeImageRegionMomentsCalculatorTest
      MIDAS{scoring-test.png.md5} )

foreach( testNum RANGE 0 37 )
  MIDAS_ADD_TEST( itktubeNJetImageFunctionTest${testNum}
    ${BASE_NUMERICS_TESTS}
      --compare MIDAS{itkTubeNJetImageFunctionTest${testNum}.mha.md5}
        ${TEMP}/itkTubeNJetImageFunctionTest${testNum}.mha
      itktubeNJetImageFunctionTest
        ${testNum}
        MIDAS{im0001.crop2.mha.md5}
        ${TEMP}/itkTubeNJetImageFunctionTest${testNum}.mha )
endforeach( testNum )

MIDAS_ADD_TEST( itktubeVotingResampleImageFunctionTest0
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest0.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest0.png
    itktubeVotingResampleImageFunctionTest
      0
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest0.png )

MIDAS_ADD_TEST( itktubeVotingResampleImageFunctionTest1
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest1.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest1.png
    itktubeVotingResampleImageFunctionTest
      1
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest1.png )

MIDAS_ADD_TEST( itktubeVotingResampleImageFunctionTest2
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest2.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest2.png
    itktubeVotingResampleImageFunctionTest
      2
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest2.png )

MIDAS_ADD_TEST( tubeSplineApproximation1DTest
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{tubeSplineApprox1DTest.mha.md5}
      ${TEMP}/tubeSplineApprox1DTest.mha
    tubeSplineApproximation1DTest
      ${TEMP}/tubeSplineApprox1DTest.mha )

MIDAS_ADD_TEST( tubeSplineNDTest
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{tubeSplineNDTest.mha.md5}
      ${TEMP}/tubeSplineNDTest.mha
    tubeSplineNDTest
      ${TEMP}/tubeSplineNDTest.mha )

add_test( itktubeRecordOptimizationParameterProgressionCommandTest
  ${BASE_NUMERICS_TESTS}
    itktubeRecordOptimizationParameterProgressionCommandTest
      ${TEMP}/itkRecordOptimizationParameterProgressionCommandTest.h5 )

MIDAS_ADD_TEST( itktubeBlurImageFunctionTest
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkTubeBlurImageFunctionTest.mha.md5}
      ${TEMP}/itkTubeBlurImageFunctionTest.mha
    itktubeBlurImageFunctionTest
      ${TEMP}/itkTubeBlurImageFunctionTest.mha )

MIDAS_ADD_TEST( itktubeNJetFeatureVectorGeneratorTest
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkTubeNJetFeatureVectorGeneratorTest_f0.mha.md5}
      ${TEMP}/itkTubeNJetFeatureVectorGeneratorTest_f0.mha
    itktubeNJetFeatureVectorGeneratorTest
      MIDAS{ES0015_Large_Subs.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
      ${TEMP}/itkTubeNJetFeatureVectorGeneratorTest_f0.mha
      ${TEMP}/itkTubeNJetFeatureVectorGeneratorTest_f1.mha )

MIDAS_ADD_TEST( itktubeNJetBasisFeatureVectorGeneratorTest
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkTubeNJetBasisFeatureVectorGeneratorTest_lda0.mha.md5}
      ${TEMP}/itkTubeNJetBasisFeatureVectorGeneratorTest_lda0.mha
    itktubeNJetBasisFeatureVectorGeneratorTest
      MIDAS{ES0015_Large_Subs.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
      255 0
      ${TEMP}/itkTubeNJetBasisFeatureVectorGeneratorTest_lda0.mha
      ${TEMP}/itkTubeNJetBasisFeatureVectorGeneratorTest_lda1.mha )

MIDAS_ADD_TEST( itktubeRidgeFeatureVectorGeneratorTest
  ${BASE_NUMERICS_TESTS}
    itktubeRidgeFeatureVectorGeneratorTest
      MIDAS{im0001.crop.contrast.mha.md5}
      ${TEMP}/itkTubeRidgeFeatureVectorGeneratorTest_f0.mha
      ${TEMP}/itkTubeRidgeFeatureVectorGeneratorTest_f1.mha )

MIDAS_ADD_TEST( itktubeRidgeBasisFeatureVectorGeneratorTest
  ${BASE_NUMERICS_TESTS}
    itktubeRidgeBasisFeatureVectorGeneratorTest
      MIDAS{im0001.crop.contrast.mha.md5}
      MIDAS{im0001.vk.mask.crop.mha.md5}
      255 127
      ${TEMP}/itkTubeRidgeBasisFeatureVectorGeneratorTest_lda0.mha
      ${TEMP}/itkTubeRidgeBasisFeatureVectorGeneratorTest_lda1.mha )
