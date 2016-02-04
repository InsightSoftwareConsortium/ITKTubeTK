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
  add_test( PrintSelf-BaseNumerics ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Numerics )
endif( PYTHON_EXECUTABLE )

set( BASE_NUMERICS_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseNumericsTests> )

set( BASE_NUMERICS_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseNumericsHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseNumerics_SRCS
  tubeBaseNumericsPrintTest.cxx
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
  tubeTubeMathTest.cxx
  tubeUserFunctionTest.cxx )

add_executable( tubeBaseNumericsHeaderTest
  tubeBaseNumericsHeaderTest.cxx )
target_link_libraries( tubeBaseNumericsHeaderTest
  TubeTKCommon TubeTKFiltering TubeTKNumerics ${ITK_LIBRARIES} )
if( TubeTK_USE_GPU_ARRAYFIRE )
  target_link_libraries( tubeBaseNumericsHeaderTest ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_GPU_ARRAYFIRE )

set( tubeBaseNumerics_External_Libraries )
if( TubeTK_USE_GPU_ARRAYFIRE )
  set( tubeBaseNumerics_External_Libraries ${ArrayFire_LIBRARIES} )
endif( TubeTK_USE_GPU_ARRAYFIRE )

SEMMacroBuildCLI(
  NAME tubeBaseNumericsTests
  ADDITIONAL_SRCS
    ${tubeBaseNumerics_SRCS}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES
    ${tubeBaseNumerics_External_Libraries}
    TubeCLI TubeTKCommon TubeTKFiltering TubeTKNumerics ${ITK_LIBRARIES}
    ITKIOMeta ITKIOHDF5 ITKOptimizers
  INSTALL_RUNTIME_DESTINATION bin
  EXECUTABLE_ONLY
  )


add_test( NAME tubeBaseNumericsHeaderTest
  COMMAND ${BASE_NUMERICS_HEADER_TEST} )

add_test( NAME tubeBaseNumericsPrintTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeBaseNumericsPrintTest )

add_test( NAME tubeMatrixMathTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeMatrixMathTest )

add_test( NAME tubeBrentOptimizer1DTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeBrentOptimizer1DTest )

# Minimization - initialized at middle of interval
add_test( NAME tubeGoldenMeanOptimizer1DTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 0.01 -0.5 -1 )

# impossible maximization - not unimodal distribution
add_test( NAME tubeGoldenMeanOptimizer1DTest1
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 0 350 100 -1 )
set_tests_properties( tubeGoldenMeanOptimizer1DTest1
  PROPERTIES WILL_FAIL true )

# minimization - initialized at xMax
add_test( NAME tubeGoldenMeanOptimizer1DTest2
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 350 1 350 100 -1 )

# minimization - initialized at xMin
add_test( NAME tubeGoldenMeanOptimizer1DTest3
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      1 1 -3 3.5 1 -3 -0.5 -1 )

# maximization - initialized at middle of interval
add_test( NAME tubeGoldenMeanOptimizer1DTest4
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 100 0 1 )

# maximization - initialized at xMax
add_test( NAME tubeGoldenMeanOptimizer1DTest5
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 300 0 1 )

# maximization - initialized at xMin
add_test( NAME tubeGoldenMeanOptimizer1DTest6
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeGoldenMeanOptimizer1DTest
      100 2 -300 300 0 -300 0 1 )

add_test( NAME tubeParabolicFitOptimizer1DTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeParabolicFitOptimizer1DTest )

add_test( NAME tubeBrentOptimizerNDTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeBrentOptimizerNDTest )

add_test( NAME tubeUserFunctionTest
  COMMAND ${BASE_NUMERICS_TESTS}
    tubeUserFunctionTest )

Midas3FunctionAddTest( NAME itktubeJointHistogramImageFunctionTest01
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest01.mha.md5}
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
    --compareNumberOfPixelsTolerance 7
    itktubeJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      0
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01StdDev.mha )

Midas3FunctionAddTest( NAME itktubeJointHistogramImageFunctionTest02
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest02.mha.md5}
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
    --compareNumberOfPixelsTolerance 7
    itktubeJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      1
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02StdDev.mha )

Midas3FunctionAddTest( NAME itktubeImageRegionMomentsCalculatorTest
  COMMAND ${BASE_NUMERICS_TESTS}
    itktubeImageRegionMomentsCalculatorTest
      MIDAS{scoring-test.png.md5} )

foreach( testNum RANGE 0 37 )
  Midas3FunctionAddTest( NAME itktubeNJetImageFunctionTest${testNum}
    COMMAND ${BASE_NUMERICS_TESTS}
      --compare MIDAS{itktubeNJetImageFunctionTest${testNum}.mha.md5}
        ${TEMP}/itktubeNJetImageFunctionTest${testNum}.mha
        --compareIntensityTolerance 0.000000000001
      itktubeNJetImageFunctionTest
        ${testNum}
        MIDAS{im0001.crop2.mha.md5}
        ${TEMP}/itktubeNJetImageFunctionTest${testNum}.mha )
endforeach( testNum )

Midas3FunctionAddTest( NAME itktubeVotingResampleImageFunctionTest0
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest0.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest0.png
    itktubeVotingResampleImageFunctionTest
      0
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest0.png )

Midas3FunctionAddTest( NAME itktubeVotingResampleImageFunctionTest1
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest1.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest1.png
    itktubeVotingResampleImageFunctionTest
      1
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest1.png )

Midas3FunctionAddTest( NAME itktubeVotingResampleImageFunctionTest2
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest2.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest2.png
    itktubeVotingResampleImageFunctionTest
      2
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest2.png )

Midas3FunctionAddTest( NAME tubeSplineApproximation1DTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{tubeSplineApprox1DTest.mha.md5}
      ${TEMP}/tubeSplineApprox1DTest.mha
    tubeSplineApproximation1DTest
      ${TEMP}/tubeSplineApprox1DTest.mha )

Midas3FunctionAddTest( NAME tubeSplineNDTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{tubeSplineNDTest.mha.md5}
      ${TEMP}/tubeSplineNDTest.mha
    tubeSplineNDTest
      ${TEMP}/tubeSplineNDTest.mha )

add_test( NAME itktubeRecordOptimizationParameterProgressionCommandTest
  COMMAND ${BASE_NUMERICS_TESTS}
    itktubeRecordOptimizationParameterProgressionCommandTest
      ${TEMP}/itkRecordOptimizationParameterProgressionCommandTest.h5 )

Midas3FunctionAddTest( NAME itktubeBlurImageFunctionTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itktubeBlurImageFunctionTest.mha.md5}
      ${TEMP}/itktubeBlurImageFunctionTest.mha
    itktubeBlurImageFunctionTest
      ${TEMP}/itktubeBlurImageFunctionTest.mha )

Midas3FunctionAddTest( NAME itktubeNJetFeatureVectorGeneratorTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itktubeNJetFeatureVectorGeneratorTest_f0.mha.md5}
      ${TEMP}/itktubeNJetFeatureVectorGeneratorTest_f0.mha
    itktubeNJetFeatureVectorGeneratorTest
      MIDAS{ES0015_Large_Subs.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
      ${TEMP}/itktubeNJetFeatureVectorGeneratorTest_f0.mha
      ${TEMP}/itktubeNJetFeatureVectorGeneratorTest_f1.mha )

Midas3FunctionAddTest( NAME itktubeNJetBasisFeatureVectorGeneratorTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itktubeNJetBasisFeatureVectorGeneratorTest_lda0.mha.md5}
      ${TEMP}/itktubeNJetBasisFeatureVectorGeneratorTest_lda0.mha
    itktubeNJetBasisFeatureVectorGeneratorTest
      MIDAS{ES0015_Large_Subs.mha.md5}
      MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
      255 0
      ${TEMP}/itktubeNJetBasisFeatureVectorGeneratorTest_lda0.mha
      ${TEMP}/itktubeNJetBasisFeatureVectorGeneratorTest_lda1.mha )

Midas3FunctionAddTest( NAME itktubeSingleValuedCostFunctionImageSourceTest
  COMMAND ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itktubeSingleValuedCostFunctionImageSourceTest.mha.md5}
      ${TEMP}/itktubeSingleValuedCostFunctionImageSourceTest.mha
    itktubeSingleValuedCostFunctionImageSourceTest
      ${TEMP}/itktubeSingleValuedCostFunctionImageSourceTest.mha )

Midas3FunctionAddTest( NAME itktubeRidgeFFTFeatureVectorGeneratorTest
  COMMAND ${BASE_NUMERICS_TESTS}
    itktubeRidgeFFTFeatureVectorGeneratorTest
      MIDAS{im0001.crop.contrast.mha.md5}
      ${TEMP}/itktubeRidgeFFTFeatureVectorGeneratorTest_f0.mha
      ${TEMP}/itktubeRidgeFFTFeatureVectorGeneratorTest_f1.mha )

Midas3FunctionAddTest( NAME itktubeRidgeBasisFeatureVectorGeneratorTest
  COMMAND ${BASE_NUMERICS_TESTS}
    itktubeRidgeBasisFeatureVectorGeneratorTest
      MIDAS{im0001.crop.contrast.mha.md5}
      MIDAS{im0001.vk.mask.crop.mha.md5}
      255 127
      ${TEMP}/itktubeRidgeBasisFeatureVectorGeneratorTest_lda0.mha
      ${TEMP}/itktubeRidgeBasisFeatureVectorGeneratorTest_lda1.mha )

add_test( NAME tubeTubeMathTest
  COMMAND ${BASE_NUMERICS_TESTS}
  tubeTubeMathTest )
