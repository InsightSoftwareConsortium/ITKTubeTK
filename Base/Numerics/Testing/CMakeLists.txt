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
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Numerics )
endif()

if( Slicer3_INSTALL_BIN_DIR )
  set( BASE_NUMERICS_TESTS
        ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsTests )
else()
  set( BASE_NUMERICS_TESTS
        ${TubeTK_BINARY_DIR}/bin/tubeBaseNumericsTests )
endif()

if( Slicer3_INSTALL_BIN_DIR )
  set( BASE_NUMERICS_HEADER_TEST
        ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsHeaderTest )
else()
  set( BASE_NUMERICS_HEADER_TEST
        ${TubeTK_BINARY_DIR}/bin/tubeBaseNumericsHeaderTest )
endif()

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseNumerics_SRCS
  tubeBaseNumericsPrintTest.cxx
  itkImageRegionMomentsCalculatorTest.cxx
  itkJointHistogramImageFunctionTest.cxx
  itkTubeNJetImageFunctionTest.cxx
  itkTubeBlurImageFunctionTest.cxx
  itkTubeLDAGeneratorTest.cxx
  itkTubeNJetLDAGeneratorTest.cxx
  itkTubeNJetLDAGenerator2FTest.cxx
  itkVotingResampleImageFunctionTest.cxx
  tubeMatrixMathTest.cxx
  tubeOptBrentNDTest.cxx
  tubeOptBrentTest.cxx
  tubeOptGoldenMeanTest.cxx
  tubeOptParabolicFitTest.cxx
  tubeSplineApprox1DTest.cxx
  tubeSplineNDTest.cxx
  tubeUserFuncTest.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  )

add_executable( tubeBaseNumericsHeaderTest
  tubeBaseNumericsHeaderTest.cxx
  )
target_link_libraries( tubeBaseNumericsHeaderTest
  TubeNumerics
  ${ITK_LIBRARIES}
  )

add_executable( tubeBaseNumericsTests
  tubeBaseNumericsTests.cxx
  ${tubeBaseNumerics_SRCS}
  )
target_link_libraries( tubeBaseNumericsTests
  TubeNumerics
  ${ITK_LIBRARIES}
  )

add_test( tubeBaseNumericsHeaderTest
  ${BASE_NUMERICS_HEADER_TEST}
  )

add_test( tubeBaseNumericsPrintTest
  ${BASE_NUMERICS_TESTS}
  tubeBaseNumericsPrintTest
  )

add_test( tubeMatrixMathTest
  ${BASE_NUMERICS_TESTS}
  tubeMatrixMathTest
  )

add_test( tubeOptBrentTest
  ${BASE_NUMERICS_TESTS}
  tubeOptBrentTest
  )

# Minimization - initialized at middle of interval
add_test( tubeOptGoldenMeanTest
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  1 1 -3 3.5 1 0.01 -0.5 -1
  )

# impossible maximization - not unimodal distribution
add_test( tubeOptGoldenMeanTest1
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  100 2 -300 350 0 350 100 -1
  )
set_tests_properties( tubeOptGoldenMeanTest1
                      PROPERTIES WILL_FAIL true )

# minimization - initialized at xMax
add_test( tubeOptGoldenMeanTest2
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  100 2 -300 350 1 350 100 -1
  )

# minimization - initialized at xMin
add_test( tubeOptGoldenMeanTest3
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  1 1 -3 3.5 1 -3 -0.5 -1  )

# maximization - initialized at middle of interval
add_test( tubeOptGoldenMeanTest4
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  100 2 -300 300 0 100 0 1
  )

# maximization - initialized at xMax
add_test( tubeOptGoldenMeanTest5
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  100 2 -300 300 0 300 0 1
  )

# maximization - initialized at xMin
add_test( tubeOptGoldenMeanTest6
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
  100 2 -300 300 0 -300 0 1
  )

add_test( tubeOptParabolicFitTest
  ${BASE_NUMERICS_TESTS}
  tubeOptParabolicFitTest
  )

add_test( tubeOptBrentNDTest
  ${BASE_NUMERICS_TESTS}
  tubeOptBrentNDTest
  )

add_test( tubeUserFuncTest
  ${BASE_NUMERICS_TESTS}
  tubeUserFuncTest
  )

MIDAS_ADD_TEST( itkJointHistogramImageFunctionTest01
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest01.mha.md5}
              ${TEMP}/itkJointHistogramImageFunctionTest01.mha
    itkJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      0
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01StdDev.mha
  )

MIDAS_ADD_TEST( itkJointHistogramImageFunctionTest02
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkJointHistogramImageFunctionTest02.mha.md5}
              ${TEMP}/itkJointHistogramImageFunctionTest02.mha
    itkJointHistogramImageFunctionTest
      MIDAS{GDS0015_1_match_Subs.mha.md5}
      MIDAS{ES0015_1_Subs.mha.md5}
      MIDAS{GDS0015_1_match_Subs.mask.mha.md5}
      1
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02StdDev.mha
  )

MIDAS_ADD_TEST( itkImageRegionMomentsCalculatorTest
  ${BASE_NUMERICS_TESTS}
    itkImageRegionMomentsCalculatorTest
      MIDAS{scoring-test.png.md5}
  )

foreach( testNum RANGE 0 37 )
  MIDAS_ADD_TEST( itkTubeNJetImageFunctionTest${testNum}
    ${BASE_NUMERICS_TESTS}
      --compare MIDAS{itkTubeNJetImageFunctionTest${testNum}.mha.md5}
                ${TEMP}/itkTubeNJetImageFunctionTest${testNum}.mha
      itkTubeNJetImageFunctionTest
        ${testNum}
        MIDAS{GDS0015_1_Subs.mha.md5}
        ${TEMP}/itkTubeNJetImageFunctionTest${testNum}.mha
    )
endforeach( testNum )

MIDAS_ADD_TEST( itkVotingResampleImageFunctionTest0
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest0.png.md5}
              ${TEMP}/itkVotingResampleImageFunctionTest0.png
    itkVotingResampleImageFunctionTest
      0
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest0.png
  )

MIDAS_ADD_TEST( itkVotingResampleImageFunctionTest1
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest1.png.md5}
              ${TEMP}/itkVotingResampleImageFunctionTest1.png
    itkVotingResampleImageFunctionTest
      1
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest1.png
  )

MIDAS_ADD_TEST( itkVotingResampleImageFunctionTest2
  ${BASE_NUMERICS_TESTS}
    --compare MIDAS{itkVotingResampleImageFunctionTest2.png.md5}
              ${TEMP}/itkVotingResampleImageFunctionTest2.png
    itkVotingResampleImageFunctionTest
      2
      MIDAS{greyscale01.png.md5}
      ${TEMP}/itkVotingResampleImageFunctionTest2.png
  )

MIDAS_ADD_TEST( tubeSplineApprox1DTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{tubeSplineApprox1DTest.mha.md5}
            ${TEMP}/tubeSplineApprox1DTest.mha
  tubeSplineApprox1DTest
    ${TEMP}/tubeSplineApprox1DTest.mha
  )

MIDAS_ADD_TEST( tubeSplineNDTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{tubeSplineNDTest.mha.md5}
            ${TEMP}/tubeSplineNDTest.mha
  tubeSplineNDTest
    ${TEMP}/tubeSplineNDTest.mha
  )

MIDAS_ADD_TEST( itkTubeBlurImageFunctionTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{itkTubeBlurImageFunctionTest.mha.md5}
            ${TEMP}/itkTubeBlurImageFunctionTest.mha
  itkTubeBlurImageFunctionTest
    ${TEMP}/itkTubeBlurImageFunctionTest.mha
  )

MIDAS_ADD_TEST( itkTubeLDAGeneratorTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{itkTubeLDAGeneratorTest_lda0.mha.md5}
            ${TEMP}/itkTubeLDAGeneratorTest_lda0.mha
  itkTubeLDAGeneratorTest
    MIDAS{ES0015_Large_Subs.mha.md5}
    MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
    ${TEMP}/itkTubeLDAGeneratorTest_lda0.mha
    ${TEMP}/itkTubeLDAGeneratorTest_lda1.mha
  )

MIDAS_ADD_TEST( itkTubeNJetLDAGeneratorTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{itkTubeNJetLDAGeneratorTest_lda0.mha.md5}
            ${TEMP}/itkTubeNJetLDAGeneratorTest_lda0.mha
  itkTubeNJetLDAGeneratorTest
    MIDAS{ES0015_Large_Subs.mha.md5}
    MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
    ${TEMP}/itkTubeNJetLDAGeneratorTest_lda0.mha
    ${TEMP}/itkTubeNJetLDAGeneratorTest_lda1.mha
  )

MIDAS_ADD_TEST( itkTubeNJetLDAGenerator2FTest
  ${BASE_NUMERICS_TESTS}
  --compare MIDAS{itkTubeNJetLDAGeneratorTest_lda0.mha.md5}
            ${TEMP}/itkTubeNJetLDAGenerator2FTest_lda0.mha
  itkTubeNJetLDAGenerator2FTest
    MIDAS{ES0015_Large_Subs.mha.md5}
    MIDAS{GDS0015_Large-TrainingMask_Subs.mha.md5}
    ${TEMP}/itkTubeNJetLDAGenerator2FTest_lda0.mha
    ${TEMP}/itkTubeNJetLDAGenerator2FTest_lda1.mha
  )
