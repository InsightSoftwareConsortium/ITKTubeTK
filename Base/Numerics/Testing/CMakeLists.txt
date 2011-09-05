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
  add_test( PrintSelf-BaseNumerics ${TCL_TCLSH}
    ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.tcl
    ${TubeTK_SOURCE_DIR}/Base/Numerics )
endif( TCL_TCLSH )

set( BASE_NUMERICS_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsTests
  )

set( BASE_NUMERICS_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsHeaderTest
  )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Numerics/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseNumerics_SRCS
  tubeBaseNumericsPrintTest.cxx
  itkImageRegionMomentsCalculatorTest.cxx
  itkJointHistogramImageFunctionTest.cxx
  itkNJetImageFunctionTest.cxx
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
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_executable( tubeBaseNumericsTests
  tubeBaseNumericsTests.cxx
  ${tubeBaseNumerics_SRCS}
  )
target_link_libraries( tubeBaseNumericsTests
  TubeNumerics
  ITKIO
  ITKNumerics
  ITKBasicFilters
  )

add_test( tubeBaseNumericsHeaderTest
  ${BASE_NUMERICS_HEADER_TEST}
  )

add_test( tubeBaseNumericsPrintTest
  ${BASE_NUMERICS_TESTS}
  tubeBaseNumericsPrintTest
  )

add_test( itkJointHistogramImageFunctionTest01
  ${BASE_NUMERICS_TESTS}
    --compare ${BASELINE}/itkJointHistogramImageFunctionTest01.mha
              ${TEMP}/itkJointHistogramImageFunctionTest01.mha
    itkJointHistogramImageFunctionTest
      ${INPUT_DATA}/GDS0015_1_match_Subs.mha
      ${INPUT_DATA}/ES0015_1_Subs.mha
      ${INPUT_DATA}/GDS0015_1_match_Subs.mask.mha
      0
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01StdDev.mha
  )

add_test( itkJointHistogramImageFunctionTest02
  ${BASE_NUMERICS_TESTS}
    --compare ${BASELINE}/itkJointHistogramImageFunctionTest02.mha
              ${TEMP}/itkJointHistogramImageFunctionTest02.mha
    itkJointHistogramImageFunctionTest
      ${INPUT_DATA}/GDS0015_1_match_Subs.mha
      ${INPUT_DATA}/ES0015_1_Subs.mha
      ${INPUT_DATA}/GDS0015_1_match_Subs.mask.mha
      1
      ${TEMP}/itkJointHistogramImageFunctionTest02.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02Mean.mha
      ${TEMP}/itkJointHistogramImageFunctionTest02StdDev.mha
  )

add_test( itkImageRegionMomentsCalculatorTest
  ${BASE_NUMERICS_TESTS}
    itkImageRegionMomentsCalculatorTest
      ${INPUT_DATA}/scoring-test.png
  )

foreach( testNum RANGE 0 37 )
  add_test( itkNJetImageFunctionTest${testNum}
    ${BASE_NUMERICS_TESTS}
      --compare ${BASELINE}/itkNJetImageFunctionTest${testNum}.mha
                ${TEMP}/itkNJetImageFunctionTest${testNum}.mha
      itkNJetImageFunctionTest
        ${testNum}
        ${INPUT_DATA}/GDS0015_1.mha
        ${TEMP}/itkNJetImageFunctionTest${testNum}.mha
    )
endforeach( testNum )

add_test( itkVotingResampleImageFunctionTest0
  ${BASE_NUMERICS_TESTS}
    --compare ${BASELINE}/itkVotingResampleImageFunctionTest0.png
              ${TEMP}/itkVotingResampleImageFunctionTest0.png
    itkVotingResampleImageFunctionTest
      0
      ${INPUT_DATA}/greyscale01.png
      ${TEMP}/itkVotingResampleImageFunctionTest0.png
  )

add_test( itkVotingResampleImageFunctionTest1
  ${BASE_NUMERICS_TESTS}
    --compare ${BASELINE}/itkVotingResampleImageFunctionTest1.png
              ${TEMP}/itkVotingResampleImageFunctionTest1.png
    itkVotingResampleImageFunctionTest
      1
      ${INPUT_DATA}/greyscale01.png
      ${TEMP}/itkVotingResampleImageFunctionTest1.png
  )

add_test( itkVotingResampleImageFunctionTest2
  ${BASE_NUMERICS_TESTS}
    --compare ${BASELINE}/itkVotingResampleImageFunctionTest2.png
              ${TEMP}/itkVotingResampleImageFunctionTest2.png
    itkVotingResampleImageFunctionTest
      2
      ${INPUT_DATA}/greyscale01.png
      ${TEMP}/itkVotingResampleImageFunctionTest2.png
  )

add_test( tubeMatrixMathTest
  ${BASE_NUMERICS_TESTS}
  tubeMatrixMathTest
  )

add_test( tubeOptBrentTest
  ${BASE_NUMERICS_TESTS}
  tubeOptBrentTest
  )

add_test( tubeOptGoldenMeanTest
  ${BASE_NUMERICS_TESTS}
  tubeOptGoldenMeanTest
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

add_test( tubeSplineApprox1DTest
  ${BASE_NUMERICS_TESTS}
  --compare ${BASELINE}/tubeSplineApprox1DTest.mha
            ${TEMP}/tubeSplineApprox1DTest.mha
  tubeSplineApprox1DTest
    ${TEMP}/tubeSplineApprox1DTest.mha
  )

add_test( tubeSplineNDTest
  ${BASE_NUMERICS_TESTS}
  --compare ${BASELINE}/tubeSplineNDTest.mha
            ${TEMP}/tubeSplineNDTest.mha
  tubeSplineNDTest
    ${TEMP}/tubeSplineNDTest.mha
  )
