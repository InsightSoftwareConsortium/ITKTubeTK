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
  itkFeatureGeneratingImageFunctionTest.cxx 
  itkStandardFeatureGeneratingImageFunctionTest.cxx 
  )

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Numerics 
  )

add_executable( tubeBaseNumericsHeaderTest 
  tubeBaseNumericsHeaderTest.cxx 
  )
target_link_libraries( tubeBaseNumericsHeaderTest
  ITKIO 
  ITKNumerics 
  ITKBasicFilters 
  )

add_executable( tubeBaseNumericsTests
  tubeBaseNumericsTests.cxx
  ${tubeBaseNumerics_SRCS} 
  )
target_link_libraries( tubeBaseNumericsTests
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
      ${INPUT_DATA}/ES0015_1.mha
      ${INPUT_DATA}/GDS0015_1_match.mha
      ${TEMP}/itkJointHistogramImageFunctionTest01.mha 
      ${TEMP}/itkJointHistogramImageFunctionTest01Mean.mha 
      ${TEMP}/itkJointHistogramImageFunctionTest01StdDev.mha 
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

