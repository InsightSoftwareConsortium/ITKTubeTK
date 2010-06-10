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
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsTests )

set( BASE_NUMERICS_HEADER_TEST 
      ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeBaseNumericsHeaderTest )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )
set( BASELINE ${TubeTK_SOURCE_DIR}/Base/Numerics/Testing/Baselines )
set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseNumerics_SRCS
       tubeBaseNumericsPrintTest.cxx
       itkNJetImageFunctionTest.cxx 
)

include_directories( 
  ${TubeTK_SOURCE_DIR}/Utilities 
  ${TubeTK_SOURCE_DIR}/Base/Numerics 
  )

add_executable( tubeBaseNumericsHeaderTest 
                  tubeBaseNumericsHeaderTest.cxx )
target_link_libraries( tubeBaseNumericsHeaderTest
                         ITKIO ITKNumerics ITKBasicFilters )

add_executable( tubeBaseNumericsTests
                  tubeBaseNumericsTests.cxx
                  ${tubeBaseNumerics_SRCS} )
target_link_libraries( tubeBaseNumericsTests
                         ITKIO ITKNumerics ITKBasicFilters )

add_test( tubeBaseNumericsHeaderTest 
            ${BASE_NUMERICS_HEADER_TEST} )

add_test( tubeBaseNumericsPrintTest 
            ${BASE_NUMERICS_TESTS}
            tubeBaseNumericsPrintTest )

add_test( itkNJetImageFunctionTest0
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest0.mha
                      ${TEMP}/itkNJetImageFunctionTest0.mha
            itkNJetImageFunctionTest
               0
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest0.mha )

add_test( itkNJetImageFunctionTest1
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest1.mha
                      ${TEMP}/itkNJetImageFunctionTest1.mha
            itkNJetImageFunctionTest
               1
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest1.mha )

add_test( itkNJetImageFunctionTest2
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest2.mha
                      ${TEMP}/itkNJetImageFunctionTest2.mha
            itkNJetImageFunctionTest
               2
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest2.mha )

add_test( itkNJetImageFunctionTest3
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest3.mha
                      ${TEMP}/itkNJetImageFunctionTest3.mha
            itkNJetImageFunctionTest
               3
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest3.mha )

add_test( itkNJetImageFunctionTest4
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest4.mha
                      ${TEMP}/itkNJetImageFunctionTest4.mha
            itkNJetImageFunctionTest
               4
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest4.mha )

add_test( itkNJetImageFunctionTest5
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest5.mha
                      ${TEMP}/itkNJetImageFunctionTest5.mha
            itkNJetImageFunctionTest
               5
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest5.mha )

add_test( itkNJetImageFunctionTest6
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest6.mha
                      ${TEMP}/itkNJetImageFunctionTest6.mha
            itkNJetImageFunctionTest
               6
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest6.mha )

add_test( itkNJetImageFunctionTest7
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest7.mha
                      ${TEMP}/itkNJetImageFunctionTest7.mha
            itkNJetImageFunctionTest
               7
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest7.mha )

add_test( itkNJetImageFunctionTest8
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest8.mha
                      ${TEMP}/itkNJetImageFunctionTest8.mha
            itkNJetImageFunctionTest
               8
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest8.mha )

add_test( itkNJetImageFunctionTest9
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest9.mha
                      ${TEMP}/itkNJetImageFunctionTest9.mha
            itkNJetImageFunctionTest
               9
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest9.mha )

add_test( itkNJetImageFunctionTest10
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest10.mha
                      ${TEMP}/itkNJetImageFunctionTest10.mha
            itkNJetImageFunctionTest
               10
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest10.mha )

add_test( itkNJetImageFunctionTest11
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest11.mha
                      ${TEMP}/itkNJetImageFunctionTest11.mha
            itkNJetImageFunctionTest
               11
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest11.mha )

add_test( itkNJetImageFunctionTest12
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest12.mha
                      ${TEMP}/itkNJetImageFunctionTest12.mha
            itkNJetImageFunctionTest
               12
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest12.mha )

add_test( itkNJetImageFunctionTest13
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest13.mha
                      ${TEMP}/itkNJetImageFunctionTest13.mha
            itkNJetImageFunctionTest
               13
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest13.mha )

add_test( itkNJetImageFunctionTest14
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest14.mha
                      ${TEMP}/itkNJetImageFunctionTest14.mha
            itkNJetImageFunctionTest
               14
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest14.mha )

add_test( itkNJetImageFunctionTest15
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest15.mha
                      ${TEMP}/itkNJetImageFunctionTest15.mha
            itkNJetImageFunctionTest
               15
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest15.mha )

add_test( itkNJetImageFunctionTest16
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest16.mha
                      ${TEMP}/itkNJetImageFunctionTest16.mha
            itkNJetImageFunctionTest
               16
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest16.mha )

add_test( itkNJetImageFunctionTest17
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest17.mha
                      ${TEMP}/itkNJetImageFunctionTest17.mha
            itkNJetImageFunctionTest
               17
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest17.mha )

