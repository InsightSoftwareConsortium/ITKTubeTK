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
       itkJointHistogramImageFunctionTest.cxx 
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

add_test( itkNJetImageFunctionTest18
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest18.mha
                      ${TEMP}/itkNJetImageFunctionTest18.mha
            itkNJetImageFunctionTest
               18
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest18.mha )

add_test( itkNJetImageFunctionTest19
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest19.mha
                      ${TEMP}/itkNJetImageFunctionTest19.mha
            itkNJetImageFunctionTest
               19
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest19.mha )

add_test( itkNJetImageFunctionTest20
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest20.mha
                      ${TEMP}/itkNJetImageFunctionTest20.mha
            itkNJetImageFunctionTest
               20
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest20.mha )

add_test( itkNJetImageFunctionTest21
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest21.mha
                      ${TEMP}/itkNJetImageFunctionTest21.mha
            itkNJetImageFunctionTest
               21
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest21.mha )

add_test( itkNJetImageFunctionTest22
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest22.mha
                      ${TEMP}/itkNJetImageFunctionTest22.mha
            itkNJetImageFunctionTest
               22
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest22.mha )

add_test( itkNJetImageFunctionTest23
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest23.mha
                      ${TEMP}/itkNJetImageFunctionTest23.mha
            itkNJetImageFunctionTest
               23
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest23.mha )

add_test( itkNJetImageFunctionTest24
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest24.mha
                      ${TEMP}/itkNJetImageFunctionTest24.mha
            itkNJetImageFunctionTest
               24
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest24.mha )

add_test( itkNJetImageFunctionTest25
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest25.mha
                      ${TEMP}/itkNJetImageFunctionTest25.mha
            itkNJetImageFunctionTest
               25
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest25.mha )

add_test( itkNJetImageFunctionTest26
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest26.mha
                      ${TEMP}/itkNJetImageFunctionTest26.mha
            itkNJetImageFunctionTest
               26
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest26.mha )

add_test( itkNJetImageFunctionTest27
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest27.mha
                      ${TEMP}/itkNJetImageFunctionTest27.mha
            itkNJetImageFunctionTest
               27
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest27.mha )

add_test( itkNJetImageFunctionTest28
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest28.mha
                      ${TEMP}/itkNJetImageFunctionTest28.mha
            itkNJetImageFunctionTest
               28
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest28.mha )

add_test( itkNJetImageFunctionTest29
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest29.mha
                      ${TEMP}/itkNJetImageFunctionTest29.mha
            itkNJetImageFunctionTest
               29
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest29.mha )

add_test( itkNJetImageFunctionTest30
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest30.mha
                      ${TEMP}/itkNJetImageFunctionTest30.mha
            itkNJetImageFunctionTest
               30
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest30.mha )

add_test( itkNJetImageFunctionTest31
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest31.mha
                      ${TEMP}/itkNJetImageFunctionTest31.mha
            itkNJetImageFunctionTest
               31
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest31.mha )

add_test( itkNJetImageFunctionTest32
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest32.mha
                      ${TEMP}/itkNJetImageFunctionTest32.mha
            itkNJetImageFunctionTest
               32
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest32.mha )

add_test( itkNJetImageFunctionTest33
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest33.mha
                      ${TEMP}/itkNJetImageFunctionTest33.mha
            itkNJetImageFunctionTest
               33
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest33.mha )

add_test( itkNJetImageFunctionTest34
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest34.mha
                      ${TEMP}/itkNJetImageFunctionTest34.mha
            itkNJetImageFunctionTest
               34
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest34.mha )

add_test( itkNJetImageFunctionTest35
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest35.mha
                      ${TEMP}/itkNJetImageFunctionTest35.mha
            itkNJetImageFunctionTest
               35
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest35.mha )

add_test( itkNJetImageFunctionTest36
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest36.mha
                      ${TEMP}/itkNJetImageFunctionTest36.mha
            itkNJetImageFunctionTest
               36
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest36.mha )

add_test( itkNJetImageFunctionTest37
            ${BASE_NUMERICS_TESTS}
            --compare ${BASELINE}/itkNJetImageFunctionTest37.mha
                      ${TEMP}/itkNJetImageFunctionTest37.mha
            itkNJetImageFunctionTest
               37
               ${INPUT_DATA}/GDS0015_1.mha
               ${TEMP}/itkNJetImageFunctionTest37.mha )

add_test( itkJointHistogramImageFunctionTest01
  ${BASE_NUMERICS_TESTS}
  --compare ${BASELINE}/itkJointHistogramImageFunctionTest01.mha
  ${TEMP}/itkJointHistogramImageFunctionTest01.mha
  itkJointHistogramImageFunctionTest
  ${INPUT_DATA}/ES0015_1.mha
  ${INPUT_DATA}/GDS0015_1.mha
  ${TEMP}/itkJointHistogramImageFunctionTest01.mha )