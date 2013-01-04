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

Include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake)
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys)

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-UtilitiesFilters ${PYTHON_EXECUTABLE}
              ${TubeTK_SOURCE_DIR}/Utilities/PrintSelfCheck.py
              ${TubeTK_SOURCE_DIR}/Utilities )
endif()

if( Slicer3_INSTALL_BIN_DIR )
  set( UTILITIES_TESTS
        ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeUtilitiesTests )
else()
  set( UTILITIES_TESTS
        ${TubeTK_BINARY_DIR}/bin/tubeUtilitiesTests )
endif()

if( Slicer3_INSTALL_BIN_DIR )
  set( UTILITIES_HEADER_TEST
        ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_BIN_DIR}/tubeUtilitiesHeaderTest )
else()
  set( UTILITIES_HEADER_TEST
        ${TubeTK_BINARY_DIR}/bin/tubeUtilitiesHeaderTest )
endif()

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

# Test CLI generation
set( tubeCLIHelperFunctionsTest_SOURCE
     tubeCLIHelperFunctionsTest.cxx )
generateclp( tubeCLIHelperFunctionsTest_SOURCE
             tubeCLIHelperFunctionsTest.xml )

# Add test here
set( tubeUtilities_SRCS
  tubeUtilitiesPrintTest.cxx
  tubeMessageTest.cxx
  tubeCLIFilterWatcherTest.cxx
  tubeCLIProgressReporterTest.cxx
  ${tubeCLIHelperFunctionsTest_SOURCE} )

include_directories( ${TubeTK_SOURCE_DIR}/Utilities )

add_executable( tubeUtilitiesHeaderTest
  tubeUtilitiesHeaderTest.cxx )
target_link_libraries( tubeUtilitiesHeaderTest
  ${ITK_LIBRARIES} )

add_executable( tubeUtilitiesTests
  tubeUtilitiesTests.cxx
  ${tubeUtilities_SRCS} )
target_link_libraries( tubeUtilitiesTests
  ${ITK_LIBRARIES} )

add_test( tubeUtilitiesHeaderTest
  ${UTILITIES_HEADER_TEST} )

add_test( tubeUtilitiesPrintTest
            ${UTILITIES_TESTS}
            tubeUtilitiesPrintTest )

## Add tests here
add_test( tubeMessageTest
            ${UTILITIES_TESTS}
            tubeMessageTest )

MIDAS_ADD_TEST( tubeCLIFilterWatcherTest
            ${UTILITIES_TESTS}
            tubeCLIFilterWatcherTest MIDAS{Branch.n010.mha.md5} )

MIDAS_ADD_TEST( tubeCLIHelperFunctionsTest
            ${UTILITIES_TESTS}
            tubeCLIHelperFunctionsTest MIDAS{Branch.n010.mha.md5} )

add_test( tubeCLIProgressReporterTest
            ${UTILITIES_TESTS}
            tubeCLIProgressReporterTest )

## Image Compare
MIDAS_ADD_TEST( ImageCompareCommand-Test1
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
    -t MIDAS{ES0015_Large.mha.md5}
    -b MIDAS{ES0015_Large.mha.md5} )

MIDAS_ADD_TEST( ImageCompareCommand-Test2
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
    -n 1000
    -r 1
    -i 50
    -t MIDAS{im0001.crop.mha.md5}
    -b MIDAS{im0001_n10.crop.mha.md5} )

MIDAS_ADD_TEST( ImageCompareCommand-Test3
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
     -t MIDAS{im0001.mha.md5}
     -b MIDAS{im0001_n10.crop.mha.md5} )
set_tests_properties( ImageCompareCommand-Test3
                      PROPERTIES WILL_FAIL true )

MIDAS_ADD_TEST( ImageCompareCommand-Test4
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
     -t MIDAS{im0001.crop.mha.md5}
     -B 2 MIDAS{im0001_n10.crop.mha.md5} MIDAS{im0001.crop.mha.md5} )

MIDAS_ADD_TEST( ImageCompareCommand-Test5
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
     -t MIDAS{im0001.crop.mha.md5}
     -B 2 MIDAS{im0001_n10.crop.mha.md5} MIDAS{im0001_n20.crop.mha.md5} )
set_tests_properties( ImageCompareCommand-Test5
                      PROPERTIES WILL_FAIL true )

MIDAS_ADD_TEST( ImageCompareCommand-Test6
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
     -t MIDAS{im0001_n10.crop.mha.md5}
     -b ${TEMP}/filedoesnotexist.mha )
set_tests_properties( ImageCompareCommand-Test6
                      PROPERTIES WILL_FAIL true )

# Cppcheck static code analysis test
set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${TubeTK_SOURCE_DIR}/Utilities/Cppcheck/FindCppcheck.cmake")
option (TubeTK_USE_CPPCHECK
	"Create a test that does static code analysis with cppcheck." ON )
if( TubeTK_USE_CPPCHECK )
  find_package( Cppcheck REQUIRED )
  set( cppcheck_paths_to_check
	"${TubeTK_SOURCE_DIR}/Applications"
	"${TubeTK_SOURCE_DIR}/Base"
	"${TubeTK_SOURCE_DIR}/Utilities/Testing"
	)
  set (cppcheck_suppressions_file "${TubeTK_SOURCE_DIR}/Utilities/Cppcheck/TubeTKCppcheck.supp" )
  set (cppcheck_cmd_prefix "${CPPCHECK_EXECUTABLE}" -q --xml --suppressions "${cppcheck_suppressions_file}" )
  add_test( CppcheckStyle ${cppcheck_cmd_prefix} --enable=style ${cppcheck_paths_to_check} )
  add_test( CppcheckPerformance ${cppcheck_cmd_prefix} --enable=performance ${cppcheck_paths_to_check} )
  add_test( CppcheckPortability ${cppcheck_cmd_prefix} --enable=portability ${cppcheck_paths_to_check} )
  add_test( CppcheckInformation ${cppcheck_cmd_prefix} --enable=information ${cppcheck_paths_to_check} )
  add_test( CppcheckUnusedFunction ${cppcheck_cmd_prefix} --enable=unusedFunction ${cppcheck_paths_to_check} )
  set_tests_properties( CppcheckStyle CppcheckPerformance CppcheckPortability
    CppcheckInformation CppcheckUnusedFunction
    PROPERTIES FAIL_REGULAR_EXPRESSION "error file="
    )
endif( TubeTK_USE_CPPCHECK )

