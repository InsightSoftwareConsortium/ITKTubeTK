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

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseCLI ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/CLI )
endif( PYTHON_EXECUTABLE )

set( BASE_CLI_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseCLITests )

set( BASE_CLI_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseCLIHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

find_package( GenerateCLP REQUIRED )
include( ${GenerateCLP_USE_FILE} )

set( tubeCLIHelperFunctionsTest_SOURCE tubeCLIHelperFunctionsTest.cxx )
generateclp( tubeCLIHelperFunctionsTest_SOURCE tubeCLIHelperFunctionsTest.xml )

set( tubeBaseCLITests_SRCS
  tubeBaseCLIPrintTest.cxx
  tubeCLIFilterWatcherTest.cxx
  tubeCLIProgressReporterTest.cxx
  ${tubeCLIHelperFunctionsTest_SOURCE} )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/CLI
  ${TubeTK_SOURCE_DIR}/Base/Common )

add_executable( tubeBaseCLIHeaderTest tubeBaseCLIHeaderTest.cxx )
target_link_libraries( tubeBaseCLIHeaderTest ${ITK_LIBRARIES} )

add_executable( tubeBaseCLITests
  tubeBaseCLITests.cxx
  ${tubeBaseCLITests_SRCS} )
target_link_libraries( tubeBaseCLITests ${ITK_LIBRARIES} ${JsonCpp_LIBRARIES} )

add_test( tubeBaseCLIHeaderTest ${BASE_CLI_HEADER_TEST} )

add_test( tubeBaseCLIPrintTest
  ${BASE_CLI_TESTS}
  tubeBaseCLIPrintTest )

Midas3FunctionAddTest( tubeCLIFilterWatcherTest
  ${BASE_CLI_TESTS}
  tubeCLIFilterWatcherTest MIDAS{Branch.n010.mha.md5} )

Midas3FunctionAddTest( tubeCLIHelperFunctionsTest
  ${BASE_CLI_TESTS}
  tubeCLIHelperFunctionsTest MIDAS{Branch.n010.mha.md5} )

add_test( tubeCLIProgressReporterTest
  ${BASE_CLI_TESTS}
  tubeCLIProgressReporterTest )
