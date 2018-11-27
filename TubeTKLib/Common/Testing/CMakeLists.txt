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

include_regular_expression( "^.*$" )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseCommon ${PYTHON_EXECUTABLE}
    ${TubeTKLib_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTKLib_SOURCE_DIR}/Common )
endif( PYTHON_EXECUTABLE )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( BASE_COMMON_TESTS
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseCommonTests> )

set( BASE_COMMON_HEADER_TEST
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseCommonHeaderTest> )

set( tubeBaseCommonTests_SRCS
  tubeBaseCommonTests.cxx
  tubeBaseCommonPrintTest.cxx
  tubeMacroTest.cxx
  tubeMessageTest.cxx
  tubeObjectTest.cxx )

include_directories(
  ${TubeTKLib_SOURCE_DIR}/Common
  ${TubeTKLib_SOURCE_DIR}/Filtering )

add_executable( tubeBaseCommonHeaderTest tubeBaseCommonHeaderTest.cxx )
target_link_libraries( tubeBaseCommonHeaderTest ${ITK_LIBRARIES} )

add_executable( tubeBaseCommonTests
  ${tubeBaseCommonTests_SRCS} )
target_link_libraries( tubeBaseCommonTests
  ${ITK_LIBRARIES}
  TubeTKLibCommon )

add_test( NAME tubeBaseCommonHeaderTest COMMAND ${BASE_COMMON_HEADER_TEST} )

add_test( NAME tubeBaseCommonPrintTest
  COMMAND ${BASE_COMMON_TESTS}
    tubeBaseCommonPrintTest )

add_test( NAME tubeMacroTest
  COMMAND ${BASE_COMMON_TESTS}
    tubeMacroTest )

add_test( NAME tubeMessageTest
  COMMAND ${BASE_COMMON_TESTS}
    tubeMessageTest )

add_test( NAME tubeObjectTest
  COMMAND ${BASE_COMMON_TESTS}
    tubeObjectTest )
