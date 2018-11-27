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
  add_test( PrintSelf-BaseObjectDocuments ${PYTHON_EXECUTABLE}
    ${TubeTKLib_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTKLib_SOURCE_DIR}/Base/ObjectDocuments )
endif( PYTHON_EXECUTABLE )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( BASE_OBJECT_DOCUMENTS_TESTS
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseObjectDocumentsTests> )

set( BASE_OBJECT_DOCUMENTS_HEADER_TEST
  ${TubeTKLib_LAUNCHER} $<TARGET_FILE:tubeBaseObjectDocumentsHeaderTest> )

set( TEMP ${TubeTKLib_BINARY_DIR}/Temporary )

set( tubeBaseObjectDocumentsTests_SRCS
  tubeBaseObjectDocumentsTests.cxx
  tubeBaseObjectDocumentsPrintTest.cxx )

include_directories(
  ${TubeTKLib_SOURCE_DIR}/Common
  ${TubeTKLib_SOURCE_DIR}/Filtering
  ${TubeTKLib_SOURCE_DIR}/ObjectDocuments )

add_executable( tubeBaseObjectDocumentsHeaderTest
  tubeBaseObjectDocumentsHeaderTest.cxx )
target_link_libraries( tubeBaseObjectDocumentsHeaderTest
  ${ITK_LIBRARIES} )

add_executable( tubeBaseObjectDocumentsTests 
  ${tubeBaseObjectDocumentsTests_SRCS} )
target_link_libraries( tubeBaseObjectDocumentsTests 
  ${ITK_LIBRARIES}
  TubeTKLibObjectDocuments )

add_test( NAME tubeBaseObjectDocumentsHeaderTest
  COMMAND ${BASE_OBJECT_DOCUMENTS_HEADER_TEST} )

add_test( NAME tubeBaseObjectDocumentsPrintTest
  COMMAND ${BASE_OBJECT_DOCUMENTS_TESTS}
  tubeBaseObjectDocumentsPrintTest )
