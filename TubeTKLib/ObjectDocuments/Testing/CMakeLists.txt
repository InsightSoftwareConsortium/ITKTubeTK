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

if( PYTHON_EXECUTABLE )
  add_test( PrintSelf-BaseObjectDocuments ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/TubeTKLib/ObjectDocuments )
endif( PYTHON_EXECUTABLE )

set( BASE_OBJECT_DOCUMENTS_TESTS
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseObjectDocumentsTests> )

set( BASE_OBJECT_DOCUMENTS_HEADER_TEST
  ${TubeTK_LAUNCHER} $<TARGET_FILE:tubeBaseObjectDocumentsHeaderTest> )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseObjectDocumentsTests_SRCS
  tubeBaseObjectDocumentsTests.cxx
  tubeBaseObjectDocumentsPrintTest.cxx )

include_directories(
  ${ITK_INCLUDE_DIRS}
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Common
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Filtering
  ${TubeTK_SOURCE_DIR}/TubeTKLib/ObjectDocuments )

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
