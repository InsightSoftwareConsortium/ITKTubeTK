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

include(${TubeTK_SOURCE_DIR}/CMake/MIDAS3.cmake )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( PythonInterp )
if( PYTHON_EXECUTABLE )
add_test( PrintSelf-BaseCommon ${PYTHON_EXECUTABLE}
    ${TubeTK_SOURCE_DIR}/CMake/PrintSelfCheck.py
    ${TubeTK_SOURCE_DIR}/Base/Common
    )
endif( PYTHON_EXECUTABLE )

set( BASE_COMMON_TESTS
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseCommonTests )

set( BASE_COMMON_HEADER_TEST
  ${TubeTK_BINARY_DIR}/${Slicer4_INSTALL_BIN_DIR}/tubeBaseCommonHeaderTest )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( tubeBaseCommonTests_SRCS
  tubeBaseCommonPrintTest.cxx
  tubeMessageTest.cxx
  )

include_directories(
  ${TubeTK_SOURCE_DIR}/Base/Common
  )

add_executable( tubeBaseCommonHeaderTest tubeBaseCommonHeaderTest.cxx )
target_link_libraries( tubeBaseCommonHeaderTest ${ITK_LIBRARIES} )

add_executable( tubeBaseCommonTests
  tubeBaseCommonTests.cxx
  ${tubeBaseCommonTests_SRCS}
  )
target_link_libraries( tubeBaseCommonTests
  ${ITK_LIBRARIES}
  ${JsonCpp_LIBRARIES}
  )

add_test( tubeBaseCommonHeaderTest ${BASE_COMMON_HEADER_TEST} )

add_test( tubeBaseCommonPrintTest
  ${BASE_COMMON_TESTS}
  tubeBaseCommonPrintTest
  )

add_test( tubeMessageTest
  ${BASE_COMMON_TESTS}
  tubeMessageTest
  )

MIDAS_ADD_TEST( ImageCompareCommand-Test1
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -t MIDAS{ES0015_Large.mha.md5}
  -b MIDAS{ES0015_Large.mha.md5}
  )

MIDAS_ADD_TEST( ImageCompareCommand-Test2
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -n 1000
  -r 1
  -i 50
  -t MIDAS{im0001.crop.mha.md5}
  -b MIDAS{im0001_n10.crop.mha.md5}
  )

MIDAS_ADD_TEST( ImageCompareCommand-Test3
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -t MIDAS{im0001.mha.md5}
  -b MIDAS{im0001_n10.crop.mha.md5}
  )
set_tests_properties( ImageCompareCommand-Test3 PROPERTIES WILL_FAIL true )

midas_add_test( ImageCompareCommand-Test4
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -t MIDAS{im0001.crop.mha.md5}
  -B 2 MIDAS{im0001_n10.crop.mha.md5} MIDAS{im0001.crop.mha.md5}
  )

midas_add_test( ImageCompareCommand-Test5
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -t MIDAS{im0001.crop.mha.md5}
  -B 2 MIDAS{im0001_n10.crop.mha.md5} MIDAS{im0001_n20.crop.mha.md5}
  )
set_tests_properties( ImageCompareCommand-Test5 PROPERTIES WILL_FAIL true )

midas_add_test( ImageCompareCommand-Test6
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand
  -t MIDAS{im0001_n10.crop.mha.md5}
  -b ${TEMP}/filedoesnotexist.mha
  )
set_tests_properties( ImageCompareCommand-Test6 PROPERTIES WILL_FAIL true )

add_test( JsonCompareCommand-Test1
  ${TubeTK_BINARY_DIR}/bin/JsonCompareCommand
  -t ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestTest1.json
  -b ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestBaseline1.json
  )

add_test( JsonCompareCommand-Test2
  ${TubeTK_BINARY_DIR}/bin/JsonCompareCommand
  -t ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestTest2.json
  -b ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestBaseline1.json
  )
set_tests_properties( JsonCompareCommand-Test2 PROPERTIES WILL_FAIL true )

add_test( JsonCompareCommand-Test3
  ${TubeTK_BINARY_DIR}/bin/JsonCompareCommand
  -t ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestTest2.json
  -b ${CMAKE_CURRENT_SOURCE_DIR}/JsonCompareCommandTestBaseline2.json
  )
