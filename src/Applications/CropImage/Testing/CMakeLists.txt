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

set( PROJ_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1
itk_add_test(
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha
    -m 10,10 -M 100,100 )

itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

# Test2
itk_add_test(
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha
    -m 10,10 -s 91,91 )

itk_add_test(
  NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

# Test3
itk_add_test(
  NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha
    -m 15,15 -s 81,81 -b 5,5 )

itk_add_test(
  NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test3 )

# Test4
itk_add_test(
  NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4
    -S 2,2 -b 5,5 )

itk_add_test(
  NAME ${MODULE_NAME}-Test4-Compare000000
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4_000000.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4_000000.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare000000 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )

itk_add_test(
  NAME ${MODULE_NAME}-Test4-Compare000001
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4_000001.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4_000001.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare000001 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )

itk_add_test(
  NAME ${MODULE_NAME}-Test4-Compare001000
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4_001000.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4_001000.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare001000 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )

itk_add_test(
  NAME ${MODULE_NAME}-Test4-Compare001001
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4_001001.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4_001001.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare001001 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )
