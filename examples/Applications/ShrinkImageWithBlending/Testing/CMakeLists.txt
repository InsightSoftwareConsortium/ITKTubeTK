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
#       https://www.apache.org/licenses/LICENSE-2.0
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
    --divideBy 2,4,2
    -p ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1_Points.mha
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha )

# Test1-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare2
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test1_Points.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1_Points.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare2 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

# Test2
itk_add_test(
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --newSize 88,64,4
    --overlap 1,1,10
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha )

# Test2-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test2.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

# Test3
itk_add_test(
  NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    --divideBy 2,4,2
    -m
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha )

# Test3-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test3.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3.mha} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test3 )

# Test4
itk_add_test(
  NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    --divideBy 2,4,2
    -g -l
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha )

# Test4-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test4-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test4.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )

# Test5
itk_add_test(
  NAME ${MODULE_NAME}-Test5
  COMMAND ${PROJ_EXE}
    -i DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test5_Points.mha}
    -n 3,3
    DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test5_In.mha}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test5_Out.mha )

# Test5-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test5-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}Test5_Out.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test5_Out.mha} )
set_tests_properties( ${MODULE_NAME}-Test5-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test5 )
