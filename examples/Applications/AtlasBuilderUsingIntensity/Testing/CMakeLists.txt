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

set( PROJ_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

configure_file(
  ${TubeTK_SOURCE_DIR}/src/Applications/${MODULE_NAME}/Testing/ListTemplate.txt.in
  ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-List.txt IMMEDIATE @ONLY
  )

# Fetch data
itk_add_test(
  NAME ${MODULE_NAME}-Test1-FetchData
  COMMAND ${CMAKE_COMMAND} -E copy
    DATA{${TubeTK_DATA_ROOT}/DDMap58.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap86.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap57.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap25.mha}
    ${ITK_TEST_OUTPUT_DIR} )

# Test1
itk_add_test(	NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-List.txt
    ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1-Mean.mha
            ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1-Variance.mha
            --outputSize 128,128,128
            --outputSpacing 2,2,2 )
set_tests_properties(${MODULE_NAME}-Test1 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1-FetchData)

# Test1-Compare
itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare00
  WORKING_DIRECTORY ${ITK_TEST_OUTPUT_DIR}
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1-Mean.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test1-Mean.mha} )
set_tests_properties(${MODULE_NAME}-Test1-Compare00 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1)

itk_add_test(
  NAME ${MODULE_NAME}-Test1-Compare01
  COMMAND ${TubeTK_CompareImages_EXE}
    CompareImages
    -t ${ITK_TEST_OUTPUT_DIR}/${MODULE_NAME}-Test1-Variance.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test1-Variance.mha} )
set_tests_properties(${MODULE_NAME}-Test1-Compare01 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1)
