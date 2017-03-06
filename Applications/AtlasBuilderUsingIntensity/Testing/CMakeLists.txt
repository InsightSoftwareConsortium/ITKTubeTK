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

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareImages_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

configure_file( ${TubeTK_SOURCE_DIR}/Applications/${MODULE_NAME}/Testing/ListTemplate.txt.in
                ${TEMP}/${MODULE_NAME}-List.txt IMMEDIATE @ONLY )

# Fetch data
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-FetchData
  COMMAND ${CMAKE_COMMAND} -E md5sum
    DATA{${TubeTK_DATA_ROOT}/DDMap58.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap86.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap57.mha}
    DATA{${TubeTK_DATA_ROOT}/DDMap25.mha} )

# Test1
add_test(	NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    ${TEMP}/${MODULE_NAME}-List.txt
    ${TEMP}/${MODULE_NAME}-Test1-Mean.mha
            ${TEMP}/${MODULE_NAME}-Test1-Variance.mha
            --outputSize 128,128,128
            --outputSpacing 2,2,2 )
set_tests_properties(${MODULE_NAME}-Test1 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1-FetchData)

# Test1-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare00
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test1-Mean.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test1-Mean.mha} )
set_tests_properties(${MODULE_NAME}-Test1-Compare00 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1)

ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare01
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test1-Variance.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}-Test1-Variance.mha} )
set_tests_properties(${MODULE_NAME}-Test1-Compare01 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1)
