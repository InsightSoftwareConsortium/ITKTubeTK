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

set( PROJ_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareImages_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

# Test1
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --divideBy 2,4,2
    -p ${TEMP}/${MODULE_NAME}Test1_Points.mha
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${TEMP}/${MODULE_NAME}Test1.mha )

# Test1-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test1.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare2
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test1_Points.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1_Points.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare2 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

# Test2
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --newSize 88,64,4
    --overlap 1,1,10
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${TEMP}/${MODULE_NAME}Test2.mha )

# Test2-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test2.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

# Test3
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    --divideBy 2,4,2
    -m
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${TEMP}/${MODULE_NAME}Test3.mha )

# Test3-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test3.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3.mha} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test3 )

# Test4
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    --divideBy 2,4,2
    -g -l
    DATA{${TubeTK_DATA_ROOT}/T1NoSkull.mha}
    ${TEMP}/${MODULE_NAME}Test4.mha )

# Test4-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test4-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test4.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test4.mha} )
set_tests_properties( ${MODULE_NAME}-Test4-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test4 )
