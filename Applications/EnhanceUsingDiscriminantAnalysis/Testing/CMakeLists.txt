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

set( CompareTextFiles_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

# Test1
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --useNumberOfPCABasis 0
    --objectId 0,255
    --saveBasisInfo ${TEMP}/${MODULE_NAME}Test1.mlda
    --labelmap DATA{${TubeTK_DATA_ROOT}/GDS0015_Large.mha}
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur1.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur2.mha}
    ${TEMP}/${MODULE_NAME}Test1 )

# Test1-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${CompareImages_EXE}
    -i 0.02
    -t ${TEMP}/${MODULE_NAME}Test1.basis00.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.basis00.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

# Test2
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 2
    --useNumberOfPCABasis 1
    --objectId 0,255
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur1.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur2.mha}
    --labelmap DATA{${TubeTK_DATA_ROOT}/GDS0015_Large.mha}
    ${TEMP}/${MODULE_NAME}Test2 )

# Test2-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${CompareImages_EXE}
    -i 0.02
    -t ${TEMP}/${MODULE_NAME}Test2.basis01.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.basis01.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

# Test3
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --objectId 0,255
    --loadBasisInfo ${TEMP}/${MODULE_NAME}Test1.mlda
    --saveBasisInfo ${TEMP}/${MODULE_NAME}Test3.mlda
    DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur1.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large_Blur2.mha}
    --labelmap DATA{${TubeTK_DATA_ROOT}/GDS0015_Large.mha}
    ${TEMP}/${MODULE_NAME}Test3 )
set_tests_properties( ${MODULE_NAME}-Test3 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

# Test3-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${CompareImages_EXE}
    -i 0.02
    -t ${TEMP}/${MODULE_NAME}Test3.basis00.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.basis00.mha} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test3 )

# Test3-Compare2
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test3-Compare2
  COMMAND ${CompareTextFiles_EXE}
    -d 0.01
    -t ${TEMP}/${MODULE_NAME}Test3.mlda
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3.mlda}
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test3-2.mlda} )
set_tests_properties( ${MODULE_NAME}-Test3-Compare2 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test3 )
