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

#######
# Test1
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --tubeId 255
    --backgroundId 127
    --unknownId 0
    --tubeScales 0.18,0.36
    --saveDiscriminantInfo ${TEMP}/${MODULE_NAME}Test1.mrs
    --labelmap DATA{${TubeTK_DATA_ROOT}/im0001.vk.mask.crop.mha}
    DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
    ${TEMP}/${MODULE_NAME}Test1_Output.mha )

#######
# Test1-Compare MRS
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-MRS-Compare
  COMMAND ${TubeTKLib_CompareTextFiles_EXE}
    -d 0.1
    -t ${TEMP}/${MODULE_NAME}Test1.mrs
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1.mrs} )
set_tests_properties( ${MODULE_NAME}-Test1-MRS-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

#######
# Test1-Compare Output
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Output-Compare
  COMMAND ${TubeTKLib_CompareImages_EXE}
    -i 0.001
    -t ${TEMP}/${MODULE_NAME}Test1_Output.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1_Output.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Output-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

#######
# Test2
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --loadDiscriminantInfo ${TEMP}/${MODULE_NAME}Test1.mrs
    --saveDiscriminantInfo ${TEMP}/${MODULE_NAME}Test2.mrs
    DATA{${TubeTK_DATA_ROOT}/im0001.crop.contrast.mha}
    ${TEMP}/${MODULE_NAME}Test2_Output.mha )
set_tests_properties( ${MODULE_NAME}-Test2 PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1-Output-Compare )

#######
# Test2-Compare Output
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-Output-Compare
  COMMAND ${TubeTKLib_CompareImages_EXE}
    -i 0.05
    -n 100
    -t ${TEMP}/${MODULE_NAME}Test2_Output.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1_Output.mha} )
set_tests_properties( ${MODULE_NAME}-Test2-Output-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )

#######
# Test2-Compare MRS
#######
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-MRS-Compare
  COMMAND ${TubeTKLib_CompareTextFiles_EXE}
    -d 0.1
    -t ${TEMP}/${MODULE_NAME}Test2.mrs
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test2.mrs} )
set_tests_properties( ${MODULE_NAME}-Test2-MRS-Compare PROPERTIES DEPENDS
  ${MODULE_NAME}-Test2 )
