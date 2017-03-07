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

set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

set( CompareImages_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

# Test1
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/Branch-truth.tre}
    ${TEMP}/${MODULE_NAME}Test1-Den.mha
    ${TEMP}/${MODULE_NAME}Test1-Rad.mha
    ${TEMP}/${MODULE_NAME}Test1-Tan.mha
    --inputTemplateImage DATA{${TubeTK_DATA_ROOT}/Branch.n010.mha}
    --useSquareDistance )

# Test1-Compare
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-Compare-Den
  COMMAND ${CompareImages_EXE}
    -i 0.001
    -t ${TEMP}/${MODULE_NAME}Test1-Den.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Den.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Den PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

ExternalData_Add_Test( TubeTKData NAME ${MODULE_NAME}-Test1-Compare-Rad
  COMMAND ${CompareImages_EXE}
    -i 0.001
    -n 50
    -t ${TEMP}/${MODULE_NAME}Test1-Rad.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Rad.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Rad PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )

ExternalData_Add_Test( TubeTKData NAME ${MODULE_NAME}-Test1-Compare-Tan
  COMMAND ${CompareImages_EXE}
    -i 0.00001
    -n 50
    -t ${TEMP}/${MODULE_NAME}Test1-Tan.mha
    -b DATA{${TubeTK_DATA_ROOT}/${MODULE_NAME}Test1-Tan.mha} )
set_tests_properties( ${MODULE_NAME}-Test1-Compare-Tan PROPERTIES DEPENDS
  ${MODULE_NAME}-Test1 )
