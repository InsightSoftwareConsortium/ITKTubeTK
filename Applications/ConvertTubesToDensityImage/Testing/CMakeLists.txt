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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

set( CompareImages_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

# Test1
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    MIDAS{Branch-truth.tre.md5}
    ${TEMP}/${MODULE_NAME}Test1-Den.mha
    ${TEMP}/${MODULE_NAME}Test1-Rad.mha
    ${TEMP}/${MODULE_NAME}Test1-Tan.mha
    --inputTemplateImage MIDAS{Branch.n010.mha.md5}
    --useSquareDistance )

# Test1-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare-Den
  COMMAND ${CompareImages_EXE}
    -i 0.001
    -t ${TEMP}/${MODULE_NAME}Test1-Den.mha
    -b MIDAS{${MODULE_NAME}Test1-Den.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare-Den
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare-Rad
  COMMAND ${CompareImages_EXE}
    -i 0.001
    -n 50
    -t ${TEMP}/${MODULE_NAME}Test1-Rad.mha
    -b MIDAS{${MODULE_NAME}Test1-Rad.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare-Rad
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare-Tan
  COMMAND ${CompareImages_EXE}
    -i 0.00001
    -n 50
    -t ${TEMP}/${MODULE_NAME}Test1-Tan.mha
    -b MIDAS{${MODULE_NAME}Test1-Tan.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare-Tan
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )
