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
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${EXECUTABLE_NAME}> )

set( IMAGECOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:ImageCompareCommand> )

set( TEXTCOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:TextCompareCommand> )

#######
# Test1
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --tubeId 255
    --backgroundId 127
    --unknownId 0
    --tubeScales 0.15,0.3
    --saveDiscriminantInfo ${TEMP}/${MODULE_NAME}Test1.mrs
    --labelmap MIDAS{im0001.vk.mask.crop.mha.md5}
    MIDAS{im0001.crop.contrast.mha.md5}
    ${TEMP}/${MODULE_NAME}Test1_Output.mha )

#######
# Test1-Compare MRS
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-MRS-Compare
  COMMAND ${TEXTCOMPARE_EXE}
    -d 0.1
    -t ${TEMP}/${MODULE_NAME}Test1.mrs
    -b MIDAS{${MODULE_NAME}Test1.mrs.md5} )
set_property(TEST ${MODULE_NAME}-Test1-MRS-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test1-Compare Output
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Output-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.001
    -t ${TEMP}/${MODULE_NAME}Test1_Output.mha
    -b MIDAS{${MODULE_NAME}Test1_Output.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test1-Output-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test2
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --loadDiscriminantInfo ${TEMP}/${MODULE_NAME}Test1.mrs
    --saveDiscriminantInfo ${TEMP}/${MODULE_NAME}Test2.mrs
    MIDAS{im0001.crop.contrast.mha.md5}
    ${TEMP}/${MODULE_NAME}Test2_Output.mha )
set_property(TEST ${MODULE_NAME}-Test2
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test2-Compare Output
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-Output-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.05
    -t ${TEMP}/${MODULE_NAME}Test2_Output.mha
    -b MIDAS{${MODULE_NAME}Test1_Output.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test2-Output-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

#######
# Test2-Compare MRS
#######
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-MRS-Compare
  COMMAND ${TEXTCOMPARE_EXE}
    -d 0.1
    -t ${TEMP}/${MODULE_NAME}Test2.mrs
    -b MIDAS{${MODULE_NAME}Test2.mrs.md5} )
set_property(TEST ${MODULE_NAME}-Test2-MRS-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

