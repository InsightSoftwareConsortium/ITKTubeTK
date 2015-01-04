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

set( IMAGECOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:ImageCompareCommand> )

set( TEXTCOMPARE_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:TextCompareCommand> )

# Test1
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --useNumberOfPCABasis 3
    --objectId 0,255
    --zeroScales 1
    --firstScales 3
    --secondScales 5
    --saveBasisInfo ${TEMP}/${MODULE_NAME}Test1.mnda
    --saveFeatureImages ${TEMP}/${MODULE_NAME}Test1
    MIDAS{ES0015_Large_Subs.mha.md5}
    --labelmap MIDAS{GDS0015_Large_Subs2.mha.md5}
    ${TEMP}/${MODULE_NAME}Test1 )

# Test1-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.01
    -t ${TEMP}/${MODULE_NAME}Test1.basis00.mha
    -b MIDAS{${MODULE_NAME}Test1.basis00.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

# Test2
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --useNumberOfPCABasis 0
    --objectId 0,255
    --zeroScales 1
    --firstScales 3
    --secondScales 5
    MIDAS{ES0015_Large_Subs.mha.md5}
    --labelmap MIDAS{GDS0015_Large_Subs2.mha.md5}
    ${TEMP}/${MODULE_NAME}Test2 )

# Test2-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.01
    -t ${TEMP}/${MODULE_NAME}Test2.basis01.mha
    -b MIDAS{${MODULE_NAME}Test2.basis01.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test2-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

# Test3
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --useNumberOfPCABasis 0
    --saveFeatureImages ${TEMP}/${MODULE_NAME}Test3
    --loadBasisInfo ${TEMP}/${MODULE_NAME}Test1.mnda
    --saveBasisInfo ${TEMP}/${MODULE_NAME}Test3.mnda
    MIDAS{ES0015_Large_Subs.mha.md5}
    ${TEMP}/${MODULE_NAME}Test3 )
set_property( TEST ${MODULE_NAME}-Test3
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1-Compare )

# Test3-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.02
    -t ${TEMP}/${MODULE_NAME}Test3.basis00.mha
    -b MIDAS{${MODULE_NAME}Test1.basis00.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test3-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )

# Test3-Compare2
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3-Compare2
  COMMAND ${TEXTCOMPARE_EXE}
    -d 0.01
    -t ${TEMP}/${MODULE_NAME}Test3.mnda
    -b MIDAS{${MODULE_NAME}Test1.mnda.md5} )
set_property( TEST ${MODULE_NAME}-Test3-Compare2
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )

# Test3-Compare3
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3-Compare3
  COMMAND ${IMAGECOMPARE_EXE}
    -i 0.01
    -t ${TEMP}/${MODULE_NAME}Test3.f01.mha
    -b ${TEMP}/${MODULE_NAME}Test1.f01.mha )
set_property( TEST ${MODULE_NAME}-Test3-Compare3
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )

# Test4
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    --useNumberOfLDABasis 1
    --useNumberOfPCABasis 3
    --objectId 127,255
    --zeroScales 1
    --firstScales 3,4
    --secondScales 5
    --loadBasisInfo ${TEMP}/${MODULE_NAME}Test1.mnda
    --saveBasisInfo ${TEMP}/${MODULE_NAME}Test4.mnda
    MIDAS{im0001_n05.crop.mha.md5}
    ${TEMP}/${MODULE_NAME}Test4 )
set_property( TEST ${MODULE_NAME}-Test4
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

# Test4-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4-Compare
  COMMAND ${TEXTCOMPARE_EXE}
    -d 0.01
    -t ${TEMP}/${MODULE_NAME}Test4.mnda
    -b MIDAS{${MODULE_NAME}Test1.mnda.md5} )
set_property( TEST ${MODULE_NAME}-Test4-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )
