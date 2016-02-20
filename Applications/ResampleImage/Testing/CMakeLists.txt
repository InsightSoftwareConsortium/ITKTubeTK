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

# Test1
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --interpolator Linear
    --resampleFactor 0.2,0.2
    MIDAS{ES0015_Large.mha.md5}
    ${TEMP}/${MODULE_NAME}Test1.mha )

# Test1-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}Test1.mha
    -b MIDAS{${MODULE_NAME}Test1.mha.md5},MIDAS{${MODULE_NAME}Test1_2.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

# Test2
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2
  COMMAND ${PROJ_EXE}
    --interpolator BSpline
    --makeIsotropic
    MIDAS{CroppedWholeLungCTScan.mhd.md5}
    ${TEMP}/${MODULE_NAME}Test2.mha
    MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )

# Test2-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    --toleranceCoordinates 1.0e-4
    -t ${TEMP}/${MODULE_NAME}Test2.mha
    -b MIDAS{${MODULE_NAME}Test2.mha.md5}
    -n 3 )
set_property( TEST ${MODULE_NAME}-Test2-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2 )

# Test3
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3
  COMMAND ${PROJ_EXE}
    --interpolator NearestNeighbor
    --origin 180,180
    --spacing 3,3
    MIDAS{ES0015_Large.mha.md5}
    ${TEMP}/${MODULE_NAME}Test3.mha )

# Test3-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test3-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}Test3.mha
    -b MIDAS{${MODULE_NAME}Test3.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test3-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test3 )

# Test4
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4
  COMMAND ${PROJ_EXE}
    --matchImage MIDAS{ES0015_1.mha.md5}
    MIDAS{ES0015_Large.mha.md5}
    ${TEMP}/${MODULE_NAME}Test4.mha )

# Test4-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test4-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}Test4.mha
    -b MIDAS{${MODULE_NAME}Test4.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test4-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test4 )

# Test5
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test5
  COMMAND ${PROJ_EXE}
    --interpolator Sinc
    --makeHighResIso
    MIDAS{CroppedWholeLungCTScan.mhd.md5}
    ${TEMP}/${MODULE_NAME}Test5.mha
    MIDAS_FETCH_ONLY{CroppedWholeLungCTScan.raw.md5} )

# Test5-Compare
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test5-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}Test5.mha
    -b MIDAS{${MODULE_NAME}Test5.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test5-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test5 )
