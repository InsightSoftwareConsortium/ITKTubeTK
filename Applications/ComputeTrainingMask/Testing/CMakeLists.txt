############################################################################
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
############################################################################
include_regular_expression( "^.*$" )

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

set( CompareImages_EXE
 ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

# Test1 - VesselMask
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-VesselMask
  COMMAND ${PROJ_EXE}
    MIDAS{ComputeTrainingMask-Test1.mha.md5}
    ${TEMP}/${MODULE_NAME}Test1.mha
    --gap 0.5
    --notVesselWidth 2)

# Test1 - Compare - VesselMask
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-VesselMask-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test1.mha
    -b MIDAS{ComputeTrainingMask-Test1_VesselMask.mha.md5}
    -i 0.01 )

set_property( TEST ${MODULE_NAME}-Test1-VesselMask-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1-VesselMask )

# Test2 - NotVesselMask
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-NotVesselMask
  COMMAND ${PROJ_EXE}
    MIDAS{ComputeTrainingMask-Test2.mha.md5}
    ${TEMP}/${MODULE_NAME}Test2.mha
    --gap 0.5
    --notVesselWidth 2
    --notVesselMask ${TEMP}/${MODULE_NAME}Test2_NotVesselMask.mha)

# Test2 - Compare - NotVesselMask
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-NotVesselMask-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test2.mha
    -b MIDAS{ComputeTrainingMask-Test2_NotVesselMask.mha.md5}
    -i 0.01 )
Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test2-NotVesselMask1-Compare
  COMMAND ${CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}Test2_NotVesselMask.mha
    -b MIDAS{ComputeTrainingMask-Test2_NotVesselMask1.mha.md5}
    -i 0.01 )

set_property( TEST ${MODULE_NAME}-Test2-NotVesselMask-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test2-NotVesselMask )
