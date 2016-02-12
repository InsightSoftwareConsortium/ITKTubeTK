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

Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1
  COMMAND ${PROJ_EXE}
    --startIndex 0
    --endIndex 2
    --incrementIndex 2
    MIDAS{reexported_tracking_data_f_Kitware_v2.txt.md5}
    ${TEMP}/${MODULE_NAME}TestOutput.mha
    MIDAS_FETCH_ONLY{ultrasound_0002392.ppm.md5}
    MIDAS_FETCH_ONLY{ultrasound_0002393.ppm.md5}
    MIDAS_FETCH_ONLY{ultrasound_0002394.ppm.md5} )

Midas3FunctionAddTest( NAME ${MODULE_NAME}-Test1-Compare
  COMMAND ${IMAGECOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}TestOutput.mha
    -b MIDAS{${MODULE_NAME}TestOutput.mha.md5} )
set_property( TEST ${MODULE_NAME}-Test1-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )
