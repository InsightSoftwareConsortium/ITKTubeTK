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
  ${SlicerExecutionModel_CLI_RUNTIME_OUTPUT_DIRECTORY}/${MODULE_NAME} )

set( IMAGECOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/ImageCompareCommand )

set( TEXTCOMPARE_EXE
  ${TubeTK_BINARY_DIR}/bin/TextCompareCommand )

#######
# Test1
#######
Midas3FunctionAddTest( ${MODULE_NAME}-Test1
  ${PROJ_EXE}
    --tubeId 255
    --backgroundId 0
    --voidId 10
    --tubeScales 0.3,1.5
    --saveDiscriminantInfo ${TEMP}/${MODULE_NAME}Test1.mrs
    --labelmap MIDAS{im0001.vk.crop.mha.md5}
    MIDAS{im0001.crop.contrast.mha.md5}
    ${TEMP}/${MODULE_NAME}Test1_Output.mha )

#######
# Test1-Compare MRS
#######
Midas3FunctionAddTest( ${MODULE_NAME}-Test1-MRS-Compare
  ${TEXTCOMPARE_EXE}
    -t ${TEMP}/${MODULE_NAME}Test1.mrs
    -b MIDAS{${MODULE_NAME}Test1.mrs.md5} )
set_property(TEST ${MODULE_NAME}-Test1-MRS-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )

#######
# Test1-Compare Output
#######
Midas3FunctionAddTest( ${MODULE_NAME}-Test1-Output-Compare
  ${IMAGECOMPARE_EXE}
    -i 0.01
    -t ${TEMP}/${MODULE_NAME}Test1_Output.mha
    -b MIDAS{${MODULE_NAME}Test1_Output.mha.md5} )
set_property(TEST ${MODULE_NAME}-Test1-Output-Compare
  APPEND PROPERTY DEPENDS ${MODULE_NAME}-Test1 )
