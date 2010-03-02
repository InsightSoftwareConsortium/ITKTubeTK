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

set( APPLICATION_PATH ${TubeTK_BINARY_DIR}/Utilities )

set( INPUT_DATA ${TubeTK_SOURCE_DIR}/Data )
set( ITK_INPUT_DATA ${ITK_DATA_ROOT}/Input )

set( BASELINE ${APPLICATION_PATH}/Testing/Baselines )

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

file( COPY ${INPUT_DATA}/slice256.png DESTINATION ${TEMP} )
file( COPY ${INPUT_DATA}/im0001.crop.mha DESTINATION ${TEMP} )

add_test( ImageCompareCommandTest1
            ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand
               -t ${TEMP}/slice256.png
               -b ${INPUT_DATA}/slice256.png )

add_test( ImageCompareCommandTest2
            ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand
               -n 1000
               -r 1
               -i 50
               -t ${TEMP}/im0001.crop.mha
               -b ${INPUT_DATA}/im0001_n10.crop.mha )

add_test( ImageCompareCommandTest3
            ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand
               -t ${TEMP}/im0001.mha
               -b ${INPUT_DATA}/im0001_n10.crop.mha )
set_tests_properties( ImageCompareCommandTest3
                      PROPERTIES WILL_FAIL true )

add_test( ImageCompareCommandTest4
            ${TubeTK_BINARY_DIR}/${Slicer3_INSTALL_PLUGINS_BIN_DIR}/ImageCompareCommand
               -t ${TEMP}/im0001.crop.mha
               -B 2 ${INPUT_DATA}/im0001_n10.crop.mha ${INPUT_DATA}/im0001.crop.mha )

