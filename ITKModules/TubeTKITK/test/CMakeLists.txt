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

set( CompareImages_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )

######################  PythonComputeTrainingMask  ##########################
if( ITK_WRAP_PYTHON )
  # Test1 - PythonComputeTrainingMask
  set( MODULE_NAME PythonComputeTrainingMask )
  Midas3FunctionAddTestWithEnv( NAME TubeTKITK.Python.ComputeTrainingMask-Test1
    COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
        MIDAS{ComputeTrainingMask-Test1.mha.md5}
        ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
        ${TEMP}/${MODULE_NAME}_NotVesselMaskTest1.mha
        0.5
        2
    ENVIRONMENT ITK_BUILD_DIR=${ITK_DIR} TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR} BUILD_TYPE=$<$<BOOL:${WIN32}>:$<CONFIG>>
    )
  # Test1 - Compare - PythonComputeTrainingMask
  Midas3FunctionAddTest( NAME
    TubeTKITK.Python.ComputeTrainingMask-Test1-VesselMask-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
      -b MIDAS{ComputeTrainingMask-Test1_VesselMask.mha.md5}
      -i 0.01
    )
  set_property( TEST
    TubeTKITK.Python.ComputeTrainingMask-Test1-VesselMask-Compare
    APPEND PROPERTY DEPENDS TubeTKITK.Python.ComputeTrainingMask-Test1 )

endif()
