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

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( PROJ_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:${MODULE_NAME}> )

# Test1 - VesselMask
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-VesselMask
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test1.mha}
    ${TEMP}/${MODULE_NAME}-Test1.mha
    --gap 0.5
    --notVesselWidth 2 )

# Test1 - Compare - VesselMask
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test1-VesselMask-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test1.mha
    -b DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test1_VesselMask.mha}
    -i 0.01 )
set_tests_properties( ${MODULE_NAME}-Test1-VesselMask-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Test1-VesselMask )

# Test2 - NotVesselMask
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-NotVesselMask
  COMMAND ${PROJ_EXE}
    DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test2.mha}
    ${TEMP}/${MODULE_NAME}-Test2.mha
    --gap 0.5
    --notVesselWidth 2
    --notVesselMask ${TEMP}/${MODULE_NAME}-Test2_NotVesselMask.mha )

# Test2 - Compare - NotVesselMask
ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-NotVesselMask-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test2.mha
    -b DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test2_NotVesselMask.mha}
    -i 0.01 )
set_tests_properties( ${MODULE_NAME}-Test2-NotVesselMask-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Test2-NotVesselMask )

ExternalData_Add_Test( TubeTKData
  NAME ${MODULE_NAME}-Test2-NotVesselMask1-Compare
  COMMAND ${TubeTK_CompareImages_EXE}
    -t ${TEMP}/${MODULE_NAME}-Test2_NotVesselMask.mha
    -b DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test2_NotVesselMask1.mha}
    -i 0.01 )
set_tests_properties( ${MODULE_NAME}-Test2-NotVesselMask1-Compare
  PROPERTIES DEPENDS ${MODULE_NAME}-Test2-NotVesselMask )
