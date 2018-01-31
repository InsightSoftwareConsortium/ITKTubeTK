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

set( TEMP ${TubeTK_BINARY_DIR}/Temporary )

set( CompareImages_EXE ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareImages> )
set( CompareTextFiles_EXE
  ${TubeTK_LAUNCHER} $<TARGET_FILE:CompareTextFiles> )

######################  PythonComputeTrainingMask  ##########################
if( ITK_WRAP_PYTHON )
  # Test1 - PythonComputeTrainingMask
  set( MODULE_NAME PythonComputeTrainingMask )
  ExternalData_Add_Test_With_Env(
    NAME
      TubeTKITK.Python.ComputeTrainingMask-Test1
    COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
        DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test1.mha}
        ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
        ${TEMP}/${MODULE_NAME}_NotVesselMaskTest1.mha
        0.5
        2
    ENVIRONMENT
      ITK_BUILD_DIR=${ITK_DIR}
      TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR}
      BUILD_TYPE=$<$<BOOL:${WIN32}>:$<CONFIG>>
    )
  # Test1 - Compare - PythonComputeTrainingMask
  ExternalData_Add_Test( TubeTKData
    NAME
      TubeTKITK.Python.ComputeTrainingMask-Test1-VesselMask-Compare
    COMMAND ${CompareImages_EXE}
      -t ${TEMP}/${MODULE_NAME}_VesselMaskTest1.mha
      -b DATA{${TubeTK_DATA_ROOT}/ComputeTrainingMask-Test1_VesselMask.mha}
      -i 0.01
    )
  set_tests_properties(
    TubeTKITK.Python.ComputeTrainingMask-Test1-VesselMask-Compare
    PROPERTIES DEPENDS TubeTKITK.Python.ComputeTrainingMask-Test1 )

######################  PythonConvertImagesToCSV  ########################
  # Test1 - PythonConvertImagesToCSV
  set( MODULE_NAME PythonConvertImagesToCSV )
  ExternalData_Add_Test_With_Env(
    NAME
      TubeTKITK.Python.ConvertImagesToCSV-Test1
    COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large-TrainingMask.mha}
        DATA{${TubeTK_DATA_ROOT}/GDS0015_Large.mha},DATA{${TubeTK_DATA_ROOT}/ES0015_Large.mha}
        ${TEMP}/${MODULE_NAME}Test1.csv
        3
    ENVIRONMENT
      ITK_BUILD_DIR=${ITK_DIR}
      TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR}
      BUILD_TYPE=$<$<BOOL:${WIN32}>:$<CONFIG>>
    )

  # Test1 - Compare - PythonConvertImagesToCSV
  ExternalData_Add_Test( TubeTKData
    NAME
      TubeTKITK.Python.ConvertImagesToCSV-Test1-Compare
    COMMAND ${CompareTextFiles_EXE}
      -t ${TEMP}/${MODULE_NAME}Test1.csv
      -b DATA{${TubeTK_DATA_ROOT}/ConvertImagesToCSVTest1.csv}
    )
  set_tests_properties( TubeTKITK.Python.ConvertImagesToCSV-Test1-Compare
    PROPERTIES DEPENDS TubeTKITK.Python.ConvertImagesToCSV-Test1 )

##############  PythonConvertShrunkenSeedImageToList  ###################
  # Test1 - PythonConvertShrunkenSeedImageToList
  set( MODULE_NAME PythonConvertShrunkenSeedImageToList )
  ExternalData_Add_Test_With_Env(
    NAME
      TubeTKITK.Python.ConvertShrunkenSeedImageToList-Test1
    COMMAND ${PYTHON_EXECUTABLE}
      ${CMAKE_CURRENT_SOURCE_DIR}/${MODULE_NAME}Test.py
        DATA{${TubeTK_DATA_ROOT}/LiverVess_Shrink_Seeds.mha}
        DATA{${TubeTK_DATA_ROOT}/LiverVess_Shrink_SeedScales.mha}
        DATA{${TubeTK_DATA_ROOT}/LiverVess_Shrink_Points.mha}
        ${TEMP}/${MODULE_NAME}Test1.sdlst
        0
    ENVIRONMENT
      ITK_BUILD_DIR=${ITK_DIR} 
      TubeTK_BUILD_DIR=${TubeTK_BINARY_DIR} 
      BUILD_TYPE=$<$<BOOL:${WIN32}>:$<CONFIG>>
    )
  # Test1 - Compare - PythonConvertShrunkenSeedImageToList
  ExternalData_Add_Test( TubeTKData
    NAME
      TubeTKITK.Python.ConvertShrunkenSeedImageToList-Test1-Compare
    COMMAND ${CompareTextFiles_EXE}
      -t ${TEMP}/${MODULE_NAME}Test1.sdlst
      -b DATA{${TubeTK_DATA_ROOT}/ConvertShrunkenSeedImageToListTest1.sdlst}
    )
  set_tests_properties(
    TubeTKITK.Python.ConvertShrunkenSeedImageToList-Test1-Compare
    PROPERTIES DEPENDS TubeTKITK.Python.ConvertShrunkenSeedImageToList-Test1 )
endif()
