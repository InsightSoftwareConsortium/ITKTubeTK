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

project( vtkSlicer${MODULE_NAME}ModuleLogic )

set( KIT ${PROJECT_NAME} )

find_package( TubeTK REQUIRED )
include( ${TubeTK_USE_FILE} )

set( ${KIT}_EXPORT_DIRECTIVE "VTK_SLICER_${MODULE_NAME_UPPER}_MODULE_LOGIC_EXPORT" )

set( ${KIT}_INCLUDE_DIRECTORIES
  ${MRML_LIBRARIES}
  ${vtkSlicer${MODULE_NAME}ModuleMRML_SOURCE_DIR}
  ${vtkSlicer${MODULE_NAME}ModuleMRML_BINARY_DIR}
  ${vtkSlicerSpatialObjectsModuleMRML_SOURCE_DIR}
  ${vtkSlicerSpatialObjectsModuleMRML_BINARY_DIR}
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  )

set( ${KIT}_SRCS
  vtkSlicerTortuosityLogic.cxx
  vtkSlicerTortuosityLogic.h
  )

set( ${KIT}_TARGET_LIBRARIES
  ${ITK_LIBRARIES}
  ${MRML_LIBRARIES}
  vtkSlicerSpatialObjectsModuleMRML
  )

SlicerMacroBuildModuleLogic(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES} )
