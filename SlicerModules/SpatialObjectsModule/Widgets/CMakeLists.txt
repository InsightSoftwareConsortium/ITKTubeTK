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
project(qSlicer${MODULE_NAME}ModuleWidgets)

set(KIT ${PROJECT_NAME})

set(${KIT}_EXPORT_DIRECTIVE "Q_SLICER_MODULE_${MODULE_NAME_UPPER}_WIDGETS_EXPORT")

set(${KIT}_INCLUDE_DIRECTORIES
  ${CMAKE_BUILD_DIR}
  ${vtkSlicerSpatialObjectsModuleMRML_SOURCE_DIR}
  ${vtkSlicerSpatialObjectsModuleMRML_BINARY_DIR}
  )

set(${KIT}_SRCS
  qMRMLSpatialObjectsTreeView.h
  qMRMLSpatialObjectsTreeView.cxx
  qMRMLSceneSpatialObjectsModel.h
  qMRMLSceneSpatialObjectsModel.cxx
  qSlicerSpatialObjectsModuleWidget.cxx
  qSlicerSpatialObjectsModuleWidget.h
  qSlicerSpatialObjectsBasicWidget.cxx
  qSlicerSpatialObjectsBasicWidget.h
  qSlicerSpatialObjectsGlyphWidget.cxx
  qSlicerSpatialObjectsGlyphWidget.h
  qSlicerSpatialObjectsWidget.cxx
  qSlicerSpatialObjectsWidget.h
  )

set(${KIT}_MOC_SRCS
  qSlicerSpatialObjectsModuleWidget.h
  qSlicerSpatialObjectsBasicWidget.h
  qSlicerSpatialObjectsGlyphWidget.h
  qSlicerSpatialObjectsWidget.h
  qMRMLSpatialObjectsTreeView.h
  qMRMLSceneSpatialObjectsModel.h
  )

set(${KIT}_UI_SRCS
  Resources/UI/qSlicerSpatialObjectsModule.ui
  Resources/UI/qSlicerSpatialObjectsBasicWidget.ui
  Resources/UI/qSlicerSpatialObjectsGlyphWidget.ui
  Resources/UI/qSlicerSpatialObjectsWidget.ui
  )

set(${KIT}_RESOURCES
  )

set(${KIT}_TARGET_LIBRARIES
  vtkSlicerSpatialObjectsModuleLogic
  ${MRML_LIBRARIES}
  vtkSlicerAnnotationsModuleMRML
  )

#-----------------------------------------------------------------------------
SlicerMacroBuildModuleWidgets(
  NAME ${KIT}
  EXPORT_DIRECTIVE ${${KIT}_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${${KIT}_INCLUDE_DIRECTORIES}
  SRCS ${${KIT}_SRCS}
  MOC_SRCS ${${KIT}_MOC_SRCS}
  UI_SRCS ${${KIT}_UI_SRCS}
  TARGET_LIBRARIES ${${KIT}_TARGET_LIBRARIES}
  RESOURCES ${${KIT}_RESOURCES}
  WRAP_PYTHONQT
  )

# Plugin
#add_subdirectory(Plugins)
