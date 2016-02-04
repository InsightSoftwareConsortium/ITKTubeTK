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

set( MODULE_NAME Tortuosity )
set( MODULE_TITLE "Tortuosity" )

string( TOUPPER ${MODULE_NAME} MODULE_NAME_UPPER )

#find_package( TubeTK REQUIRED )
#include( ${TubeTK_USE_FILE} )

add_subdirectory( Logic )
add_subdirectory( Widgets )

set( MODULE_EXPORT_DIRECTIVE "Q_SLICER_QTMODULES_${MODULE_NAME_UPPER}_EXPORT" )

# Additional includes - Current_{source,binary} and Slicer_{Libs,Base} already included
set( MODULE_INCLUDE_DIRECTORIES
  ${CMAKE_CURRENT_SOURCE_DIR}/Logic
  ${CMAKE_CURRENT_BINARY_DIR}/Logic
  ${CMAKE_CURRENT_SOURCE_DIR}/Widgets
  ${CMAKE_CURRENT_BINARY_DIR}/Widgets
  )

set( MODULE_SRCS
  qSlicer${MODULE_NAME}Module.cxx
  qSlicer${MODULE_NAME}Module.h )

set( MODULE_MOC_SRCS
  qSlicer${MODULE_NAME}Module.h )

set( MODULE_UI_SRCS )

set( MODULE_TARGET_LIBRARIES
  vtkSlicer${MODULE_NAME}ModuleLogic
  qSlicer${MODULE_NAME}ModuleWidgets
  ${MRML_LIBRARIES}
  TubeTKCommon
  TubeTKFiltering
  TubeTKNumerics )

set( MODULE_RESOURCES
  Resources/qSlicer${MODULE_NAME}Module.qrc )

slicerMacroBuildQtModule(
  NAME ${MODULE_NAME}
  TITLE ${MODULE_TITLE}
  EXPORT_DIRECTIVE ${MODULE_EXPORT_DIRECTIVE}
  INCLUDE_DIRECTORIES ${MODULE_INCLUDE_DIRECTORIES}
  SRCS ${MODULE_SRCS}
  MOC_SRCS ${MODULE_MOC_SRCS}
  UI_SRCS ${MODULE_UI_SRCS}
  TARGET_LIBRARIES ${MODULE_TARGET_LIBRARIES}
  RESOURCES ${MODULE_RESOURCES} )

# Testing is delegated to Applications/App/Testing/Python
#if( BUILD_TESTING )
#  add_subdirectory( Testing )
#endif( BUILD_TESTING )
