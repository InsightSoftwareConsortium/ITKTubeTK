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

project( ImageEditor )

set( MODULE_NAME ${PROJECT_NAME} )

find_package( Qt4 COMPONENTS QtCore QtGui QtNetwork QtOpenGL QtWebkit REQUIRED )
include( ${QT_USE_FILE} )

find_package( ITK REQUIRED )
if( TubeTK_BUILD_WITHIN_SLICER )
  set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
endif( TubeTK_BUILD_WITHIN_SLICER )
include( ${ITK_USE_FILE} )

find_package(OpenGL REQUIRED )

find_package( ImageViewer REQUIRED )
include( ${ImageViewer_USE_FILE} )

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

set( ${MODULE_NAME}_SRCS
  QtImageEditor.cxx
  QtOverlayControlsWidget.cxx
)

set( ${MODULE_NAME}_MOC_SRCS
  QtImageEditor.h
  QtOverlayControlsWidget.h
)

set( ${MODULE_NAME}_GUI_SRCS
  Resources/UI/QtOverlayControlsWidget.ui
  )

# \todo: expose resource file into ImageViewer package
set( ${MODULE_NAME}_RESOURCES
#  Resources/qtImageEditorResources.qrc
  ${ImageViewer_SOURCE_DIR}/QtImageViewer/Resources/qtImageViewerResources.qrc
  )

QT4_WRAP_UI( ${MODULE_NAME}_UI_CXX ${${MODULE_NAME}_GUI_SRCS} )
QT4_WRAP_CPP( ${MODULE_NAME}_SRCS ${${MODULE_NAME}_MOC_SRCS} )
QT4_ADD_RESOURCES( ${MODULE_NAME}_QRC_SRCS ${${MODULE_NAME}_RESOURCES} )

#set( ${MODULE_NAME}_H_FILES
#  itkGaussianDerivativeImageSource.h
#  )

#set( ${MODULE_NAME}_HXX_FILES
#  itkGaussianDerivativeImageSource.hxx
#  )

add_custom_target( ${MODULE_NAME}Include SOURCES
  ${${MODULE_NAME}_H_FILES}
  ${${MODULE_NAME}_HXX_FILES} )

set( ${MODULE_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${ImageViewer_INCLUDE_DIRS}
  ${QT_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
)

set(BUILD_EXEC_ONLY EXECUTABLE_ONLY)
if(BUILD_SHARED_LIBS)
  set(BUILD_EXEC_ONLY)
endif(BUILD_SHARED_LIBS)

SEMMacroBuildCLI(
  ${BUILD_EXEC_ONLY}
  NAME ${MODULE_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  ADDITIONAL_SRCS
    ${${MODULE_NAME}_SRCS}
    ${${MODULE_NAME}_GUI_SRCS}
    ${${MODULE_NAME}_QRC_SRCS}
    ${${MODULE_NAME}_UI_CXX}
  TARGET_LIBRARIES
    ${ImageViewer_LIBRARIES}
    ${ITK_LIBRARIES}
    ${QT_LIBRARIES}
    ${OPENGL_LIBRARIES}
    TubeTKFiltering
  INCLUDE_DIRECTORIES
    ${${MODULE_NAME}_INCLUDE_DIRS}
  INSTALL_RUNTIME_DESTINATION bin )

# Configure ImageEditorConfigure.h.
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/ImageEditorConfigure.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/ImageEditorConfigure.h @ONLY )
