##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 ( the "License" );
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
PROJECT( tubetkImageViewer )

INCLUDE ( ${CMAKE_ROOT}/Modules/FindOpenGL.cmake )

SET( QT_USE_QTOPENGL TRUE )
INCLUDE( ${QT_USE_FILE} )

MESSAGE( ${QT_INCLUDE_DIR} )

INCLUDE_DIRECTORIES(
 ${QT_INCLUDE_DIR}
 ${QT_INCLUDE_DIR}/Qt
 ${tubetkImageViewer_BINARY_DIR}
 ${tubetkImageViewer_SOURCE_DIR}
 )

LINK_LIBRARIES(
 ${ITK_LIBRARIES}
 ${QT_QT_LIBRARY}
 ${QT_GL_LIBRARY}
 ${OPENGL_LIBRARY}
 )

SET( tubetkImageViewer_SRCS
  QtGlSliceView.cxx
  tubetkImageViewer.cxx
  tubetkImageViewerTest.cxx
 )

QT4_AUTOMOC( ${tubetkImageViewer_SRCS} )

QT4_WRAP_UI( tubetkImageViewer_UIS_H
  tubetkImageViewerGUI.ui
  tubetkImageViewerHelpGUI.ui
  )

QT4_WRAP_CPP( tubetkImageViewer_MOC_SRCS
  QtGlSliceView.h
  )

source_group( "Generated" FILES
  ${tubetkImageViewer_UIS_H}
  ${tubetkImageViewer_MOC_SRCS}
  )

ADD_DEFINITIONS( -DQT_DLL )

ADD_DEFINITIONS( -DQT_THREAD_SUPPORT )

ADD_EXECUTABLE( ${PROJECT_NAME}
  ${tubetkImageViewer_SRCS}
  ${tubetkImageViewer_UIS_H}
  ${tubetkImageViewer_MOC_SRCS}
  )

TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
  ${QT_LIBRARIES}
  )

slicer3_set_plugins_output_path( ${PROJECT_NAME} )
slicer3_install_plugins( ${PROJECT_NAME} )
