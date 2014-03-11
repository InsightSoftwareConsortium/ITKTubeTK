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

cmake_minimum_required( VERSION 2.8.7 )

project( ImageEditor )

set( MODULE_NAME ${PROJECT_NAME} )

find_package( Qt4 COMPONENTS QtCore QtGui QtNetwork QtOpenGL QtWebkit REQUIRED )
include( ${QT_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

find_package(OpenGL REQUIRED )

find_package( ImageViewer REQUIRED )
include( ${ImageViewer_USE_FILE} )

set(${MODULE_NAME}_INCLUDE_DIRS
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${ImageViewer_INCLUDE_DIRS}
  ${QT_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
)
include_directories(
  ${${MODULE_NAME}_INCLUDE_DIRS}
)

set(${MODULE_NAME}_MOC_SRCS
  QtControlView.h
)

set( ${MODULE_NAME}_SRCS
#  ImageEditor.cxx
  QtControlView.cxx
)

QT4_WRAP_CPP( ${MODULE_NAME}_SRCS ${${MODULE_NAME}_MOC_SRCS} )

add_library( ${MODULE_NAME}Lib SHARED ${${MODULE_NAME}_SRCS} )
target_link_libraries( ${MODULE_NAME}Lib
  ${ImageViewer_LIBRARIES}
  ${ITK_LIBRARIES}
  ${QT_LIBRARIES}
  ${OPENGL_LIBRARIES}
)

add_executable( ${MODULE_NAME} ImageEditor.cxx )
target_link_libraries( ${MODULE_NAME}
  ${MODULE_NAME}Lib
)
