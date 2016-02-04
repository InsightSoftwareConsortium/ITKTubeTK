##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0(the "License");
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

project( ConvertToMetaImage )

cmake_minimum_required( VERSION 2.4 )
if( COMMAND CMAKE_POLICY )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND CMAKE_POLICY )

include_directories( ${TubeTK_SOURCE_DIR}/Base/Common )

if( TubeTK_USE_QT )

  # Find QT
  find_package( Qt4 REQUIRED )
  include( ${QT_USE_FILE} )

  include_directories(
    ${QT_INCLUDE_DIR}
    ${QT_QTGUI_INCLUDE_DIR}
    ${QT_QTCORE_INCLUDE_DIR}
    ${CMAKE_CURRENT_BINARY_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR} )

  set( MetaImageImporter_SRCS
    main.cxx
    MetaImageImporterUI.h
    MetaImageImporterUI.cxx )

  set( MetaImageImporter_MOC_SRCS
    MetaImageImporterUI.h )
  set( MetaImageImporter_UI_SRCS MetaImageImporterUI.ui )

  QT4_WRAP_CPP( MetaImageImporter_SRCS ${MetaImageImporter_MOC_SRCS} )
  QT4_WRAP_UI( MetaImageImporter_UI_CXX ${MetaImageImporter_UI_SRCS} )

  option( HIDE_CONSOLE "Hide the console" OFF )
  if( HIDE_CONSOLE )
    set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows" )
    link_libraries( ${QT_QTMAIN_LIBRARY} )
  endif( HIDE_CONSOLE )

  add_executable( ConvertToMetaImage
    ${MetaImageImporter_SRCS} ${MetaImageImporter_UI_CXX} )

  link_libraries(
    ${QT_QTCORE_LIBRARY}
    ${QT_QTGUI_LIBRARY} )

  target_link_libraries( ConvertToMetaImage
    ${QT_QTCORE_LIBRARY}
    ${QT_QTGUI_LIBRARY}
    ${QT_QTDESIGNER_LIBRARY} )

  # Installation file
  install( TARGETS
    ConvertToMetaImage
    DESTINATION bin
    COMPONENT RUNTIME )

endif( TubeTK_USE_QT )

add_executable( ConvertToMetaImageCLI ConvertToMetaImageCLI.cxx )
target_link_libraries(  ConvertToMetaImageCLI TubeTKCommon )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
