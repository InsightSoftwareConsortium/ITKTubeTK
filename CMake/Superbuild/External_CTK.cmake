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

set( proj CTK )

# Sanity checks.
if( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )
  message( FATAL_ERROR "${proj}_DIR variable is defined but corresponds to a nonexistent directory" )
endif( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )

# Set dependency list
set( ${proj}_DEPENDENCIES "ITK" )
if( TubeTK_USE_VTK )
  set( ${proj}_DEPENDENCIES ${${proj}_DEPENDENCIES} "VTK" )
endif( TubeTK_USE_VTK )

# Include dependent projects, if any.
ExternalProject_Include_Dependencies( ${proj}
  PROJECT_VAR proj
  DEPENDS_VAR ${proj}_DEPENDENCIES
  USE_SYSTEM_VAR USE_SYSTEM_${proj}
  SUPERBUILD_VAR TubeTK_USE_SUPERBUILD )

if( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )
  set( ${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj} )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )

  set( ${proj}_use_git_protocol ON )
  if( GIT_PROTOCOL_HTTP )
    set( ${proj}_use_git_protocol OFF )
  endif( GIT_PROTOCOL_HTTP )

  set( TubeTK_CTK_VTK_ARGS )
  if( TubeTK_USE_VTK )
    list( APPEND TubeTK_CTK_VTK_ARGS
      -DVTK_DIR:PATH=${VTK_DIR}
      -DCTK_USE_SYSTEM_VTK:BOOL=OFF
      -DCTK_LIB_Visualization/VTK/Widgets:BOOL=ON
      -DCTK_LIB_Visualization/VTK/Widgets_USE_TRANSFER_FUNCTION_CHARTS:BOOL=ON )
  endif( TubeTK_USE_VTK )

  if( CMAKE_VERSION VERSION_LESS 3.4 )
    set( CMAKE_CXX_COMPILER_LAUNCHER_FLAG )
    set( CMAKE_C_COMPILER_LAUNCHER_FLAG )
  else()
    set( CMAKE_CXX_COMPILER_LAUNCHER_FLAG
      -DCMAKE_CXX_COMPILER_LAUNCHER:FILEPATH=${CMAKE_CXX_COMPILER_LAUNCHER} )
    set( CMAKE_C_COMPILER_LAUNCHER_FLAG
      -DCMAKE_C_COMPILER_LAUNCHER:FILEPATH=${CMAKE_C_COMPILER_LAUNCHER} )
  endif()

  ExternalProject_Add( ${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY ${${proj}_URL}
    GIT_TAG ${${proj}_HASH_OR_TAG}
    DOWNLOAD_DIR ${${proj}_SOURCE_DIR}
    SOURCE_DIR ${${proj}_SOURCE_DIR}
    BINARY_DIR ${${proj}_DIR}
    INSTALL_DIR ${${proj}_DIR}
    LOG_DOWNLOAD 1
    LOG_UPDATE 0
    LOG_CONFIGURE 0
    LOG_BUILD 0
    LOG_TEST 0
    LOG_INSTALL 0
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      ${CMAKE_CXX_COMPILER_LAUNCHER_FLAG}
      ${CMAKE_C_COMPILER_LAUNCHER_FLAG}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_TESTING:BOOL=OFF
      -DCTK_LIB_ImageProcessing/ITK/Core:BOOL=ON
      -DCTK_LIB_PluginFramework:BOOL=OFF
      -DCTK_LIB_Widgets:BOOL=ON
      -DCTK_PLUGIN_org.commontk.eventbus:BOOL=OFF
      -DCTK_USE_GIT_PROTOCOL:BOOL=${${proj}_use_git_protocol}
      -DITK_DIR:PATH=${ITK_DIR}
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      ${TubeTK_CTK_VTK_ARGS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES} )

else( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )
  if( ${USE_SYSTEM_${proj}} )
    find_package( ${proj} REQUIRED )
  endif( ${USE_SYSTEM_${proj}} )

  ExternalProject_Add_Empty( ${proj} DEPENDS ${${proj}_DEPENDENCIES} )
endif( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )

list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS -D${proj}_DIR:PATH=${${proj}_DIR} )
