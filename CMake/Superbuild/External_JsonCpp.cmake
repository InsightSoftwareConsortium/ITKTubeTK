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

set( proj JsonCpp )

# Sanity checks.
if( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )
  message( FATAL_ERROR
    "${proj}_DIR variable is defined as ${${proj}_DIR} which corresponds to a nonexistent directory" )
endif( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )

set( ${proj}_DEPENDENCIES "" )

# Include dependent projects, if any.
ExternalProject_Include_Dependencies( ${proj}
  PROJECT_VAR proj
  DEPENDS_VAR ${proj}_DEPENDENCIES
  USE_SYSTEM_VAR USE_SYSTEM_${proj}
  SUPERBUILD_VAR TubeTK_USE_SUPERBUILD )

if( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_JSONCPP} )
  set( ${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj} )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )

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
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DJSONCPP_WITH_TESTS:BOOL=OFF
      -DCMAKE_INSTALL_PREFIX:PATH=${${proj}_DIR}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS} )

else( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_JSONCPP} )
  if( ${USE_SYSTEM_JSONCPP} )
    find_package( ${proj} REQUIRED )
  endif( ${USE_SYSTEM_JSONCPP} )

  ExternalProject_Add_Empty( ${proj} DEPENDS ${${proj}_DEPENDENCIES} )
endif( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_JSONCPP} )

list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS -D${proj}_DIR:PATH=${${proj}_DIR} )
