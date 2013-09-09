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

# Make sure this file is included only once.
get_filename_component( CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE}
  NAME_WE )
if( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED )
  return()
endif( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED )
set( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1 )

set( proj KWStyle )

# Sanity checks.
if( KWSTYLE_EXECUTABLE AND NOT EXISTS ${KWSTYLE_EXECUTABLE} )
  message( FATAL_ERROR "The KWSTYLE_EXECUTABLE variable is defined, but corresponds to a nonexistent file" )
endif( KWSTYLE_EXECUTABLE AND NOT EXISTS ${KWSTYLE_EXECUTABLE} )

set( ${proj}_DEPENDENCIES "" )

# Include dependent projects, if any.
TubeTKMacroCheckExternalProjectDependency( ${proj} )

if( NOT KWSTYLE_EXECUTABLE AND NOT ${USE_SYSTEM_KWSTYLE} )
  set( ${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj} )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )
  set( KWSTYLE_EXECUTABLE ${${proj}_DIR}/bin/KWStyle )

  ExternalProject_Add( ${proj}
    GIT_REPOSITORY ${${proj}_URL}
    GIT_TAG ${${proj}_URL_HASH}
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
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES} )

else( NOT KWSTYLE_EXECUTABLE AND NOT ${USE_SYSTEM_KWSTYLE} )
  if( ${USE_SYSTEM_KWSTYLE} )
    find_package( ${proj} REQUIRED )
  endif( ${USE_SYSTEM_KWSTYLE} )

  TubeTKMacroEmptyExternalProject( ${proj} "${${proj}_DEPENDENCIES}" )
endif( NOT KWSTYLE_EXECUTABLE AND NOT ${USE_SYSTEM_KWSTYLE} )

list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS -DKWSTYLE_EXECUTABLE:FILEPATH=${KWSTYLE_EXECUTABLE} )
