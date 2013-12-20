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

set( proj ITK )

# Sanity checks.
if( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )
  message( FATAL_ERROR "${proj}_DIR variable is defined but corresponds to a nonexistent directory" )
endif( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )

set( ${proj}_DEPENDENCIES "" )

# Include dependent projects, if any.
TubeTKMacroCheckExternalProjectDependency( ${proj} )

if( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )
  set( ${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj} )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )

  set( TubeTK_ITKHDF5_VALGRIND_ARGS )
  if( TubeTK_USE_VALGRIND )
    list( APPEND TubeTK_ITKHDF5_VALGRIND_ARGS
      -DHDF5_ENABLE_USING_MEMCHECKER:BOOL=ON )
  endif( TubeTK_USE_VALGRIND )

  ExternalProject_Add( ${proj}
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
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DITKV3_COMPATIBILITY:BOOL=ON
      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
      -DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
      -DITK_LEGACY_REMOVE:BOOL=OFF
      -DKWSYS_USE_MD5:BOOL=ON
      -DModule_ITKReview:BOOL=ON
      ${TubeTK_ITKHDF5_VALGRIND_ARGS}
    INSTALL_COMMAND "" )

else( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )
  if( ${USE_SYSTEM_${proj}} )
    find_package( ${proj} REQUIRED )
  endif( ${USE_SYSTEM_${proj}} )

  TubeTKMacroEmptyExternalProject( ${proj} "${${proj}_DEPENDENCIES}" )
endif( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )

list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS -D${proj}_DIR:PATH=${${proj}_DIR} )
