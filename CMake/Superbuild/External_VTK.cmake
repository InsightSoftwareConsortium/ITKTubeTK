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

set( proj VTK )

# Sanity checks.
if( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )
  message( FATAL_ERROR "${proj}_DIR variable is defined but corresponds to a nonexistent directory" )
endif( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )

if( TubeTK_USE_JsonCpp )
  set( ${proj}_DEPENDENCIES "JsonCpp" )
  set( ${proj}_DEPENDENCIES_ARGS "-DJsonCpp_DIR:PATH=${JsonCpp_DIR}"
    "-DJsonCpp_LIBRARIES:PATH=${JsonCpp_LIBRARIES}"
    "-DJsonCpp_INCLUDE_DIRS:PATH=${JsonCpp_INCLUDE_DIRS}" )
endif( TubeTK_USE_JsonCpp )

# Include dependent projects, if any.
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

if( NOT DEFINED ${proj}_DIR AND NOT ${USE_SYSTEM_${proj}} )
  set( ${proj}_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj} )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )

  set( ${proj}_QT_ARGS )
  if( TubeTK_USE_QT )
    set( ${proj}_QT_ARGS
      -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
      -DVTK_Group_Qt:BOOL=ON )
  endif( TubeTK_USE_QT )

  set( TubeTK_VTKHDF5_VALGRIND_ARGS )
  if( TubeTK_USE_VALGRIND )
    list( APPEND TubeTK_VTKHDF5_VALGRIND_ARGS
      -DH5_USING_MEMCHECKER:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=ON )
  endif( TubeTK_USE_VALGRIND )

  if( UNIX )
    # Ignore in VTK5 that has limited maintenance.
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-strict-overflow" )
  endif( UNIX )

  if(APPLE)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_ARGS
      -DVTK_USE_CARBON:BOOL=OFF
      -DVTK_USE_COCOA:BOOL=ON
      -DVTK_USE_X:BOOL=OFF
      )
  endif()

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
      -Wno-dev
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      ${CMAKE_CXX_COMPILER_LAUNCHER_FLAG}
      ${CMAKE_C_COMPILER_LAUNCHER_FLAG}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DVTK_LEGACY_REMOVE:BOOL=ON
      ${${proj}_QT_ARGS}
      ${TubeTK_VTKHDF5_VALGRIND_ARGS}
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
