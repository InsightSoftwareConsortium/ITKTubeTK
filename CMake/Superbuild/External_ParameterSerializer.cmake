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

#
# ParameterSerializer
#

set( proj ParameterSerializer )

# Make sure this file is included only once
get_filename_component( CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE )
if( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED )
  return()
endif( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED )
set( ${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1 )

# Sanity checks
if( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )
  message( FATAL_ERROR "${proj}_DIR variable is defined but corresponds to non-existing directory" )
endif( DEFINED ${proj}_DIR AND NOT EXISTS ${${proj}_DIR} )

# Set dependency list
set( ${proj}_DEPENDENCIES "" )
if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  if( NOT USE_SYSTEM_ITK )
    # Depends on ITK if ITK was build using superbuild
    set( ${proj}_DEPENDENCIES ${${proj}_DEPENDENCIES} "Insight")
  endif( NOT USE_SYSTEM_ITK )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )
if( NOT USE_SYSTEM_JsonCpp )
  set( ${proj}_DEPENDENCIES ${${proj}_DEPENDENCIES} "JsonCpp")
endif( NOT USE_SYSTEM_JsonCpp )

# Include dependent projects if any
SlicerMacroCheckExternalProjectDependency( ${proj} )

if( NOT DEFINED ${proj}_DIR )
  set( ${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build )

  ExternalProject_Add( ${proj}
    GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/ParameterSerializer.git"
    GIT_TAG "23911582d78ced102946c95dd67294d4057bb645"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/${proj}"
    BINARY_DIR ${${proj}_DIR}
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
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DITK_DIR:PATH=${ITK_DIR}
      -DJsonCpp_DIR:PATH=${JsonCpp_DIR}
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES} )

else( NOT DEFINED ${proj}_DIR )
  # The project is provided using ${proj}_DIR, nevertheless since other project may depend on ${proj},
  # let's add an 'empty' one
  SlicerMacroEmptyExternalProject( ${proj} "${${proj}_DEPENDENCIES}" )
endif( NOT DEFINED ${proj}_DIR )
