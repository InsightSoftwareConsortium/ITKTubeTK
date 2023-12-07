##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#       https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
##############################################################################
cmake_minimum_required(VERSION 3.16.3)

if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

project( TubeTK )

set( TubeTK_VERSION 1.3.5 )

####
# Setup CMake
####
list( APPEND CMAKE_MODULE_PATH ${TubeTK_SOURCE_DIR}/CMake )

set( TubeTK_LIBRARIES TubeTK )

####
# Setup include directories
####
set( TubeTK_INCLUDE_DIRS
  $<INSTALL_INTERFACE:include>
  $<INSTALL_INTERFACE:include/Common>
  $<INSTALL_INTERFACE:include/Filtering>
  $<INSTALL_INTERFACE:include/IO>
  $<INSTALL_INTERFACE:include/MetaIO>
  $<INSTALL_INTERFACE:include/Numerics>
  $<INSTALL_INTERFACE:include/ObjectDocuments>
  $<INSTALL_INTERFACE:include/Registration>
  $<INSTALL_INTERFACE:include/Segmentation> )

####
# Build base library and optional components
####
set( TubeTK_SYSTEM_INCLUDE_DIRS
  ${TubeTK_BINARY_DIR}/include
  ${TubeTK_SOURCE_DIR}/src
  ${TubeTK_SOURCE_DIR}/src/Common
  ${TubeTK_SOURCE_DIR}/src/Filtering
  ${TubeTK_SOURCE_DIR}/src/IO
  ${TubeTK_SOURCE_DIR}/src/MetaIO
  ${TubeTK_SOURCE_DIR}/src/Numerics
  ${TubeTK_SOURCE_DIR}/src/ObjectDocuments
  ${TubeTK_SOURCE_DIR}/src/Registration
  ${TubeTK_SOURCE_DIR}/src/Segmentation )

include_directories( ${TubeTK_SYSTEM_INCLUDE_DIRS} )

####
# Build TubeTK module
####
if( NOT ITK_SOURCE_DIR )
  find_package( ITK REQUIRED )
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )
  include( ITKModuleExternal )
  #if( WASI OR EMSCRIPTEN )
    #set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
    #add_subdirectory( wasm )
  #endif()
else()
  set(ITK_DIR ${CMAKE_BINARY_DIR})
  itk_module_impl()
endif()

itk_module_examples()
