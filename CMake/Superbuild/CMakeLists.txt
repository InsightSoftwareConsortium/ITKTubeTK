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

include( ExternalProject )
include( SlicerMacroCheckExternalProjectDependency )

set( base "${CMAKE_BINARY_DIR}" )
set_property( DIRECTORY PROPERTY EP_BASE ${base} )

if( DEFINED BUILD_SHARED_LIBS )
  set( shared ${BUILD_SHARED_LIBS} )
else( DEFINED BUILD_SHARED_LIBS )
  set( shared ON ) # use for BUILD_SHARED_LIBS on all subsequent projects
endif( DEFINED BUILD_SHARED_LIBS )
set( testing OFF ) # use for BUILD_TESTING on all subsequent projects
set( build_type "Debug" )
if( CMAKE_BUILD_TYPE )
  set( build_type "${CMAKE_BUILD_TYPE}" )
endif( CMAKE_BUILD_TYPE )

set( gen "${CMAKE_GENERATOR}" )

##
## Find GIT and determine proper protocol for accessing GIT repos.
##  - Users may need to choose HTTP is they are behind a firewall.
##
if( NOT GIT_EXECUTABLE )
  find_package( Git REQUIRED )
endif( NOT GIT_EXECUTABLE )

option( GIT_PROTOCOL_HTTP
  "Use HTTP for git access (useful if behind a firewall)" OFF )
if( GIT_PROTOCOL_HTTP )
  set( GIT_PROTOCOL "http" CACHE STRING "Git protocol for file transfer" )
else( GIT_PROTOCOL_HTTP )
  set( GIT_PROTOCOL "git" CACHE STRING "Git protocol for file transfer" )
endif( GIT_PROTOCOL_HTTP )
mark_as_advanced( GIT_PROTOCOL )

#
# Use for CMAKE_OSX_* in external projects.
#
set( CMAKE_OSX_EXTERNAL_PROJECT_ARGS )
if( APPLE )
  list( APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET} )
endif( APPLE )

## Set TubeTK Dependencies
set( TubeTK_DEPENDENCIES "" )

## Check if system JsonCpp or superbuild JsonCpp
if( NOT USE_SYSTEM_JsonCpp )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "JsonCpp" )
else( NOT USE_SYSTEM_JsonCpp )
  find_package( JsonCpp REQUIRED )
endif( NOT USE_SYSTEM_JsonCpp )

## Check if system ITK or superbuild ITK
if( NOT USE_SYSTEM_ITK )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "Insight" )

  if( TubeTK_USE_SimpleITK )
    set(TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "SimpleITK" )
  endif( TubeTK_USE_SimpleITK )
endif( NOT USE_SYSTEM_ITK )

## Check if system VTK or superbuild VTK
if( TubeTK_USE_VTK AND NOT USE_SYSTEM_VTK )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "VTK" )
endif( TubeTK_USE_VTK AND NOT USE_SYSTEM_VTK )

## Check if system ParameterSerializer or superbuild ParameterSerializer
if( NOT USE_SYSTEM_ParameterSerializer )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "ParameterSerializer" )
else( NOT USE_SYSTEM_ParameterSerializer )
  find_package( ParameterSerializer REQUIRED )
endif( NOT USE_SYSTEM_ParameterSerializer )

## SlicerExecutionModel
if( NOT USE_SYSTEM_SlicerExecutionModel )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "SlicerExecutionModel" )
endif( NOT USE_SYSTEM_SlicerExecutionModel )

if( TubeTK_USE_QT )
  ## CTK
  if( TubeTK_USE_CTK AND NOT USE_SYSTEM_CTK )
      set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "CTK" )
  endif( TubeTK_USE_CTK AND NOT USE_SYSTEM_CTK )
endif( TubeTK_USE_QT )

if( TubeTK_BUILD_ImageViewer )
  include( "External_ImageViewer.cmake" )
endif( TubeTK_BUILD_ImageViewer )

## LibSVM
##
if( TubeTK_USE_LIBSVM )
 include( "External_LibSVM.cmake" )
endif( TubeTK_USE_LIBSVM )

##
## TubeTK - Check Dependencies
##

set( EXTERNAL_PROJECT_DIR ${TubeTK_SOURCE_DIR}/CMake/Superbuild )
SlicerMacroCheckExternalProjectDependency( TubeTK )

##
## TubeTK - Normal Build
##

set( proj TubeTK )
if( TubeTK_USE_KWSTYLE )
  set( kwstyle_dashboard_submission_arg
    "-DKWSTYLE_DASHBOARD_SUBMISSION:BOOL=${KWSTYLE_DASHBOARD_SUBMISSION}" )
endif( TubeTK_USE_KWSTYLE )

ExternalProject_Add( ${proj}
  DEPENDS ${TubeTK_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR "${TubeTK_SOURCE_DIR}"
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build
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
    -DBUILDNAME:STRING=${BUILDNAME}
    -DSITE:STRING=${SITE}
    -DMAKECOMMAND:STRING=${MAKECOMMAND}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DTubeTK_CONFIG_BINARY_DIR:PATH=${TubeTK_CONFIG_BINARY_DIR}
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DTubeTK_BUILD_SLICER_EXTENSION:BOOL=${TubeTK_BUILD_SLICER_EXTENSION}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    ${kwstyle_dashboard_submission_arg}
    -DTubeTK_USE_Slicer:BOOL=${TubeTK_USE_Slicer}
    -DSlicer_DIR:PATH=${Slicer_DIR}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DCTK_DIR:PATH=${CTK_DIR}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
    -DTubeTK_USE_VTK:BOOL=${TubeTK_USE_VTK}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DTubeTK_USE_QT:BOOL=${TubeTK_USE_QT}
    -DTubeTK_USE_Boost:BOOL=${TubeTK_USE_Boost}
    -DTubeTK_USE_LIBSVM:BOOL=${TubeTK_USE_LIBSVM}
    -DJsonCpp_DIR:PATH=${JsonCpp_DIR}
    -DParameterSerializer_DIR:PATH=${ParameterSerializer_DIR}
    -DTubeTK_EXECUTABLE_DIRS:STRING=${TubeTK_EXECUTABLE_DIRS}
    -DTubeTK_REQUIRED_QT_VERSION=${TubeTK_REQUIRED_QT_VERSION}
    -DTubeTK_USE_CPPCHECK:BOOL=${TubeTK_USE_CPPCHECK}
    ${TubeTK_SimpleITK_Def}
    -DTubeTK_BUILD_ALL_MODULES=${TubeTK_BUILD_ALL_MODULES}
  INSTALL_COMMAND "" )
