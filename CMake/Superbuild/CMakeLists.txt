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

set( CMAKE_MODULE_PATH ${TubeTK_SOURCE_DIR}/CMake/Superbuild ${CMAKE_MODULE_PATH} )

include( ExternalProject )
include( TubeTKMacroCheckExternalProjectDependency )

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

# Find Git and determine proper protocol for accessing Git repositories. Users
# may need to choose HTTP is they are behind a firewall.
if( NOT GIT_EXECUTABLE )
  find_package( Git REQUIRED )
endif( NOT GIT_EXECUTABLE )

option( GIT_PROTOCOL_HTTP
  "Use HTTP for git access (useful if behind a firewall)" OFF )
mark_as_advanced( GIT_PROTOCOL_HTTP )

set( git_protocol "git" )
if( GIT_PROTOCOL_HTTP )
  set( git_protocol "http" )
endif( GIT_PROTOCOL_HTTP )

option( SVN_PROTOCOL_HTTP
  "Use HTTP for svn access (useful if behind a firewall)" OFF )
mark_as_advanced( SVN_PROTOCOL_HTTP )

set( svn_protocol "svn" )
if( SVN_PROTOCOL_HTTP )
  set( svn_protocol "http" )
endif( SVN_PROTOCOL_HTTP )

include( ExternalProjectsConfig )

# Use for CMAKE_OSX_* in external projects.
set( CMAKE_OSX_EXTERNAL_PROJECT_ARGS )
if( APPLE )
  list( APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET} )
endif( APPLE )

# Set dependencies.
set( TubeTK_DEPENDENCIES "ITK" "JsonCpp" "ParameterSerializer"
  "SlicerExecutionModel" )

if( TubeTK_USE_SIMPLEITK )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "SimpleITK" )
endif( TubeTK_USE_SIMPLEITK )

if( TubeTK_USE_VTK )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "VTK" )
endif( TubeTK_USE_VTK )
if( TubeTK_USE_CPPCHECK )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "Cppcheck" )
endif( TubeTK_USE_CPPCHECK )

set( TubeTK_QT_VARIABLES )
if( TubeTK_USE_QT )
  # CTK.
  if( TubeTK_USE_CTK )
      set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "CTK" )
  endif( TubeTK_USE_CTK )
  set( TubeTK_QT_VARIABLES
    -DTubeTK_REQUIRED_QT_VERSION=${TubeTK_REQUIRED_QT_VERSION} )
endif( TubeTK_USE_QT )

set( TubeTK_SLICER_DIRS
    -DTubeTK_SLICER_SOURCE_DIR:PATH=${TubeTK_SLICER_SOURCE_DIR} )
if( TubeTK_USE_SLICER )
  set( TubeTK_SLICER_DIRS
    ${TubeTK_SLICER_DIRS}
    -DSlicer_DIR:PATH=${Slicer_DIR} )
endif( TubeTK_USE_SLICER )
if( TubeTK_USE_KWSTYLE )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "KWStyle" )
endif( TubeTK_USE_KWSTYLE )

if( TubeTK_USE_LIBSVM )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "LIBSVM" )
endif( TubeTK_USE_LIBSVM )

if( TubeTK_USE_IMAGE_VIEWER )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "ImageViewer" )
endif( TubeTK_USE_IMAGE_VIEWER )

set( TubeTK_EXTERNAL_PROJECTS_ARGS )

if( TubeTK_USE_BOOST )
  list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS -DBOOST_ROOT=${BOOST_ROOT} )
endif( TubeTK_USE_BOOST )

if( TubeTK_USE_KWSTYLE )
  set( TubeTK_DEPENDENCIES ${TubeTK_DEPENDENCIES} "KWStyle" )
endif( TubeTK_USE_KWSTYLE )

if( TubeTK_USE_PYTHON )
  list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS
    -DTubeTK_USE_PYQTGRAPH:BOOL=${TubeTK_USE_PYQTGRAPH}
    -DTubeTK_USE_NUMPY:BOOL=${TubeTK_USE_NUMPY}
    -DTubeTK_USE_NOTEBOOKS:BOOL=${TubeTK_USE_NOTEBOOKS} )
endif( TubeTK_USE_PYTHON )

# Check dependencies.
set( EXTERNAL_PROJECT_DIR ${TubeTK_SOURCE_DIR}/CMake/Superbuild )
TubeTKMacroCheckExternalProjectDependency( TubeTK )


# Normal build.
set( proj TubeTK )

ExternalProject_Add( ${proj}
  DEPENDS ${TubeTK_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  DOWNLOAD_DIR ${TubeTK_SOURCE_DIR}
  SOURCE_DIR ${TubeTK_SOURCE_DIR}
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build
  INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    -DBUILDNAME:STRING=${BUILDNAME}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DMAKECOMMAND:STRING=${MAKECOMMAND}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DSITE:STRING=${SITE}
    -DTubeTK_BUILD_ALL_MODULES=${TubeTK_BUILD_ALL_MODULES}
    -DTubeTK_BUILD_MODULES=${TubeTK_BUILD_MODULES}
    -DTubeTK_BUILD_SLICER_EXTENSION:BOOL=${TubeTK_BUILD_SLICER_EXTENSION}
    -DTubeTK_CONFIG_BINARY_DIR:PATH=${TubeTK_CONFIG_BINARY_DIR}
    -DTubeTK_EXECUTABLE_DIRS:STRING=${TubeTK_EXECUTABLE_DIRS}
    -DTubeTK_USE_PYTHON:BOOL=${TubeTK_USE_PYTHON}
    -DTubeTK_USE_BOOST:BOOL=${TubeTK_USE_BOOST}
    -DTubeTK_USE_CPPCHECK:BOOL=${TubeTK_USE_CPPCHECK}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DTubeTK_USE_IMAGE_VIEWER:BOOL=${TubeTK_USE_IMAGE_VIEWER}
    -DTubeTK_USE_LIBSVM:BOOL=${TubeTK_USE_LIBSVM}
    -DTubeTK_USE_QT:BOOL=${TubeTK_USE_QT}
    -DTubeTK_USE_SLICER:BOOL=${TubeTK_USE_SLICER}
    -DTubeTK_USE_SUPERBUILD:BOOL=OFF
    -DTubeTK_USE_VTK:BOOL=${TubeTK_USE_VTK}
    ${TubeTK_EXTERNAL_PROJECTS_ARGS}
    ${TubeTK_QT_VARIABLES}
    ${TubeTK_SLICER_DIRS}
  INSTALL_COMMAND "" )
