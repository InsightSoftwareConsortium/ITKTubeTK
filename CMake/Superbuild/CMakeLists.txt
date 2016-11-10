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

set( base "${CMAKE_BINARY_DIR}" )
set_property( DIRECTORY PROPERTY EP_BASE ${base} )

if( DEFINED BUILD_SHARED_LIBS )
  set( shared ${BUILD_SHARED_LIBS} )
else( DEFINED BUILD_SHARED_LIBS )
  set( shared ON ) # use for BUILD_SHARED_LIBS on all subsequent projects
endif( DEFINED BUILD_SHARED_LIBS )
set( testing OFF ) # use for BUILD_TESTING on all subsequent projects
set( build_type "Release" )
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
set( github_protocol "git" )
if( GIT_PROTOCOL_HTTP )
  set( git_protocol "http" )
  set( github_protocol "https" )
endif( GIT_PROTOCOL_HTTP )

include( ExternalProjectsConfig )

# Use for CMAKE_OSX_* in external projects.
set( CMAKE_OSX_EXTERNAL_PROJECT_ARGS )
if( APPLE )
  list( APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET} )
endif( APPLE )

#
# Add applications and modules build options
#
if( TubeTK_BUILD_APPLICATIONS )
  # Expose the TubeTK_BUILD_* variables to Superbuild
  file( GLOB_RECURSE TubeTK_modules_files
    "${TubeTK_SOURCE_DIR}/TubeTKModules.cmake" )
  foreach( f ${TubeTK_modules_files} )
    ExternalProject_Message( ${PROJECT_NAME} "Including ${f}" )
    include( ${f} )
  endforeach( f ${TubeTK_modules_files} )
endif( TubeTK_BUILD_APPLICATIONS )


macro(list_conditional_append cond list)
  if(${cond})
    list(APPEND ${list} ${ARGN})
  endif()
endmacro()

set( TubeTK_DEPENDENCIES "ITK" "SlicerExecutionModel" )

list_conditional_append( TubeTK_USE_VTK TubeTK_DEPENDENCIES "VTK" )

list_conditional_append( TubeTK_USE_CPPCHECK TubeTK_DEPENDENCIES "Cppcheck" )

list_conditional_append( TubeTK_BUILD_IMAGE_VIEWER TubeTK_DEPENDENCIES
  "ImageViewer" )

list_conditional_append( TubeTK_USE_KWSTYLE TubeTK_DEPENDENCIES "KWStyle" )

list_conditional_append( TubeTK_USE_JsonCpp TubeTK_DEPENDENCIES "JsonCpp" )

list_conditional_append( TubeTK_USE_LIBSVM TubeTK_DEPENDENCIES "LIBSVM" )

list_conditional_append( TubeTK_USE_RandomForest TubeTK_DEPENDENCIES
  "RandomForest" )

list_conditional_append( TubeTK_USE_CTK TubeTK_DEPENDENCIES "CTK" )

set( TubeTK_EXTERNAL_PROJECTS_ARGS )

list_conditional_append( TubeTK_USE_BOOST
  TubeTK_EXTERNAL_PROJECTS_ARGS
    -DBOOST_ROOT=${BOOST_ROOT} )

list_conditional_append( TubeTK_BUILD_USING_SLICER
  TubeTK_EXTERNAL_PROJECTS_ARGS
    -DSlicer_DIR=${Slicer_DIR} )

list_conditional_append( TubeTK_USE_QT
  TubeTK_EXTERNAL_PROJECTS_ARGS
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DTubeTK_REQUIRED_QT_VERSION=${TubeTK_REQUIRED_QT_VERSION} )

list_conditional_append( TubeTK_USE_PYTHON
  TubeTK_EXTERNAL_PROJECTS_ARGS
    -DTubeTK_USE_NUMPY_STACK:BOOL=${TubeTK_USE_NUMPY_STACK}
    -DTubeTK_USE_PYQTGRAPH:BOOL=${TubeTK_USE_PYQTGRAPH}
    -DPythonVirtualEnvDir:PATH=${PythonVirtualEnvDir}
    -DPYTHON_EXECUTABLE:FILEPATH=${PYTHON_EXECUTABLE}
    -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
    -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
 )

# Append TubeTK_BUILD_MODULE_DEFAULT is has been specified by the user
if( DEFINED TubeTK_BUILD_MODULE_DEFAULT )
  list( APPEND TubeTK_EXTERNAL_PROJECTS_ARGS
    -DTubeTK_BUILD_MODULE_DEFAULT:BOOL=${TubeTK_BUILD_MODULE_DEFAULT} )
endif()

# Check dependencies.
set( EXTERNAL_PROJECT_DIR ${TubeTK_SOURCE_DIR}/CMake/Superbuild )

ExternalProject_Include_Dependencies(TubeTK
  DEPENDS_VAR TubeTK_DEPENDENCIES
  SUPERBUILD_VAR TubeTK_USE_SUPERBUILD )

# Normal build.
set( proj TubeTK )

ExternalProject_Add( ${proj}
  ${${proj}_EP_ARGS}
  DEPENDS ${TubeTK_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  DOWNLOAD_DIR ${TubeTK_SOURCE_DIR}
  SOURCE_DIR ${TubeTK_SOURCE_DIR}
  BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build
  INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-build
  CMAKE_GENERATOR ${gen}
  BUILD_AWAYS 1
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
    -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
    -Dgit_protocol:STRING=${git_protocol}
    -Dgithub_protocol:STRING=${github_protocol}
    -DSITE:STRING=${SITE}
    -DTubeTK_BUILD_APPLICATIONS:BOOL=${TubeTK_BUILD_APPLICATIONS}
    -DTubeTK_BUILD_IMAGE_VIEWER:BOOL=${TubeTK_BUILD_IMAGE_VIEWER}
    -DTubeTK_BUILD_USING_SLICER:BOOL=${TubeTK_BUILD_USING_SLICER}
    -DTubeTK_BUILD_WITHIN_SLICER:BOOL=${TubeTK_BUILD_WITHIN_SLICER}
    -DTubeTK_CONFIG_BINARY_DIR:PATH=${TubeTK_CONFIG_BINARY_DIR}
    -DTubeTK_EXECUTABLE_DIRS:STRING=${TubeTK_EXECUTABLE_DIRS}
    -DTubeTK_SUPPORT_2D_IMAGES:BOOL=${TubeTK_SUPPORT_2D_IMAGES}
    -DTubeTK_USE_GPU_ARRAYFIRE:BOOL=${TubeTK_USE_GPU_ARRAYFIRE}
    -DTubeTK_USE_BOOST:BOOL=${TubeTK_USE_BOOST}
    -DTubeTK_USE_CPPCHECK:BOOL=${TubeTK_USE_CPPCHECK}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DTubeTK_USE_EXAMPLES_AS_TESTS:BOOL=${TubeTK_USE_EXAMPLES_AS_TESTS}
    -DTubeTK_USE_JsonCpp:BOOL=${TubeTK_USE_JsonCpp}
    -DTubeTK_USE_LIBSVM:BOOL=${TubeTK_USE_LIBSVM}
    -DTubeTK_USE_RandomForest:BOOL=${TubeTK_USE_RandomForest}
    -DTubeTK_USE_PYTHON:BOOL=${TubeTK_USE_PYTHON}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    -DTubeTK_USE_QT:BOOL=${TubeTK_USE_QT}
    -DTubeTK_USE_SUPERBUILD:BOOL=OFF
    -DTubeTK_USE_VTK:BOOL=${TubeTK_USE_VTK}
    ${TubeTK_EXTERNAL_PROJECTS_ARGS}
  INSTALL_COMMAND "" )
