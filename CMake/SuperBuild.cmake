##############################################################################
#
# Library:   TubeTK
#
# Copyright 2010 Kitware Inc. 28 Corporate Drive,
# Clifton Park, NY, 12065, USA.
#
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 ( the "License" );
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

set( base "${CMAKE_BINARY_DIR}" )
set_property( DIRECTORY PROPERTY EP_BASE ${base} )

set( shared ON ) # use for BUILD_SHARED_LIBS on all subsequent projects
set( testing OFF ) # use for BUILD_TESTING on all subsequent projects
set( build_type "Debug" )
if( CMAKE_BUILD_TYPE )
  set( build_type "${CMAKE_BUILD_TYPE}" )
endif()

set( TubeTK_DEPENDS "" )

set( gen "${CMAKE_GENERATOR}" )

##
## Check if sytem ITK or superbuild ITK (or ITKv4)
##
if( NOT USE_SYSTEM_ITK )

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

  ##
  ## Insight
  ##
  if( TubeTK_USE_ITKV4 )

    set( proj Insight )
    ExternalProject_Add( ${proj}
      GIT_REPOSITORY "${GIT_PROTOCOL}://itk.org/ITK.git"
      GIT_TAG "origin/master"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/Insight"
      BINARY_DIR Insight-Build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        -Dgit_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
        -DCMAKE_BUILD_TYPE:STRING=${build_type}
        -DBUILD_SHARED_LIBS:BOOL=${shared}
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DITK_USE_REVIEW:BOOL=ON
        -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
      INSTALL_COMMAND ""
      )

    set( ITK_DIR "${base}/Insight-Build" )

    # Also get SimpleITK
    set( proj SimpleITK )
    ExternalProject_Add( ${proj}
      GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/SimpleITK/SimpleITK.git"
      GIT_TAG "origin/master"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/SimpleITK"
      BINARY_DIR "SimpleITK-Build"
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
        -DCMAKE_BUILD_TYPE:STRING=${build_type}
        -DITK_DIR:STRING=${CMAKE_BINARY_DIR}/Insight-Build
        -DBUILD_SHARED_LIBS:BOOL=${shared}
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DUSE_TESTING:BOOL=OFF
        -DWRAP_JAVA:BOOL=OFF
        -DWRAP_PYTHON:BOOL=OFF
        -DWRAP_LUA:BOOL=OFF
        -DWRAP_CSHARP:BOOL=OFF
        -DWRAP_TCL:BOOL=OFF
        -DWRAP_R:BOOL=OFF
        -DWRAP_RUBY:BOOL=OFF
        -DUSE_SYSTEM_LUA:BOOL=OFF
      INSTALL_COMMAND ""
      DEPENDS
        "Insight"
      )

    set( SimpleITK_DIR "${base}/SimpleITK-Build" )
    set( TubeTK_SimpleITK_Def "-DSimpleITK_DIR:PATH=${SimpleITK_DIR}" )
    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "SimpleITK" )

  else( TubeTK_USE_ITKV4 )

    set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
    if(APPLE)
      list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
        -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
        -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
        -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
    endif()

    set( proj Insight )
    ExternalProject_Add( ${proj}
      GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/Kitware/ITK.git"
      # release-3.20 branch on 2012-11-27.
      GIT_TAG "787200e5250da02069964063548a39d5e46662e2"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/Insight"
      BINARY_DIR Insight-Build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
        -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
        ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
        -DCMAKE_BUILD_TYPE:STRING=${build_type}
        -DBUILD_SHARED_LIBS:BOOL=${shared}
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DITK_USE_REVIEW:BOOL=ON
        -DITK_USE_REVIEW_STATISTICS:BOOL=ON
        -DITK_USE_OPTIMIZED_REGISTRATION_METHODS:BOOL=ON
        -DITK_USE_CENTERED_PIXEL_COORDINATES_CONSISTENTLY:BOOL=ON
        -DITK_USE_TRANSFORM_IO_FACTORIES:BOOL=ON
        -DITK_LEGACY_REMOVE:BOOL=ON
        -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      INSTALL_COMMAND ""
      )

  endif( TubeTK_USE_ITKV4 )

  set( ITK_DIR "${base}/Insight-Build" )
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "Insight" )

  set( SimpleITK_DIR "" )
  set( TubeTK_SimpleITK_Def "" )

endif( NOT USE_SYSTEM_ITK )

##
## VTK
##
if( TubeTK_USE_VTK )

  ##
  ## Check if sytem VTK or superbuild VTK
  ##
  if( NOT USE_SYSTEM_VTK )

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

    if( TubeTK_USE_QT )

      ##
      ## VTK
      ##
      set( proj VTK )
      ExternalProject_Add( VTK
        GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/Slicer/VTK.git"
        GIT_TAG "3138cd11a9eb76a741ce5e3a86417449acfe7eba"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
        BINARY_DIR VTK-Build
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS
          -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
          -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
          -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
          -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
          -DCMAKE_BUILD_TYPE:STRING=${build_type}
          -DBUILD_SHARED_LIBS:BOOL=${shared}
          -DBUILD_EXAMPLES:BOOL=OFF
          -DBUILD_TESTING:BOOL=OFF
          -DVTK_USE_GUISUPPORT:BOOL=ON
          -DVTK_USE_QVTK_QTOPENGL:BOOL=ON
          -DVTK_USE_QT:BOOL=ON
          -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        INSTALL_COMMAND ""
        )
      set( VTK_DIR "${base}/VTK-Build" )

    else( TubeTK_USE_QT )

      ##
      ## VTK
      ##
      set( proj VTK )
      ExternalProject_Add( VTK
        GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/Slicer/VTK.git"
        GIT_TAG "origin/slicer-4.0"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
        BINARY_DIR VTK-Build
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS
          -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
          -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
          -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
          -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
          -DCMAKE_BUILD_TYPE:STRING=${build_type}
          -DBUILD_SHARED_LIBS:BOOL=${shared}
          -DBUILD_EXAMPLES:BOOL=OFF
          -DBUILD_TESTING:BOOL=OFF
          -DVTK_USE_GUISUPPORT:BOOL=ON
        INSTALL_COMMAND ""
        )
      set( VTK_DIR "${base}/VTK-Build" )

    endif( TubeTK_USE_QT )

    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "VTK" )

  endif( NOT USE_SYSTEM_VTK )

endif( TubeTK_USE_VTK )


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

set( proj SlicerExecutionModel )

# Set dependency list
if( NOT USE_SYSTEM_ITK )
  # Depends on ITK if ITK was build using superbuild
  set(SlicerExecutionModel_DEPENDS "Insight")
else( NOT USE_SYSTEM_ITK )
  set(SlicerExecutionModel_DEPENDS "" )
endif( NOT USE_SYSTEM_ITK )

ExternalProject_Add(${proj}
  GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/Slicer/SlicerExecutionModel.git"
  GIT_TAG "origin/master"
  SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
  BINARY_DIR ${proj}-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}

  INSTALL_COMMAND ""
  DEPENDS
    ${SlicerExecutionModel_DEPENDS}
  )
set( SlicerExecutionModel_DIR ${CMAKE_BINARY_DIR}/${proj}-Build )
set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "SlicerExecutionModel" )

if( TubeTK_USE_QT )

  ##
  ## CTK
  ##
  if( TubeTK_USE_CTK )

    if( NOT USE_SYSTEM_CTK )

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

      if( TubeTK_USE_VTK )
        if( NOT USE_SYSTEM_VTK )
          set( CTK_DEPENDS "VTK" )
        endif( NOT USE_SYSTEM_VTK )
      else( TubeTK_USE_VTK )
        set( CTK_DEPENDS "" )
      endif( TubeTK_USE_VTK )

      set( proj CTK )
      ExternalProject_Add( CTK
        GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/commontk/CTK.git"
        GIT_TAG "11a59a395deb8ff50ddb07b427e47e6665447586"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/CTK"
        BINARY_DIR CTK-Build
        CMAKE_GENERATOR ${gen}
        CMAKE_ARGS
          -DCMAKE_BUILD_TYPE:STRING=${build_type}
          -DBUILD_SHARED_LIBS:BOOL=${shared}
          -DBUILD_TESTING:BOOL=OFF
          -DCTK_USE_GIT_PROTOCOL:BOOL=TRUE
          -DCTK_LIB_Widgets:BOOL=ON
          -DCTK_LIB_Visualization/VTK/Widgets:BOOL=OFF
          -DCTK_LIB_PluginFramework:BOOL=OFF
          -DCTK_PLUGIN_org.commontk.eventbus:BOOL=OFF
          -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        INSTALL_COMMAND ""
        DEPENDS
          ${CTK_DEPENDS}
        )
      set( CTK_DIR "${CMAKE_BINARY_DIR}/CTK-Build" )

      set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "CTK" )

    endif( NOT USE_SYSTEM_CTK )

  endif( TubeTK_USE_CTK )

endif( TubeTK_USE_QT )


##
## TubeTK - Normal Build
##
set( proj TubeTK )
if( TubeTK_USE_KWSTYLE )
  set( kwstyle_dashboard_submission_arg
    "-DKWSTYLE_DASHBOARD_SUBMISSION:BOOL=${KWSTYLE_DASHBOARD_SUBMISSION}" )
endif()
ExternalProject_Add( ${proj}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  BINARY_DIR TubeTK-Build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    -DBUILDNAME:STRING=${BUILDNAME}
    -DSITE:STRING=${SITE}
    -DMAKECOMMAND:STRING=${MAKECOMMAND}
    -DCMAKE_BUILD_TYPE:STRING=${build_type}
    -DBUILD_SHARED_LIBS:BOOL=${shared}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DTubeTK_USE_VTK:BOOL=${TubeTK_USE_VTK}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    ${kwstyle_dashboard_submission_arg}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DTubeTK_USE_QT:BOOL=${TubeTK_USE_QT}
    -DTubeTK_USE_ITKV4:BOOL=${TubeTK_USE_ITKV4}
    -DTubeTK_EXECUTABLE_DIRS:BOOL=${TubeTK_EXECUTABLE_DIRS}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
    -DCTK_DIR:PATH=${CTK_DIR}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    ${TubeTK_SimpleITK_Def}
  INSTALL_COMMAND ""
  DEPENDS
    ${TubeTK_DEPENDS}
  )
