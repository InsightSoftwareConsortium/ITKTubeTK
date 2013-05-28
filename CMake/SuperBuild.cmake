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

set( TubeTK_DEPENDS "" )

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
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET}
  )
endif( APPLE )

##
## Check if system JsonCpp or superbuild JsonCpp
##
if( NOT USE_SYSTEM_JsonCpp )
  ##
  ## JsonCpp
  ##
  set( proj JsonCpp )
  ExternalProject_Add( JsonCpp
    GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/jsoncpp-cmake.git"
    GIT_TAG "838583b52ba3fac143ec4e09bd189498a7ff9686"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/JsonCpp"
    BINARY_DIR JsonCpp-Build
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
    INSTALL_COMMAND ""
    )
  set( JsonCpp_DIR "${base}/JsonCpp-Build" )
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "JsonCpp" )
endif( NOT USE_SYSTEM_JsonCpp )

if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  ##
  ## Check if system ITK or superbuild ITK
  ##
  if( NOT USE_SYSTEM_ITK )
    ##
    ## Insight
    ##
    set( proj Insight )
    ExternalProject_Add( ${proj}
      GIT_REPOSITORY "${GIT_PROTOCOL}://itk.org/ITK.git"
      # master on 2013-02-28
      GIT_TAG "9fe3eaaac2321f"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/Insight"
      BINARY_DIR Insight-Build
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
        -DITK_USE_REVIEW:BOOL=ON
        -DITKV3_COMPATIBILITY:BOOL=ON
      INSTALL_COMMAND ""
      )

    set( ITK_DIR "${base}/Insight-Build" )

    if( TubeTK_USE_SimpleITK )
      # Also get SimpleITK
      set( proj SimpleITK )
      ExternalProject_Add( ${proj}
        GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/SimpleITK/SimpleITK.git"
        GIT_TAG "v0.6.1"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/SimpleITK"
        BINARY_DIR "SimpleITK-Build"
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
          -DITK_DIR:STRING=${CMAKE_BINARY_DIR}/Insight-Build
          -DBUILD_SHARED_LIBS:BOOL=${shared}
          -DBUILD_TESTING:BOOL=OFF
          -DWRAP_JAVA:BOOL=OFF
          -DWRAP_PYTHON:BOOL=OFF
          -DWRAP_LUA:BOOL=OFF
          -DWRAP_CSHARP:BOOL=OFF
          -DWRAP_TCL:BOOL=OFF
          -DWRAP_R:BOOL=OFF
          -DWRAP_RUBY:BOOL=OFF
        INSTALL_COMMAND ""
        DEPENDS
          "Insight"
        )

      set( SimpleITK_DIR "${base}/SimpleITK-Build" )
      set( TubeTK_SimpleITK_Def "-DSimpleITK_DIR:PATH=${SimpleITK_DIR}" )
    endif( TubeTK_USE_SimpleITK )
    set( SimpleITK_DIR "" )
    set( TubeTK_SimpleITK_Def "" )

    set( ITK_DIR "${base}/Insight-Build" )
    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "Insight" )
  endif( NOT USE_SYSTEM_ITK )

  ##
  ## VTK
  ##
  if( TubeTK_USE_VTK )

    ##
    ## Check if system VTK or superbuild VTK
    ##
    if( NOT USE_SYSTEM_VTK )

      if( TubeTK_USE_QT )

        ##
        ## VTK
        ##
        set( proj VTK )
        ExternalProject_Add( VTK
          GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/Slicer/VTK.git"
          GIT_TAG "03ddda8cd503b957268fee138c7b57465e5842a4"
          SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
          BINARY_DIR VTK-Build
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
          GIT_TAG "03ddda8cd503b957268fee138c7b57465e5842a4"
          SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
          BINARY_DIR VTK-Build
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
            -DVTK_USE_GUISUPPORT:BOOL=ON
          INSTALL_COMMAND ""
          )
        set( VTK_DIR "${base}/VTK-Build" )

      endif( TubeTK_USE_QT )

      set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "VTK" )

    endif( NOT USE_SYSTEM_VTK )

  endif( TubeTK_USE_VTK )

  ##
  ## Check if system ParameterSerializer or superbuild ParameterSerializer
  ##
  if( NOT USE_SYSTEM_ParameterSerializer )
    ##
    ## ParameterSerializer
    ##
    set( proj ParameterSerializer )
    # Set dependency list
    if( NOT USE_SYSTEM_ITK )
      # Depends on ITK if ITK was build using superbuild
      set( ${proj}_DEPENDS "Insight")
    else()
      set( ${proj}_DEPENDS "" )
    endif()
    if( NOT USE_SYSTEM_JsonCpp )
      # Depends on ITK if ITK was build using superbuild
      set( ${proj}_DEPENDS ${${proj}_DEPENDS} "JsonCpp")
    else()
      set( ${proj}_DEPENDS "" )
    endif()
    ExternalProject_Add( ${proj}
      GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/ParameterSerializer.git"
      GIT_TAG "23911582d78ced102946c95dd67294d4057bb645"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/ParameterSerializer"
      BINARY_DIR ParameterSerializer-Build
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
        ${${proj}_DEPENDS}
      )
    set( ParameterSerializer_DIR "${base}/ParameterSerializer-Build" )
    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "ParameterSerializer" )
  endif( NOT USE_SYSTEM_ParameterSerializer )

  set( proj SlicerExecutionModel )
  # Set dependency list
  if( NOT USE_SYSTEM_ITK )
    # Depends on ITK if ITK was build using superbuild
    set( SlicerExecutionModel_DEPENDS "Insight")
  else()
    set( SlicerExecutionModel_DEPENDS "" )
  endif()
  if( NOT USE_SYSTEM_ParameterSerializer )
    # Depends on ITK if ITK was build using superbuild
    set( SlicerExecutionModel_DEPENDS "ParameterSerializer")
  else()
    set( SlicerExecutionModel_DEPENDS "" )
  endif()
  ExternalProject_Add( ${proj}
    # Has enhancements required for the SEM ParameterSerializer.
    # Change back to upstream after pull request merged.
    GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/SlicerExecutionModel.git"
    GIT_TAG "8ebf568a7b90a1d0e5ec7cf4cf2e8e490a24d0ab"
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    BINARY_DIR ${proj}-Build
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
      -DBUILD_TESTING:BOOL=OFF
      -DITK_DIR:PATH=${ITK_DIR}
      -DSlicerExecutionModel_USE_SERIALIZER:BOOL=ON
      -DJsonCpp_DIR:PATH=${JsonCpp_DIR}
      -DParameterSerializer_DIR:PATH=${ParameterSerializer_DIR}
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
          GIT_TAG "1e823001cb41c92667299635643f1007876d09f6"
          SOURCE_DIR "${CMAKE_BINARY_DIR}/CTK"
          BINARY_DIR CTK-Build
          CMAKE_GENERATOR ${gen}
          CMAKE_ARGS
            -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
            -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
            -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
            -DCMAKE_BUILD_TYPE:STRING=${build_type}
            ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
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

else( NOT TubeTK_BUILD_SLICER_EXTENSION )

  unset( CTK_DIR )
  unset( ITK_DIR )
  unset( VTK_DIR )
  unset( QT_QMAKE_EXECUTABLE )
  unset( SlicerExecutionModel_DIR )

endif( NOT TubeTK_BUILD_SLICER_EXTENSION )

if( TubeTK_BUILD_ImageViewer )
  ##
  ## A convenient 2D/3D image viewer that can handle anisotropic spacing.
  ##
  set( ImageViewer_DEPENDS )
  if( NOT TubeTK_BUILD_SLICER_EXTENSION )
    if( NOT USE_SYSTEM_ITK )
      set( ImageViewer_DEPENDS Insight )
    endif( NOT USE_SYSTEM_ITK )
  else( NOT TubeTK_BUILD_SLICER_EXTENSION )
  endif( NOT TubeTK_BUILD_SLICER_EXTENSION )
  set( proj ImageViewer )
  ExternalProject_Add( ImageViewer
    GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/ImageViewer.git"
    GIT_TAG "b2ac8ef25da64f17b1e2e08ddb5415ec95783515"
    SOURCE_DIR "${CMAKE_BINARY_DIR}/ImageViewer"
    BINARY_DIR ImageViewer-Build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DITK_DIR:PATH=${ITK_DIR}
    INSTALL_COMMAND ""
    DEPENDS
      ${ImageViewer_DEPENDS}
    )
endif( TubeTK_BUILD_ImageViewer )

## LibSVM
##
if( TubeTK_USE_LIBSVM )
  set( proj LIBSVM )
  ExternalProject_Add( ${proj}
    SOURCE_DIR "${CMAKE_BINARY_DIR}/${proj}"
    BINARY_DIR ${proj}-Build
    GIT_REPOSITORY "${GIT_PROTOCOL}://github.com/TubeTK/cmake-libsvm.git"
    GIT_TAG "a802a4224a6c3d7458e46887e77d75bf305a105b"
    CMAKE_ARGS
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
      -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
      -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
      -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_BINARY_DIR}/${proj}-Build
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      -DBUILD_SHARED_LIBS:BOOL=${shared}
  )
  set( LIBSVM_DIR "${CMAKE_BINARY_DIR}/${proj}-Build" )
endif( TubeTK_USE_LIBSVM )

##
## TubeTK - Normal Build
##
set( proj TubeTK )

if( TubeTK_USE_KWSTYLE )
  set( kwstyle_dashboard_submission_arg
    "-DKWSTYLE_DASHBOARD_SUBMISSION:BOOL=${KWSTYLE_DASHBOARD_SUBMISSION}" )
endif( TubeTK_USE_KWSTYLE )

set( TubeTK_cmake_args )
if( NOT TubeTK_BUILD_SLICER_EXTENSION )
  list( APPEND TubeTK_cmake_args
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DCTK_DIR:PATH=${CTK_DIR}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
    -DSlicerExecutionModel_DIR:PATH=${SlicerExecutionModel_DIR}
    )
endif( NOT TubeTK_BUILD_SLICER_EXTENSION )

ExternalProject_Add( ${proj}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}"
  BINARY_DIR TubeTK-Build
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
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DTubeTK_BUILD_SLICER_EXTENSION:BOOL=${TubeTK_BUILD_SLICER_EXTENSION}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    ${kwstyle_dashboard_submission_arg}
    ${TubeTK_cmake_args}
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
  INSTALL_COMMAND ""
  DEPENDS
    ${TubeTK_DEPENDS}
  )
