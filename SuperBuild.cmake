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
## Check if sytem ITK or superbuild ITK
##
if( NOT USE_SYSTEM_ITK )

  ##
  ## Insight
  ##
  set( proj Insight )
  ExternalProject_Add( ${proj}
    GIT_REPOSITORY "http://github.com/Slicer/ITK.git"
    GIT_TAG "origin/slicer-4.0"
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
  
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "Insight" )

endif( NOT USE_SYSTEM_ITK )

##
## VTK
##
if( TubeTK_USE_VTK )

  ##
  ## Check if sytem VTK or superbuild VTK
  ##
  if( NOT USE_SYSTEM_VTK )

    if( TubeTK_USE_QT )

      ##
      ## VTK
      ##
      set( proj VTK )
      ExternalProject_Add( VTK
        GIT_REPOSITORY "http://github.com/Slicer/VTK.git"
        GIT_TAG "origin/slicer-4.0"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
        BINARY_DIR VTK-Build
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
        GIT_REPOSITORY "http://github.com/Slicer/VTK.git"
        GIT_TAG "origin/slicer-4.0"
        SOURCE_DIR "${CMAKE_BINARY_DIR}/VTK"
        BINARY_DIR VTK-Build
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
          -DVTK_USE_GUISUPPORT:BOOL=ON
        INSTALL_COMMAND ""
        )
      set( VTK_DIR "${base}/VTK-Build" )

    endif( TubeTK_USE_QT )

    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "VTK" )
    
  endif( NOT USE_SYSTEM_VTK )

endif( TubeTK_USE_VTK )
  
  
##
## TCLAP
##
set( proj tclap )
ExternalProject_Add( ${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/tclap"
  SOURCE_DIR tclap
  BINARY_DIR tclap-Build
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
  INSTALL_COMMAND ""
  )
set( TCLAP_DIR "${base}/tclap-Build" )


##
## ModuleDescriptionParser 
##
set( proj ModuleDescriptionParser )

if( NOT USE_SYSTEM_ITK )
  # Depends on ITK if ITK was build using superbuild
  set( ModuleDescriptionParser_DEPENDS "Insight" )
else( NOT USE_SYSTEM_ITK )
  set( ModuleDescriptionParser_DEPENDS "" )
endif( NOT USE_SYSTEM_ITK )

ExternalProject_Add( ${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/ModuleDescriptionParser"
  SOURCE_DIR ModuleDescriptionParser
  BINARY_DIR ModuleDescriptionParser-Build
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
    -DITK_DIR:PATH=${ITK_DIR}
  INSTALL_COMMAND ""
  DEPENDS ${ModuleDescriptionParser_DEPENDS}
  )
set( ModuleDescriptionParser_DIR "${base}/ModuleDescriptionParser-Build" )


##
## GenerateCLP 
##
set( proj GenerateCLP )
ExternalProject_Add( ${proj}
  SVN_REPOSITORY 
    "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel/GenerateCLP"
  SOURCE_DIR GenerateCLP
  BINARY_DIR GenerateCLP-Build
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
    -DITK_DIR:PATH=${ITK_DIR}
    -DTCLAP_DIR:PATH=${TCLAP_DIR}
    -DModuleDescriptionParser_DIR:PATH=${ModuleDescriptionParser_DIR}
  INSTALL_COMMAND ""
  DEPENDS
    "tclap"
    "ModuleDescriptionParser"
  )
set( GenerateCLP_DIR "${base}/GenerateCLP-Build" )
set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "GenerateCLP" )


##
## OpenIGTLink 
##
if( TubeTK_USE_OpenIGTLink )

  set( proj OpenIGTLink )
  ExternalProject_Add( ${proj}
    SVN_REPOSITORY "http://svn.na-mic.org/NAMICSandBox/trunk/OpenIGTLink"
    SOURCE_DIR OpenIGTLink
    BINARY_DIR OpenIGTLink-Build
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
    INSTALL_COMMAND ""
    )
  set( OpenIGTLink_DIR "${CMAKE_BINARY_DIR}/OpenIGTLink-Build" )
  set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "OpenIGTLink" )

else( TubeTK_USE_OpenIGTLink )

  set( OpenIGTLink_DIR "" )

endif( TubeTK_USE_OpenIGTLink )


if( TubeTK_USE_QT )

  ##
  ## CTK 
  ##
  if( TubeTK_USE_CTK )

    if( TubeTK_USE_VTK )
      if( NOT USE_SYSTEM_VTK )
        set( CTK_DEPENDS "VTK" )
      endif( NOT USE_SYSTEM_VTK )
    else( TubeTK_USE_VTK )
      set( CTK_DEPENDS "" )
    endif( TubeTK_USE_VTK )

    set( proj CTK )
    ExternalProject_Add( CTK
      GIT_REPOSITORY "http://github.com/commontk/CTK.git"
      GIT_TAG "origin/master"
      SOURCE_DIR "${CMAKE_BINARY_DIR}/CTK"
      BINARY_DIR CTK-Build
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
        -DCTK_USE_GIT_PROTOCOL:BOOL=TRUE
        -DCTK_LIB_Widgets:BOOL=ON
        -DCTK_LIB_Visualization/VTK/Widgets:BOOL=OFF
        -DCTK_LIB_PluginFramework:BOOL=OFF
        -DCTK_PLUGIN_org.commontk.eventbus:BOOL=OFF
        -Dgit_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
        -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
        -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
        -DVTK_DIR:PATH=${VTK_DIR}
      INSTALL_COMMAND ""
      DEPENDS
        ${CTK_DEPENDS}
      )
    set( CTK_DIR "${CMAKE_BINARY_DIR}/CTK-Build" )

    set( TubeTK_DEPENDS ${TubeTK_DEPENDS} "CTK" )

  endif( TubeTK_USE_CTK )

endif( TubeTK_USE_QT )


##
## TubeTK - Normal Build
##
set( proj TubeTK )
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
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DBUILD_EXAMPLES:BOOL=${BUILD_EXAMPLES}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DBUILD_DOCUMENTATION:BOOL=${BUILD_DOCUMENTATION}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DTubeTK_USE_SUPERBUILD:BOOL=FALSE
    -DTubeTK_USE_VTK:BOOL=${TubeTK_USE_VTK}
    -DTubeTK_USE_KWSTYLE:BOOL=${TubeTK_USE_KWSTYLE}
    -DTubeTK_USE_CTK:BOOL=${TubeTK_USE_CTK}
    -DTubeTK_USE_FANN:BOOL=${TubeTK_USE_FANN}
    -DTubeTK_USE_QT:BOOL=${TubeTK_USE_QT}
    -DTubeTK_USE_OpenIGTLink:BOOL=${TubeTK_USE_OpenIGTLink}
    -DITK_DIR:PATH=${ITK_DIR}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
    -DOpenIGTLink_DIR:PATH=${OpenIGTLink_DIR}
    -DCTK_DIR:PATH=${CTK_DIR}
    -DQT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
  INSTALL_COMMAND ""
  DEPENDS
    ${TubeTK_DEPENDS}
 )

