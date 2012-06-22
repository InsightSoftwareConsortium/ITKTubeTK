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
set( GITCOMMAND "${SITE_GIT_COMMAND}" )
set( SVNCOMMAND "${SITE_SVN_COMMAND}" )
set( GIT_EXECUTABLE "${SITE_GIT_COMMAND}" )
set( QT_QMAKE_EXECUTABLE "${SITE_QMAKE_COMMAND}" }

set( CMAKE_GENERATOR ${SITE_CMAKE_GENERATOR} )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SITE_CXX_FLAGS}" )
set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SITE_C_FLAGS}" )
set( CMAKE_EXE_LINKER_FLAGS
  "${CMAKE_EXE_LINKER_FLAGS} ${SITE_EXE_LINKER_FLAGS}" )
set( CMAKE_SHARED_LINKER_FLAGS
  "${CMAKE_SHARED_LINKER_FLAGS} ${SITE_SHARED_LINKER_FLAGS}" )

message( "Running ${SCRIPT_NAME} script..." )

file( WRITE "${TUBETK_BINARY_DIR}/CMakeCache.txt" "
  SITE:STRING=${CTEST_SITE}
  BUILDNAME:STRING=${CTEST_BUILD_NAME}
  BUILD_TESTING:BOOL=ON
  BUILD_SHARED_LIBS:BOOL=ON
  CMAKE_PROJECT_NAME:STRING=${CMAKE_PROJECT_NAME}
  CMAKE_GENERATOR:STRING=${SITE_CMAKE_GENERATOR}
  CMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}
  CMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  CMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  CMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  CMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  MAKECOMMAND:FILEPATH=${CTEST_BUILD_COMMAND}
  GITCOMMAND:FILEPATH=${GITCOMMAND}
  SVNCOMMAND:FILEPATH=${SVNCOMMAND}
  GIT_EXECUTABLE:FILEPATH=${GITCOMMAND}
  QT_QMAKE_EXECUTABLE:FILEPATH=${QT_QMAKE_EXECUTABLE}
  TubeTK_USE_SUPERBUILD:BOOL=ON
  TubeTK_USE_KWSTYLE:BOOL=ON
  " )

if( EXISTS ${TUBETK_BINARY_DIR}/Insight-Build/CMakeCache.txt )
  file( APPEND "${TUBETK_BINARY_DIR}/CMakeCache.txt"
    "
    ITK_DIR:PATH=${TUBETK_BINARY_DIR}/Insight-Build
    " )
endif()

if( EXISTS ${TUBETK_BINARY_DIR}/GenerateCLP-Build/CMakeCache.txt)
  file( APPEND "${TUBETK_BINARY_DIR}/CMakeCache.txt"
    "
    GenerateCLP_DIR:PATH=${TUBETK_BINARY_DIR}/GenerateCLP-Build
    " )
endif()

if( EXISTS ${TUBETK_BINARY_DIR}/SlicerExecutionModel-Build/CMakeCache.txt)
  file( APPEND "${TUBETK_BINARY_DIR}/CMakeCache.txt"
    "
    SlicerExecutionModel_DIR:PATH=${TUBETK_BINARY_DIR}/SlicerExecutionModel-Build
    " )
endif()

if( EXISTS ${TUBETK_BINARY_DIR}/CTK-Build/CMakeCache.txt )
  file( APPEND "${TUBETK_BINARY_DIR}/CMakeCache.txt"
    "
    CTK_DIR:PATH=${TUBETK_BINARY_DIR}/CTK-Build
    " )
endif()

if( EXISTS ${TUBETK_BINARY_DIR}/VTK-Build/CMakeCache.txt )
  file( APPEND "${TUBETK_BINARY_DIR}/CMakeCache.txt"
    "
    VTK_DIR:PATH=${TUBETK_BINARY_DIR}/VTK-Build
    " )
endif()
