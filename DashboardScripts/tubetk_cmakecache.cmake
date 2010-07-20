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

set( CTEST_SITE "${SITE_NAME}" )
set( CTEST_BUILD_NAME 
  "${SITE_BUILD_NAME}-${SCRIPT_NAME}-${RUN_MODEL}" )
set( CTEST_BUILD_CONFIGURATION "${SITE_BUILD_TYPE}" )

set( CTEST_TEST_TIMEOUT 1500 )

set( CTEST_CMAKE_GENERATOR ${SITE_CMAKE_GENERATOR} )
set( CTEST_NOTES_FILES "${SITE_SCRIPT_DIR}/${CTEST_SCRIPT_NAME}" )
set( CTEST_SOURCE_DIRECTORY "${SITE_SOURCE_DIR}" )
set( CTEST_BINARY_DIRECTORY "${SITE_BINARY_DIR}/${SITE_BINARY_SUBDIR}" )

set( CTEST_BUILD_COMMAND "${SITE_MAKE_COMMAND}" )
set( CTEST_CMAKE_COMMAND "${SITE_CMAKE_COMMAND}" )
set( CTEST_CTEST_COMMAND "${SITE_CTEST_COMMAND}" )
set( CTEST_QMAKE_COMMAND "${SITE_QMAKE_COMMAND}" )

set( GITCOMMAND "${SITE_GIT_COMMAND}" )
set( CVSCOMMAND "${SITE_CVS_COMMAND}" )
set( SVNCOMMAND "${SITE_SVN_COMMAND}" )

set( CTEST_UPDATE_COMMAND "${SITE_GIT_COMMAND}" )
set( CTEST_UPDATE_OPTIONS "${GIT_UPDATE_OPTIONS}" )

message( "Running ${SCRIPT_NAME} script..." )

# force a build if this is the first run and the build dir is empty
if( NOT EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" )

  message( "First time build!" )

  file( WRITE "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" "
    SITE:STRING=${CTEST_SITE}
    BUILDNAME:STRING=${CTEST_BUILD_NAME}
    BUILD_TESTING:BOOL=ON
    BUILD_SHARED_LIBS:BOOL=ON
    CMAKE_BUILD_TYPE:STRING=${CTEST_BUILD_CONFIGURATION}
    CMAKE_CXX_FLAGS:STRING=${SITE_CXX_FLAGS}
    CMAKE_C_FLAGS:STRING=${SITE_C_FLAGS}
    MAKECOMMAND:FILEPATH=${CTEST_BUILD_COMMAND}
    GITCOMMAND:FILEPATH=${GITCOMMAND}
    SVNCOMMAND:FILEPATH=${SVNCOMMAND}
    CVSCOMMAND:FILEPATH=${CVSCOMMAND}
    QT_QMAKE_EXECUTABLE:FILEPATH=${CTEST_QMAKE_COMMAND}
    TubeTK_USE_SUPERBUILD:BOOL=${SCRIPT_TubeTK_USE_SUPERBUILD}
    TubeTK_USE_KWSTYLE:BOOL=ON
    CTK_DIR:PATH=${SITE_BINARY_DIR}/CTK-Build
    ITK_DIR:PATH=${SITE_BINARY_DIR}/Insight-Build
    GenerateCLP_DIR:PATH=${SITE_BINARY_DIR}/GenerateCLP-Build
    OpenIGTLink_DIR:PATH=${SITE_BINARY_DIR}/OpenIGTLink-Build
    " )
  
endif( NOT EXISTS "${CTEST_BINARY_DIRECTORY}/CMakeCache.txt" )

