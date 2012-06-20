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

##############################################################################
#
# Configure the following variables and move this file to the directory above
#   the TubeTK source directory.
#
set( SITE_NAME "Ginger.Aylward.Org" )
set( SITE_PLATFORM "WindowsXP-VS2010_32" )
set( SITE_BUILD_TYPE "Release" )
set( SITE_CTEST_MODE "Experimental" ) # Experimental, Continuous, or Nightly
set( SITE_CMAKE_GENERATOR "Visual Studio 10" )

set( TUBETK_GIT_REPOSITORY "git://tubetk.org/TubeTK.git" )
set( TUBETK_SOURCE_DIR
  "C:\\Users\\aylward\\src\\TubeTK" )
set( TUBETK_BINARY_DIR
  "C:\\Users\\aylward\\src\\TubeTK-${SITE_BUILD_TYPE}" )

set( SITE_MAKE_COMMAND "${CTEST_BUILD_COMMAND}" )
set( SITE_CMAKE_COMMAND "C:\\Program Files\\CMake 2.8\\bin\\cmake" )
set( SITE_QMAKE_COMMAND "C:\\Qt\\4.6.3\\bin\\qmake" )
set( SITE_CTEST_COMMAND "C:\\Program Files\\CMake 2.8\\bin\\ctest" )

set( SITE_MEMORYCHECK_COMMAND "" )
set( SITE_COVERAGE_COMMAND "" )
set( SITE_KWSTYLE_DIR "" )

set( SITE_GIT_COMMAND "C:\\Program Files\\Git\\bin\\git" )
set( SITE_SVN_COMMAND "C:\\Program Files\\CollabNet\\Subversion Client\\svn" )

set( SITE_EXPERIMENTAL_BUILD ON )
set( SITE_EXPERIMENTAL_TEST ON )
set( SITE_EXPERIMENTAL_COVERAGE OFF ) # Linux only
set( SITE_EXPERIMENTAL_MEMORY OFF )   # Linux only
set( SITE_EXPERIMENTAL_PACKAGE OFF )
set( SITE_EXPERIMENTAL_UPLOAD OFF )
set( SITE_EXPERIMENTAL_STYLE OFF )    # Linux only

set( SITE_CONTINUOUS_BUILD ON )
set( SITE_CONTINUOUS_TEST ON )
set( SITE_CONTINUOUS_COVERAGE OFF ) # Linux only
set( SITE_CONTINUOUS_MEMORY OFF )   # Linux only
set( SITE_CONTINUOUS_PACKAGE ON )
set( SITE_CONTINUOUS_UPLOAD ON )
set( SITE_CONTINUOUS_STYLE OFF )    # Linux only

set( SITE_NIGHTLY_BUILD ON )
set( SITE_NIGHTLY_TEST ON )
set( SITE_NIGHTLY_COVERAGE OFF ) # Linux only
set( SITE_NIGHTLY_MEMORY OFF )   # Linux only
set( SITE_NIGHTLY_PACKAGE ON )
set( SITE_NIGHTLY_UPLOAD ON )
set( SITE_NIGHTLY_STYLE OFF )    # Linux only


###########################################################################
#
# The following advanced variables should only be changed by experts
#
set( CTEST_SOURCE_DIRECTORY "${TUBETK_SOURCE_DIR}" )
set( CTEST_BINARY_DIRECTORY "${TUBETK_BINARY_DIR}" )
set( CTEST_COMMAND "${SITE_CTEST_COMMAND}" )

set( SITE_BUILD_NAME "${SITE_PLATFORM}-${SITE_BUILD_TYPE}" )
set( TUBETK_SCRIPT_DIR "${TUBETK_SOURCE_DIR}/CMake/DashboardScripts" )
set( SITE_UPDATE_COMMAND "${SITE_GIT_COMMAND}" )
set( GIT_EXECUTABLE "${SITE_GIT_COMMAND}" )

set( SITE_EXECUTABLE_DIRS "${TUBETK_BINARY_DIR}/ModuleDescriptionParser-Build/${SITE_BUILD_TYPE};${TUBETK_BINARY_DIR}/GenerateCLP-Build/${SITE_BUILD_TYPE};${TUBETK_BINARY_DIR}/Insight-Build/bin/${SITE_BUILD_TYPE};${TUBETK_BINARY_DIR}/VTK-Build/bin/${SITE_BUILD_TYPE};${TUBETK_BINARY_DIR}/TubeTK-Build/bin/${SITE_BUILD_TYPE};${TUBETK_BINARY_DIR}/TubeTK-Build/lib/TubeTK/Plugins/${SITE_BUILD_TYPE}" )
set( ENV{PATH} "${SITE_EXECUTABLE_DIRS};$ENV{PATH}" )

set( SITE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /GR /MP /EHsc" )
set( SITE_C_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zm1000 /MP" )
set( SITE_EXE_LINKER_FLAGS "" )
set( SITE_SHARED_LINKER_FLAGS "" )

set( SITE_MEMORYCHECK_COMMAND_OPTIONS "" )
set( SITE_MEMORYCHECK_SUPPRESSIONS_FILE
  "${CTEST_SCRIPT_DIRECTORY}/valgrind_suppressions.txt" )

if( NOT EXISTS "${TUBETK_SOURCE_DIR}/CMakeLists.txt" )
  execute_process( COMMAND
    "${SITE_GIT_COMMAND}"
    clone "${TUBETK_GIT_REPOSITORY}" "${TUBETK_SOURCE_DIR}" )
  ctest_run_script()
endif()

if( "${SITE_CTEST_MODE}" STREQUAL "Experimental" )
  include( "${TUBETK_SCRIPT_DIR}/run_experimental.cmake" )
elseif( "${SITE_CTEST_MODE}" STREQUAL "Continuous" )
  include( "${TUBETK_SCRIPT_DIR}/run_continuous.cmake" )
elseif( "${SITE_CTEST_MODE}" STREQUAL "Nightly" )
  include( "${TUBETK_SCRIPT_DIR}/run_nightly.cmake" )
endif()
