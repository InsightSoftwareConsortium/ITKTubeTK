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
#   the tubetk source directory.
#
set( SITE_NAME "Abydos" )
set( SITE_PLATFORM "Linux" )
set( SITE_BUILD_TYPE "Release" )
set( SITE_CMAKE_GENERATOR "Unix Makefiles" )

set( ENV{DISPLAY} ":0" )

set( SITE_MAKE_COMMAND "make" )
set( SITE_CMAKE_COMMAND "/usr/local/bin/cmake" )
set( SITE_QMAKE_COMMAND "/usr/bin/qmake" )
set( SITE_CTEST_COMMAND "/usr/local/bin/ctest" )
set( SITE_MEMORYCHECK_COMMAND "/usr/bin/valgrind" )
set( SITE_COVERAGE_COMMAND "/usr/bin/gcov" )
set( SITE_STYLE_COMMAND "/usr/local/bin/kwstyle" )
set( GIT_COMMAND "/usr/bin/git" )
set( SITE_UPDATE_COMMAND "/usr/bin/git" )

set( SITE_SOURCE_DIR "/home/aylward/src/tubetk" )
set( SITE_BINARY_DIR "/home/aylward/src/tubetk-Release" )

set( SITE_NIGHTLY ON )
set( SITE_NIGHTLY_STYLECHECK ON )
set( SITE_NIGHTLY_COVERAGE ON )
set( SITE_NIGHTLY_MEMORYCHECK ON )
set( SITE_CONTINUOUS ON )
set( SITE_CONTINUOUS_STYLECHECK ON )
set( SITE_CONTINUOUS_COVERAGE ON )
set( SITE_CONTINUOUS_MEMORYCHECK ON )


###########################################################################
# 
# The following advanced variables should only be changed by experts
#
set( SITE_BUILD_NAME "${SITE_PLATFORM}-${SITE_BUILD_TYPE}" )
set( SITE_SCRIPT_DIR "${SITE_SOURCE_DIR}/DashboardScripts" )

set( SITE_CXX_FLAGS "-fPIC -fdiagnostics-show-option -W -Wall -Wextra -Wshadow -Wno-system-headers -Wwrite-strings -Wno-deprecated -Woverloaded-virtual" )
set( SITE_C_FLAGS "-fPIC -fdiagnostics-show-option -W -Wall -Wextra -Wshadow -Wno-system-headers -Wwrite-strings" )

set( COVERAGE_OPTIONS "-fprofile-arcs -ftest-coverage -lgcov" )
if( SITE_NIGHTLY_COVERAGE OR SITE_CONTINUOUS_COVERAGE )
  set( SITE_CXX_FLAGS "${SITE_CXX_FLAGS} ${COVERAGE_OPTIONS}" )
endif( SITE_NIGHTLY_COVERAGE OR SITE_CONTINUOUS_COVERAGE )

set( MEMORYCHECK_OPTIONS "--gen-suppressions=all --trace-children=yes -q --leak-check=yes --show-reachable=yes --num-callers=50" )

