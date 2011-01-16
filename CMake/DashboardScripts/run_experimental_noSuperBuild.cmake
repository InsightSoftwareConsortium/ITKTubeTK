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

set( CTEST_CTEST_COMMAND ${SITE_CTEST_COMMAND} )

set( ENV{TUBETK_RUN_MODEL} "Experimental" )
set( ENV{TUBETK_FORCE_BUILD} "0" )

if( SITE_EXPERIMENTAL_BUILD_TEST )
  include( "${SITE_SCRIPT_DIR}/build_test_noSuperBuild.cmake" )
ENDif( SITE_EXPERIMENTAL_BUILD_TEST )

if( SITE_EXPERIMENTAL_STYLE )
  include( "${SITE_SCRIPT_DIR}/style.cmake" )
endif( SITE_EXPERIMENTAL_STYLE )

if( SITE_EXPERIMENTAL_COVERAGE )
  include( "${SITE_SCRIPT_DIR}/coverage.cmake" )
endif( SITE_EXPERIMENTAL_COVERAGE )

if( SITE_EXPERIMENTAL_MEMORY )
  include( "${SITE_SCRIPT_DIR}/memory.cmake" )
endif( SITE_EXPERIMENTAL_MEMORY )

set(CTEST_RUN_CURRENT_SCRIPT 0)
