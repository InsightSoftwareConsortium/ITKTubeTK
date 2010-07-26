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

include( ${CTEST_SCRIPT_DIRECTORY}/../../tubetk_config.cmake )

set( CTEST_CTEST_COMMAND ${SITE_CTEST_COMMAND} )

if( SITE_CONTINUOUS_BUILD_TEST )

  ctest_empty_binary_directory( "${SITE_BINARY_DIR}" )

  set( ENV{TUBETK_RUN_MODEL} "Continuous" )
  set( ENV{TUBETK_FORCE_BUILD} "1" )
  
  ###########################################################################
  # run some "inside-the-loop" continuous scripts for a while
  #
  while( ${CTEST_ELAPSED_TIME} LESS 56000 )
  
    set( START_TIME ${CTEST_ELAPSED_TIME} )
  
    ctest_run_script( "${SITE_SCRIPT_DIR}/tubetk_build_test.cmake" )

    if( "$ENV{TUBETK_FORCE_BUILD}" STREQUAL "1" )
  
      if( SITE_CONTINUOUS_STYLE )
        ctest_run_script( "${SITE_SCRIPT_DIR}/tubetk_style.cmake" )
      endif( SITE_CONTINUOUS_STYLE )
    
      if( SITE_CONTINUOUS_COVERAGE )
        ctest_run_script( "${SITE_SCRIPT_DIR}/tubetk_coverage.cmake" )
      endif( SITE_CONTINUOUS_COVERAGE )
    
      if( SITE_CONTINUOUS_MEMORY )
        ctest_run_script( "${SITE_SCRIPT_DIR}/tubetk_memory.cmake" )
      endif( SITE_CONTINUOUS_MEMORY )
  
      set( ENV{TUBETK_FORCE_BUILD} "0" )

    endif()

    # loop no faster than once every 2 minutes
    ctest_sleep( 120 )
  
  endwhile()

endif()

# do not run this script as a dashboard
set(CTEST_RUN_CURRENT_SCRIPT 0)

