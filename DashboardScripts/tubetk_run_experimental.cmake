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

set( RUN_DASHBOARD_MODEL "Experimental" )

ctest_empty_binary_directory( "${SITE_BINARY_DIR}" )

ctest_run_script( 
  "${SITE_SCRIPT_DIR}/tubetk_build_test.cmake" )

if( SITE_CONTINUOUS_STYLECHECK )
  ctest_run_script( 
    "${SITE_SCRIPT_DIR}/tubetk_style.cmake" )
endif( SITE_CONTINUOUS_STYLECHECK )
  
if( SITE_CONTINUOUS_COVERAGE )
  ctest_run_script( 
    "${SCRIPT_DIR}/tubetk_coverage.cmake" )
endif( SITE_CONTINUOUS_COVERAGE )
 
if( SITE_CONTINUOUS_MEMORYCHECK )
  ctest_run_script( 
    "${SCRIPT_DIR}/tubetk_memory.cmake" )
endif( SITE_CONTINUOUS_MEMORYCHECK )

