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

cmake_minimum_required(VERSION 2.8)

include( ${CTEST_SCRIPT_DIRECTORY}/../../tubetk_config.cmake )

set( SCRIPT_NAME "Style" )
set( SCRIPT_BINARY_SUBDIR "TubeTK-Build" )
set( SCRIPT_TubeTK_USE_SUPERBUILD OFF )

include( ${CTEST_SCRIPT_DIRECTORY}/tubetk_cmakecache.cmake )

set( ENV{PATH} "$ENV{PATH}:/usr/local/bin" )

ctest_start( "$ENV{TUBETK_RUN_MODEL}" )
ctest_configure( BUILD "${CTEST_BINARY_DIRECTORY}" )
ctest_read_custom_files( "${CTEST_BINARY_DIRECTORY}" )
EXECUTE_PROCESS( COMMAND make -C "${CTEST_BINARY_DIRECTORY}" StyleCheck )
ctest_submit()

