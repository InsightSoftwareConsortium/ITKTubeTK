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

project( TubeCLI )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_CLI_H_Files
  tubeCLIFilterWatcher.h
  tubeCLIHelperFunctions.h
  tubeCLIProgressReporter.h
  TubeTKLogo.h )

add_library(TubeCLI INTERFACE)
target_include_directories( TubeCLI INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR} )
target_link_libraries( TubeCLI INTERFACE TubeTKCommon )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_CLI_H_Files}
  DESTINATION include )
