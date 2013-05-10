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

# If you follow the SampleCLIApplication format, you only need to change the
# following line to configure this CMakeLists.txt file.
project( ImageMath )

cmake_minimum_required( VERSION 2.8 )
if( COMMAND CMAKE_POLICY )
  cmake_policy( SET CMP0003 NEW )
endif( COMMAND CMAKE_POLICY )

# Disable MSVC 8 warnings
if( WIN32 )
  option( DISABLE_MSVC8_DEPRECATED_WARNINGS
    "Disable Visual Studio 8 deprecated warnings" ON )
  mark_as_advanced( FORCE DISABLE_MSVC8_DEPRECATED_WARNINGS )
  if( DISABLE_MSVC8_DEPRECATED_WARNINGS )
    add_definitions( -D_CRT_SECURE_NO_DEPRECATE )
  endif( DISABLE_MSVC8_DEPRECATED_WARNINGS )
endif( WIN32)

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

# Include Utilities to access tubeCLIHelperFunctions.h
include_directories(
  ${TubeTK_SOURCE_DIR}/Utilities
  ${TubeTK_SOURCE_DIR}/Base/Numerics
  ${TubeTK_SOURCE_DIR}/Base/Filtering
  ${TubeTK_SOURCE_DIR}/Base/Segmentation
  )

set( PROJECT_SOURCE
  ${PROJECT_NAME}.cxx
  )

add_executable( ${PROJECT_NAME} ${PROJECT_SOURCE} )
target_link_libraries( ${PROJECT_NAME} ${ITK_LIBRARIES} TubeNumerics )

# Add the install rule
install( FILES DESTINATION bin COMPONENT RUNTIME )

# Installation file
install(
  TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  COMPONENT RUNTIME )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
