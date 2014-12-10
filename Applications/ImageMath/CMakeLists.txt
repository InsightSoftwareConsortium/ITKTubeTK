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
endif( WIN32 )

if( NOT TubeTK_SOURCE_DIR )
  find_package( TubeTK REQUIRED )
  include( ${TubeTK_USE_FILE} )
endif( NOT TubeTK_SOURCE_DIR )

# Find ITK
find_package( ITK REQUIRED )
set( ITK_NO_IO_FACTORY_REGISTER_MANAGER 1 )
include( ${ITK_USE_FILE} )

# Find SEM
find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

SEMMacroBuildCLI(
  NAME ${PROJECT_NAME}
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} TubeTKNumerics ITKIOMeta ITKStatistics
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/CLI
    ${TubeTK_SOURCE_DIR}/Base/Common
    ${TubeTK_SOURCE_DIR}/Base/Numerics
    ${TubeTK_SOURCE_DIR}/Base/Filtering
    ${TubeTK_SOURCE_DIR}/Base/Segmentation
    ${TubeTK_SOURCE_DIR}/Base/Registration
  INSTALL_RUNTIME_DESTINATION lib
  INSTALL_LIBRARY_DESTINATION lib
  INSTALL_ARCHIVE_DESTINATION bin )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
