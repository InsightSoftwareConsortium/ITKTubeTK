##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
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

if( WIN32 )
  set( TubeTK_EXECUTABLE_DIRS CACHE INTERNAL
    "Bin and Lib dirs for running apps." FORCE )
  foreach( _build_type "" Debug Release )
    list( APPEND TubeTK_EXECUTABLE_DIRS
      ${ITK_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/bin/${_build_type}
      ${TubeTK_BINARY_DIR}/lib/${_build_type} )
    if( TubeTK_USE_VTK )
      list( APPEND TubeTK_EXECUTABLE_DIRS ${VTK_DIR}/bin/${_build_type} )
    endif()
  endforeach()
else( WIN32 )
  set( TubeTK_EXECUTABLE_DIRS
    ${ITK_DIR}/bin
    ${TubeTK_BINARY_DIR}/bin
    ${TubeTK_BINARY_DIR}/lib
    CACHE INTERNAL "Bin and Lib dirs for running apps." FORCE )
  if( TubeTK_USE_VTK )
    list( APPEND TubeTK_EXECUTABLE_DIRS ${VTK_DIR}/bin )
  endif()
endif( WIN32 )

message( STATUS "Configuring Launcher script" )
set( TubeTK_LAUNCHER )
if( WIN32 )
  set( _launcher_platform "windows" )
  foreach( _build_type "" Debug Release )
    configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.bat.in
      ${TubeTK_BINARY_DIR}/bin/${_build_type}/TubeTKLauncher.bat @ONLY )
  endforeach()

  find_program( CMD_EXECUTABLE "cmd" )
  if( NOT CMD_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'cmd' executable required to run tests and launcher" )
  endif()
  mark_as_advanced( CMD_EXECUTABLE )

  set( TubeTK_LAUNCHER ${TubeTK_BINARY_DIR}/Release/TubeTKLauncher.bat )
elseif( UNIX )
  set( _launcher_platform "unix" )
  configure_file( ${TubeTK_SOURCE_DIR}/CMake/TubeTKLauncher.sh.in
    ${TubeTK_BINARY_DIR}/bin/TubeTKLauncher.sh @ONLY )

  find_program( SH_EXECUTABLE "sh" )
  if( NOT SH_EXECUTABLE )
    message( FATAL_ERROR
      "Could not find 'sh' executable required to test using the launcher" )
  endif()
  set( TubeTK_LAUNCHER ${SH_EXECUTABLE}
    ${TubeTK_BINARY_DIR}/bin/TubeTKLauncher.sh )
elseif( NOT UNIX )
  message( FATAL_ERROR
    "Configuring Launcher script - failed [unknown platform]" )
endif()
  message( STATUS "Configuring Launcher script - ok [${_launcher_platform}]" )
