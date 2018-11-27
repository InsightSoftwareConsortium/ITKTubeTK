##############################################################################
#
# Library:   TubeTKLib
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

find_program( KWSTYLE_EXECUTABLE NAMES KWStyle
  PATHS /usr/local/bin
  /usr/bin
  ${KWStyle_DIR}/bin
  ${KWStyle_DIR}/bin/Release
  ${KWStyle_DIR}/bin/MinSizeRel
  ${KWStyle_DIR}/bin/RelWithDebInfo
  ${KWStyle_DIR}/bin/Debug )
mark_as_advanced( KWSTYLE_EXECUTABLE )

find_package( Git )
if( GIT_FOUND )
  execute_process( COMMAND ${GIT_EXECUTABLE} config hooks.KWStyle.path
    "${KWSTYLE_EXECUTABLE}"
    WORKING_DIRECTORY ${TubeTKLib_SOURCE_DIR} )
endif( GIT_FOUND )

option( KWSTYLE_USE_VIM_FORMAT
  "Set KWStyle to generate errors in a Vim compatible format." OFF )
mark_as_advanced( KWSTYLE_USE_VIM_FORMAT )

option( KWSTYLE_USE_MSVC_FORMAT
  "Set KWStyle to generate errors in Microsoft Visual Studio format."
  OFF )
mark_as_advanced( KWSTYLE_USE_MSVC_FORMAT )

set( KWSTYLE_DASHBOARD_ARGUMENTS
  -dart ${TubeTKLib_BINARY_DIR} -1 1
  -xml ${TubeTKLib_BINARY_DIR}/KWStyle.kws.xml
  -o ${TubeTKLib_BINARY_DIR}/KWStyle.Overwrite.txt
  -R -D ${TubeTKLib_BINARY_DIR}/KWStyle.Files.txt )
set( KWSTYLE_ARGUMENTS -v -xml ${TubeTKLib_BINARY_DIR}/KWStyle.kws.xml
  -o ${TubeTKLib_BINARY_DIR}/KWStyle.Overwrite.txt
  -R -D ${TubeTKLib_BINARY_DIR}/KWStyle.Files.txt )
if( KWSTYLE_USE_VIM_FORMAT )
  set( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
else( KWSTYLE_USE_VIM_FORMAT )
  if( KWSTYLE_USE_MSVC_FORMAT )
    set( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
  endif( KWSTYLE_USE_MSVC_FORMAT )
endif( KWSTYLE_USE_VIM_FORMAT )

configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/KWStyle/KWStyle.kws.xml.in
  ${TubeTKLib_BINARY_DIR}/KWStyle.kws.xml )

configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/KWStyle/KWStyle.Files.txt.in
  ${TubeTKLib_BINARY_DIR}/KWStyle.Files.txt )

configure_file( ${TubeTKLib_SOURCE_DIR}/CMake/KWStyle/KWStyle.Overwrite.txt
  ${TubeTKLib_BINARY_DIR}/KWStyle.Overwrite.txt )

add_custom_command( OUTPUT KWStyle_Dashboard_Log.txt
  COMMAND ${KWSTYLE_EXECUTABLE}
  ARGS ${KWSTYLE_DASHBOARD_ARGUMENTS} )
add_custom_target( StyleCheckDashboard
  DEPENDS KWStyle_Dashboard_Log.txt )

add_custom_command( OUTPUT KWStyle_Log.txt
  COMMAND ${KWSTYLE_EXECUTABLE}
  ARGS ${KWSTYLE_ARGUMENTS} )
add_custom_target( StyleCheck
  DEPENDS KWStyle_Log.txt )
