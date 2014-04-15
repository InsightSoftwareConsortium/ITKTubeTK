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
find_package( KWStyle REQUIRED )
find_package( Git )
if( GIT_FOUND )
  execute_process( COMMAND ${GIT_EXECUTABLE} config hooks.KWStyle.path
    "${KWSTYLE_EXECUTABLE}"
    WORKING_DIRECTORY ${TubeTK_SOURCE_DIR} )
endif( GIT_FOUND )

option( KWSTYLE_USE_VIM_FORMAT
  "Set KWStyle to generate errors in a Vim compatible format." OFF )
mark_as_advanced( KWSTYLE_USE_VIM_FORMAT )

option( KWSTYLE_USE_MSVC_FORMAT
  "Set KWStyle to generate errors in Microsoft Visual Studio format."
  OFF )
mark_as_advanced( KWSTYLE_USE_MSVC_FORMAT )

option( KWSTYLE_DASHBOARD_SUBMISSION
  "Set KWStyle to generate a report for CDash dashboard submission."
  OFF )
mark_as_advanced( KWSTYLE_DASHBOARD_SUBMISSION )

if( KWSTYLE_DASHBOARD_SUBMISSION )
  set( KWSTYLE_ARGUMENTS -lesshtml
    -xml ${TubeTK_BINARY_DIR}/KWStyle.kws.xml
    -o ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt
    -dart ${TubeTK_BINARY_DIR} -1 1
    -D ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )
else( KWSTYLE_DASHBOARD_SUBMISSION )
  set( KWSTYLE_ARGUMENTS -xml ${TubeTK_BINARY_DIR}/KWStyle.kws.xml
    -v -o ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt
    -D ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )
endif( KWSTYLE_DASHBOARD_SUBMISSION )

if( KWSTYLE_USE_VIM_FORMAT )
  set( KWSTYLE_ARGUMENTS -vim ${KWSTYLE_ARGUMENTS} )
endif( KWSTYLE_USE_VIM_FORMAT )

if( KWSTYLE_USE_MSVC_FORMAT )
  set( KWSTYLE_ARGUMENTS -msvc ${KWSTYLE_ARGUMENTS} )
endif( KWSTYLE_USE_MSVC_FORMAT )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.kws.xml.in
  ${TubeTK_BINARY_DIR}/KWStyle.kws.xml )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.Files.txt.in
  ${TubeTK_BINARY_DIR}/KWStyle.Files.txt )

configure_file( ${TubeTK_SOURCE_DIR}/CMake/KWStyle/KWStyle.Overwrite.txt.in
  ${TubeTK_BINARY_DIR}/KWStyle.Overwrite.txt )

add_custom_target( StyleCheck
  COMMAND ${KWSTYLE_EXECUTABLE} ${KWSTYLE_ARGUMENTS}
  WORKING_DIRECTORY ${TubeTK_BINARY_DIR}
  COMMENT "Coding style checker" VERBATIM )
