##############################################################################
#
# Initially developed for:   TubeTK
#
# Copyright 2013 Kitware Inc. 28 Corporate Drive,
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

if( JsonCpp_DIR )
  if( EXISTS ${JsonCpp_DIR}/CMakeCache.txt )
    file( STRINGS ${JsonCpp_DIR}/CMakeCache.txt _source_dir_def REGEX "jsoncpp_SOURCE_DIR" )
    string( REGEX REPLACE "[^=]+[=]" "" _source_dir "${_source_dir_def}" )
    set( _jsoncpp_include_dir "${_source_dir}/include" )
  endif()
  set( _jsoncpp_library
    ${JsonCpp_DIR}/src/lib_json
    ${JsonCpp_DIR}/src/lib_json/Release
    ${JsonCpp_DIR}/src/lib_json/MinSizeRel
    ${JsonCpp_DIR}/src/lib_json/RelWithDebInfo
    ${JsonCpp_DIR}/src/lib_json/Debug
    ${JsonCpp_DIR}/lib
    ${JsonCpp_DIR}/lib/Release
    ${JsonCpp_DIR}/lib/MinSizeRel
    ${JsonCpp_DIR}/lib/RelWithDebInfo
    ${JsonCpp_DIR}/lib/Debug )
  set( _jsoncpp_options NO_DEFAULT_PATH )
endif()

find_path( JsonCpp_INCLUDE_DIR NAMES json/json.h
  HINTS ${_jsoncpp_include_dir}
  PATH_SUFFIXES jsoncpp
  ${_jsoncpp_options} )

find_library( JsonCpp_LIBRARY NAMES jsoncpp libjsoncpp
  HINTS ${_jsoncpp_library}
  ${_jsoncpp_options} )

set( JsonCpp_INCLUDE_DIRS ${JsonCpp_INCLUDE_DIR} )
set( JsonCpp_LIBRARIES ${JsonCpp_LIBRARY} )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( JsonCpp DEFAULT_MSG JsonCpp_LIBRARIES
  JsonCpp_INCLUDE_DIRS )

# Do not use because of https://gitlab.kitware.com/cmake/cmake/issues/15448
# mark_as_advanced( JsonCpp_INCLUDE_DIR JsonCpp_LIBRARY )
