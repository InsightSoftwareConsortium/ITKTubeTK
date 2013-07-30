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

if( JsonCpp_DIR )
  set( JsonCpp_INCLUDEDIR ${JsonCpp_DIR}/include )
  set( JsonCpp_LIBDIR ${JsonCpp_DIR}/lib )
endif( JsonCpp_DIR )

find_package( PkgConfig QUIET )
pkg_check_modules( PC_JsonCpp QUIET json jsoncpp )

find_path( JsonCpp_INCLUDE_DIR NAMES json/json.h
  HINTS ${JsonCpp_INCLUDEDIR} ${PC_JsonCpp_INCLUDEDIR} ${PC_JsonCpp_INCLUDE_DIRS}
  PATH_SUFFIXES jsoncpp )

find_library( JsonCpp_LIBRARIES NAMES json jsoncpp libjson libjsoncpp
  HINTS ${JsonCpp_LIBDIR} ${PC_JsonCpp_LIBDIR} ${PC_JsonCpp_LIBRARY_DIRS} )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( JsonCpp DEFAULT_MSG JsonCpp_LIBRARIES
  JsonCpp_INCLUDE_DIR )

mark_as_advanced( JsonCpp_INCLUDE_DIR JsonCpp_LIBRARIES )
