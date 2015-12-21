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

if( ArrayFire_DIR )
  set( _arrayfire_include_dir
    "${ArrayFire_DIR}/include" )
  set( _arrayfire_library_dir
    "${ArrayFire_DIR}/lib"
    "${ArrayFire_DIR}/lib/Release"
    "${ArrayFire_DIR}/lib/MinSizeRel"
    "${ArrayFire_DIR}/lib/RelWithDebInfo"
    "${ArrayFire_DIR}/lib/Debug" )
endif( ArrayFire_DIR )

find_path( ArrayFire_INCLUDE_DIR NAMES arrayfire.h
  HINTS
    ${_arrayfire_include_dir}
    "C:/Program Files/ArrayFire/v3"
    "C:/Program Files/ArrayFire/v3/include"
    /usr/local/include
    /usr/include
    /opt/local/include )

find_library( ArrayFire_LIBRARY
  NAMES af afcpu libaf libafcpu libafcuda libafopencl
  HINTS
    ${_arrayfire_library_dir}
    "C:/Program Files/ArrayFire/v3"
    "C:/Program Files/ArrayFire/v3/lib"
    /usr/local/lib )

set( ArrayFire_INCLUDE_DIRS ${ArrayFire_INCLUDE_DIR} )
set( ArrayFire_LIBRARIES ${ArrayFire_LIBRARY} )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( ArrayFire DEFAULT_MSG ArrayFire_LIBRARIES
  ArrayFire_INCLUDE_DIRS )

mark_as_advanced( ArrayFire_INCLUDE_DIR ArrayFire_LIBRARY )
