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

if( LIBSVM_DIR )
  set( _LIBSVM_include_dir "${LIBSVM_DIR}/include" )
  set( _LIBSVM_library "${LIBSVM_DIR}/lib"
    "${LIBSVM_DIR}/lib/Release"
    "${LIBSVM_DIR}/lib/MinSizeRel"
    "${LIBSVM_DIR}/lib/RelWithDebInfo"
    "${LIBSVM_DIR}/lib/Debug" )
endif( LIBSVM_DIR )

find_path( LIBSVM_INCLUDE_DIR NAMES svm.h
  HINTS ${_LIBSVM_include_dir} )

find_library( LIBSVM_LIBRARY NAMES svm libsvm
  HINTS ${_LIBSVM_library} )

set( LIBSVM_INCLUDE_DIRS ${LIBSVM_INCLUDE_DIR} )
set( LIBSVM_LIBRARIES ${LIBSVM_LIBRARY} )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( LIBSVM DEFAULT_MSG LIBSVM_LIBRARIES
  LIBSVM_INCLUDE_DIRS )

mark_as_advanced( LIBSVM_INCLUDE_DIR LIBSVM_LIBRARY )
