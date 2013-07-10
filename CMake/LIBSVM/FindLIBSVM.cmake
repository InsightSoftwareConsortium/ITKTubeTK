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

# Sets LIBSVM_FOUND, LIBSVM_INCLUDE_DIR, and LIBSVM_LIBRARY.

if( LIBSVM_DIR )
  set( _svm_lib_dir "${LIBSVM_DIR}/lib" )
  set( _svm_include_dir "${LIBSVM_DIR}/include" )
endif( LIBSVM_DIR )

find_path( LIBSVM_INCLUDE_DIR svm.h PATHS ${_svm_include_dir} DOC "Path to the LIBSVM include directory" )
find_library( LIBSVM_LIBRARY svm PATHS ${_svm_lib_dir} DOC "The LIBSVM library." )

mark_as_advanced( LIBSVM_INCLUDE_DIR )
mark_as_advanced( LIBSVM_LIBRARY )

set( LIBSVM_FOUND 0 )

if ( LIBSVM_INCLUDE_DIR AND LIBSVM_LIBRARY )
  set( LIBSVM_FOUND 1 )
  set( LIBSVM_LIBRARIES ${LIBSVM_LIBRARY} )
endif( LIBSVM_INCLUDE_DIR AND LIBSVM_LIBRARY )
