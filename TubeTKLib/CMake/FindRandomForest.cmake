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

if( RandomForest_DIR )
  set( _RandomForest_include_dir "${RandomForest_DIR}/include" )
endif( RandomForest_DIR )

find_path( RandomForest_INCLUDE_DIR NAMES andres/marray.hxx
  HINTS ${_RandomForest_include_dir} )

set( RandomForest_INCLUDE_DIR ${RandomForest_INCLUDE_DIR} )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( RandomForest DEFAULT_MSG
  RandomForest_INCLUDE_DIR )

mark_as_advanced( RandomForest_INCLUDE_DIR )
