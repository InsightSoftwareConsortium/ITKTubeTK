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

# Build the documentation
include( Documentation OPTIONAL )

if( BUILD_DOCUMENTATION AND TubeTK_USE_DOXYGEN )
  option( DOXYGEN_SHORT_NAMES "Build Doxygen of CURRENT using short file names" OFF )
  mark_as_advanced( DOXYGEN_SHORT_NAMES )
  if( DOXYGEN_SHORT_NAMES )
    set( DOXYGEN_SHORT_NAMES YES )
  else( DOXYGEN_SHORT_NAMES )
    set( DOXYGEN_SHORT_NAMES NO )
  endif( DOXYGEN_SHORT_NAMES )

  # Configure the script and the doxyfile, then add target
  if( NOT DOT_PATH )
    get_filename_component( DOT_PATH ${DOT} PATH )
  endif( NOT DOT_PATH )

  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/doxyfile )

endif( BUILD_DOCUMENTATION AND TubeTK_USE_DOXYGEN )
