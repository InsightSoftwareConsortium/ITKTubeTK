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
#
# Configure tubetk.pth.
# Add a symlink to this file in the directory that results from::
#
#   python -m site --user-site
#
# to make TubeTK Python packages and modules importable.
#
set( TubeTK_C_PYTHON_EXTENSION_DIR ${CMAKE_CURRENT_BINARY_DIR} )
configure_file( tubetk.pth.in tubetk.pth )
if( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )
  # Install into the testing virtualenv.
  execute_process( COMMAND ${PYTHON_TESTING_EXECUTABLE}
      -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
    RESULT_VARIABLE _result_variable
    ERROR_VARIABLE _error_variable
    OUTPUT_VARIABLE site_packages )
  if( ${_result_variable} )
    message( SEND_ERROR "Could not get virtual packages directory:
      ${_result_variable}
      ${_error_variable}" )
  else( ${_result_variable} )
    string( STRIP "${site_packages}" site_packages )
    message( STATUS "Detected testing virtualenv packages directory: ${site_packages}" )
    configure_file( tubetk.pth.in "${site_packages}/tubetk.pth" )
  endif( ${_result_variable} )
endif( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )

include_directories( ${TubeTK_SOURCE_DIR}/Base/Filtering )

add_subdirectory( tubetk )

if( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )
  add_subdirectory( Testing )
endif( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )
