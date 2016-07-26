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

if( BUILD_TESTING_VIRTUALENV AND TubeTK_BUILD_TESTING AND
  NOT TubeTK_BUILD_USING_SLICER )
  # Install into the testing virtualenv.
  execute_process( COMMAND "${PYTHON_EXECUTABLE}" "-c"
    "import virtualenv;print(virtualenv.__file__)"
    RESULT_VARIABLE virtualEnv_failed
    OUTPUT_VARIABLE virtualEnv_output
    ERROR_VARIABLE virtualEnv_error )
  if( virtualEnv_failed )
    execute_process(
      COMMAND ${PYTHON_TESTING_EXECUTABLE} -c
        "from distutils.sysconfig import get_python_lib; print(get_python_lib())"
      RESULT_VARIABLE pythonLib_result
      ERROR_VARIABLE pythonLib_error
      OUTPUT_VARIABLE site_packages )
    if( ${pythonLib_result} )
      message( SEND_ERROR "Could not get virtual packages directory:
        ${pythonLib_result}
        ${pythonLib_error}" )
    else()
      string( STRIP "${site_packages}" site_packages )
      message( STATUS
        "Detected testing virtualenv packages directory: ${site_packages}" )
      configure_file( tubetk.pth.in "${site_packages}/tubetk.pth" )
    endif()
  endif()
endif()

add_subdirectory( tubetk )

if( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )
  add_subdirectory( Testing )
endif( BUILD_TESTING AND NOT TubeTK_BUILD_USING_SLICER )
