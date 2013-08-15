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

include( Midas3FunctionAddTest )
set( MIDAS_REST_URL http://midas3.kitware.com/midas/api/rest )
set( MIDAS_KEY_DIR ${TubeTK_SOURCE_DIR}/MIDAS_Keys )

find_package( PythonInterp )

# This sets the virtual environment directory
set( PythonVirtualenvHome ${TubeTK_BINARY_DIR}/Temporary/PythonVirtualenv )

# Make the original Python executable explicit
set( ORIGINAL_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE} )

# Temporarily change the Python executable to the virtual environment one
if( WIN32 )
  set( PYTHON_TESTING_EXECUTABLE ${PythonVirtualenvHome}/Scripts/python )
else( )
  set( PYTHON_TESTING_EXECUTABLE ${PythonVirtualenvHome}/bin/python )
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/PythonVirtualEnvInstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/PythonVirtualEnvInstall.cmake"
  @ONLY )
  
# Install the virtual environment
execute_process( COMMAND ${CMAKE_COMMAND} -P 
  "${CMAKE_CURRENT_BINARY_DIR}/PythonVirtualEnvInstall.cmake" )

# Setup IPython notebook driver  
set( NOTEBOOK_TEST_DRIVER 
  ${TubeTK_SOURCE_DIR}/Utilities/Python/EvaluateIPythonNotebook.py)

# We need SimpleITK (almost surely) in all examples
execute_process( COMMAND ${PythonVirtualenvHome}/bin/easy_install SimpleITK )

# These are the most commonly required Python packages
set( PYTHON_TESTING_MODULES 
  tornado 
  ipython 
  ipython[zmq]
  jinja2
  numpy
  matplotlib )

# Install Python packages via pip install in the virtual environment 
foreach ( MODULE ${PYTHON_TESTING_MODULES} )
  execute_process( COMMAND ${PythonVirtualenvHome}/bin/pip install ${MODULE} )
endforeach()
  
# Add notebook/example directories
add_subdirectory( MergeAdjacentImages )




