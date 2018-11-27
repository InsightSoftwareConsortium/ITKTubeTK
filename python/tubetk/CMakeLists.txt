##############################################################################
#
# Library:   TubeTK
#
# Copyright Kitware Inc.
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
find_package( PythonInterp )
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE _numpy_include
  OUTPUT_STRIP_TRAILING_WHITESPACE )
find_path( NUMPY_INCLUDE_DIR numpy/arrayobject.h
  HINTS ${_numpy_include} )

find_package( PythonLibs REQUIRED )
include_directories( ${NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS} )

if( PYTHON_VERSION_STRING VERSION_GREATER 3.0 )
  add_definitions( "-DPYTHON3")
endif()

python_add_module( _tubetk_numpy _tubetk_numpy.cxx )
set_target_properties( _tubetk_numpy
  PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )

target_link_libraries( _tubetk_numpy
  TubeTK ${ITK_LIBRARIES} ${PYTHON_LIBRARIES} )

add_custom_target( CopyPurePythonTubeTK ALL "${CMAKE_COMMAND}" -P
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake/CopyTubeTK.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}"
  "*.py"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" )

add_custom_target( CopyPythonLibTubeTK ALL "${CMAKE_COMMAND}" -P
  "${CMAKE_CURRENT_SOURCE_DIR}/CMake/CopyTubeTK.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/../../lib"
  "*pyd*"
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}" )
add_dependencies( CopyPythonLibTubeTK _tubetk_numpy )
