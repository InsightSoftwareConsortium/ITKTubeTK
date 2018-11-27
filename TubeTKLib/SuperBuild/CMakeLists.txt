##############################################################################
#
# Library:   TubeTKLib
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

include( ExternalProject )

set( superbuild_base_dir "${TubeTKLib_BINARY_DIR}" )
set_property( DIRECTORY PROPERTY EP_BASE ${superbuild_base_dir} )

set( gen "${CMAKE_GENERATOR}" )

include( MinimalPathExtraction.cmake )

ExternalProject_Add( TubeTKLib
  ${TubeTKLib_EP_ARGS}
  DOWNLOAD_COMMAND ""
  DOWNLOAD_DIR ${TubeTKLib_SOURCE_DIR}
  SOURCE_DIR ${TubeTKLib_SOURCE_DIR}
  BINARY_DIR ${TubeTKLib_BINARY_DIR}
  INSTALL_COMMAND ""
  INSTALL_DIR ${TubeTKLib_BINARY_DIR}
  CMAKE_GENERATOR ${gen}
  BUILD_ALWAYS 1
  CMAKE_CACHE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_SHARED_LIBS:STRING=${CMAKE_BUILD_SHARED_LIBS}
    -DCMAKE_INSTALL_PREFIX:PATH=${TubeTKLib_BINARY_DIR}
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
    -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
    -DBUILD_TESTING:BOOL=${BUILD_TESTING}
    -DITK_DIR:PATH=${ITK_DIR}
    -DITK_INSTALL_RUNTIME_DIR:PATH=${ITK_INSTALL_RUNTIME_DIR}
    -DITK_INSTALL_LIBRARY_DIR:PATH=${ITK_INSTALL_LIBRARY_DIR}
    -DMinimalPathExtraction_DIR:PATH=${MinimalPathExtraction_DIR}
    -DMinimalPathExtraction_SOURCE_DIR:PATH=${MinimalPathExtraction_SOURCE_DIR}
    -DMinimalPathExtraction_BINARY_DIR:PATH=${MinimalPathExtraction_BINARY_DIR}
    -DTubeTKLib_SUPERBUILD:BOOL=OFF
    -DTubeTKLib_USE_VTK:BOOL=${TubeTKLib_USE_VTK}
    -DVTK_DIR:PATH=${VTK_DIR}
    -DTubeTKLib_USE_KWStyle:BOOL=${TubeTKLib_USE_KWStyle}
    -DTubeTKLib_USE_CPPCHECK:BOOL=${TubeTKLib_USE_CPPCHECK}
    -DTubeTKLib_USE_ARRAYFIRE:BOOL=${TubeTKLib_USE_ARRAYFIRE}
    -DTubeTKLib_USE_LIBSVM:BOOL=${TubeTKLib_USE_LIBSVM}
    -DTubeTKLib_USE_RANDOMFOREST:BOOL=${TubeTKLib_USE_RANDOMFOREST}
    -DTubeTKLib_USE_VALGRIND:BOOL=${TubeTKLib_USE_VALGRIND}
  DEPENDS
    "MinimalPathExtraction"
  )
