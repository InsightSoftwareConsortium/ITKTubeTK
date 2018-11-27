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

set( gen "${CMAKE_GENERATOR}" )

ExternalProject_Add( MinimalPathExtraction 
  ${MinimalPathExtraction_EP_ARGS}
  GIT_REPOSITORY https://github.com/InsightSoftwareConsortium/ITKMinimalPathExtraction
  DOWNLOAD_DIR ${TubeTKLib_BINARY_DIR}/MinimalPathExtraction
  SOURCE_DIR ${TubeTKLib_BINARY_DIR}/MinimalPathExtraction
  BINARY_DIR ${TubeTKLib_BINARY_DIR}/MinimalPathExtraction-build
  INSTALL_DIR ${TubeTKLib_BINARY_DIR}/MinimalPathExtraction-build
  INSTALL_COMMAND ""
  CMAKE_GENERATOR ${gen}
  BUILD_ALWAYS 1
  CMAKE_CACHE_ARGS
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS_INIT}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS_INIT}
    -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
    -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
    -DCMAKE_BUILD_SHARED_LIBS:BOOL=OFF
    -DCMAKE_INSTALL_PREFIX:PATH=${MinimalPathExtraction_DIR}
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
    -DITK_INSTALL_RUNTIME_DIR:PATH=${ITK_INSTALL_RUNTIME_DIR}
    -DITK_INSTALL_LIBRARY_DIR:PATH=${ITK_INSTALL_LIBRARY_DIR}
  DEPENDS
     ""
  )

set( MinimalPathExtraction_DIR 
  ${TubeTKLib_BINARY_DIR}/MinimalPathExtraction-build )
