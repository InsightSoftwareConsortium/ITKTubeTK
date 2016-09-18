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
cmake_minimum_required( VERSION 2.8.9 )
project( TubeTKITK )
set(TubeTKITK_LIBRARIES TubeTKITK)

if( NOT ITK_SOURCE_DIR )

  find_package( ITK REQUIRED )
  list( APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR} )
  include( ITKModuleExternal )

else()

  set( TubeTK_INSTALL_DEVELOPMENT 1 )
  set( TubeTK_INSTALL_INCLUDE_DIR ${ITK_INSTALL_INCLUDE_DIR} )
  set( TubeTK_BUILD_TESTING 0 )

  # Include depends file to set TubeTKITK_DEPENDS in current scope
  include( tubetk-depends.cmake )

  # Explicitly enable TubeTKBase libraries needed by TubeTKITK
  foreach( name ${${PROJECT_NAME}_DEPENDS} )
    set( TubeTK_BUILD_BASE_${name} 1 )
  endforeach()

  # Build TubeTKBase
  add_subdirectory(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Base
    ${CMAKE_CURRENT_BINARY_DIR}/TubeTKBase-build
    )

  itk_module_impl()

endif()
