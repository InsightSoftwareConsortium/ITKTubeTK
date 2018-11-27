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

project( TubeTKLibCommon )

if(NOT ITK_SOURCE_DIR)
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
endif()

set( TubeTKLib_Base_Common_H_Files
  tubeIndent.h
  tubeMacro.h
  tubeMessage.h
  tubeObject.h
  tubeStringUtilities.h
  tubeTestMain.h )

set( TubeTKLib_Base_Common_HXX_Files )

set( TubeTKLib_Base_Common_CXX_Files
  tubeIndent.cxx
  tubeObject.cxx )

add_library( ${PROJECT_NAME} STATIC
  ${TubeTKLib_Base_Common_H_Files}
  ${TubeTKLib_Base_Common_HXX_Files}
  ${TubeTKLib_Base_Common_CXX_Files} )
target_link_libraries( ${PROJECT_NAME}
  ${ITK_LIBRARIES} )
target_include_directories( ${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )

if( TubeTKLib_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTKLib_Base_Common_H_Files}
    ${TubeTKLib_Base_Common_HXX_Files}
    DESTINATION ${TubeTKLib_INSTALL_INCLUDE_DIR}
    COMPONENT Development )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS ${PROJECT_NAME} )
