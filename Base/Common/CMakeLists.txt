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

project( TubeTKCommon )

if(NOT ITK_SOURCE_DIR)
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
endif()

set( TubeTK_Base_Common_H_Files
  itktubeDifferenceImageFilter.h
  itktubeImageMath.h
  tubeIndent.h
  tubeMacro.h
  tubeMessage.h
  tubeObject.h
  tubeStringUtilities.h
  tubeTestMain.h )

set( TubeTK_Base_Common_HXX_Files
  itktubeDifferenceImageFilter.hxx
  itktubeImageMath.hxx )

set( TubeTK_Base_Common_CXX_Files
  tubeIndent.cxx
  tubeObject.cxx )

add_library( ${PROJECT_NAME} STATIC
  ${TubeTK_Base_Common_H_Files}
  ${TubeTK_Base_Common_HXX_Files}
  ${TubeTK_Base_Common_CXX_Files} )
target_link_libraries( ${PROJECT_NAME}
  ${ITK_LIBRARIES} )
target_include_directories( ${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )

if( TubeTK_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTK_BUILD_TESTING )

if( TubeTK_INSTALL_DEVELOPMENT )
  install( FILES
    ${TubeTK_Base_Common_H_Files}
    ${TubeTK_Base_Common_HXX_Files}
    DESTINATION ${TubeTK_INSTALL_INCLUDE_DIR}
    COMPONENT Development )

  TubeTKMacroInstallPlugins( ${PROJECT_NAME} )
endif()

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS ${PROJECT_NAME} )
