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

add_library( TubeTKLibCommon
  ${TubeTKLib_Base_Common_H_Files}
  ${TubeTKLib_Base_Common_HXX_Files}
  ${TubeTKLib_Base_Common_CXX_Files} )
target_link_libraries( TubeTKLibCommon
  ${ITK_LIBRARIES} )
target_include_directories( TubeTKLibCommon PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> )
include_directories(
  ${ITK_INCLUDE_DIRS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_Common_H_Files}
  ${TubeTKLib_Base_Common_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS TubeTKLibCommon )
