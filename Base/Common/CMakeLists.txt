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

find_package( SlicerExecutionModel REQUIRED )
include( ${SlicerExecutionModel_USE_FILE} )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_Common_H_Files
  itktubeDifferenceImageFilter.h
  tubeIndent.h
  tubeMacro.h
  tubeMessage.h
  tubeObject.h
  tubeStringUtilities.h
  tubeTestMain.h )

set( TubeTK_Base_Common_HXX_Files
  itktubeDifferenceImageFilter.hxx )

set( TubeTK_Base_Common_CXX_Files
  tubeIndent.cxx
  tubeObject.cxx )

add_library( TubeTKCommon STATIC
  ${TubeTK_Base_Common_H_Files}
  ${TubeTK_Base_Common_HXX_Files}
  ${TubeTK_Base_Common_CXX_Files} )
target_link_libraries( TubeTKCommon
  ${ITK_LIBRARIES} )
target_include_directories( TubeTKCommon PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

TubeTKMacroInstallPlugins( TubeTKCommon )

SEMMacroBuildCLI(
  NAME "ImageCompareCommand"
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES}
  INCLUDE_DIRECTORIES
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration/ITKRegistrationHelper )

add_executable( JsonCompareCommand JsonCompareCommand.cxx )
target_link_libraries( JsonCompareCommand
  ${ITK_LIBRARIES}
  ITKIOMeta
  ITKMetaIO
  ${JsonCpp_LIBRARIES} )

add_executable( TextCompareCommand TextCompareCommand.cxx )
target_link_libraries( TextCompareCommand ${ITK_LIBRARIES} ITKIOMeta ITKMetaIO)

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_Common_H_Files}
  ${TubeTK_Base_Common_HXX_Files}
  DESTINATION include )

install( TARGETS
  ImageCompareCommand
  JsonCompareCommand
  TextCompareCommand
  DESTINATION bin
  COMPONENT RUNTIME )

TubeTKMacroInstallPlugins( TubeTKCommon )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKCommon )
