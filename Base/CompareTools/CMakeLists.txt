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

project(TukeTKCompareTools)

if(NOT ITK_SOURCE_DIR)
  find_package( SlicerExecutionModel REQUIRED )
  include( ${SlicerExecutionModel_USE_FILE} )
endif()

if(NOT ITK_SOURCE_DIR)
  find_package( ITK REQUIRED )
  include( ${ITK_USE_FILE} )
endif()

set(_compare_command_targets
  ImageCompareCommand
  TextCompareCommand )

SEMMacroBuildCLI(
  NAME "ImageCompareCommand"
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} TubeTKCommon
  INCLUDE_DIRECTORIES
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration
    ${TubeTK_SLICER_SOURCE_DIR}/Modules/CLI/ExpertAutomatedRegistration/ITKRegistrationHelper )

if( TubeTK_USE_JsonCpp )
  list(APPEND _compare_command_targets
    JsonCompareCommand )
  add_executable( JsonCompareCommand JsonCompareCommand.cxx )
  target_link_libraries( JsonCompareCommand
    ${ITK_LIBRARIES}
    ITKIOMeta
    ITKMetaIO
    ${JsonCpp_LIBRARIES} )
endif( TubeTK_USE_JsonCpp )

add_executable( TextCompareCommand TextCompareCommand.cxx )
target_link_libraries( TextCompareCommand ${ITK_LIBRARIES} ITKIOMeta ITKMetaIO)

install( TARGETS
  ${_compare_command_targets}
  DESTINATION ${TubeTK_INSTALL_BIN_DIR}
  COMPONENT Runtime )
