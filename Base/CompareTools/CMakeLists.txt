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

#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME CompareImages
  LOGO_HEADER ${TubeTK_SOURCE_DIR}/Base/CLI/TubeTKLogo.h
  TARGET_LIBRARIES ${ITK_LIBRARIES} TubeTKCommon
  INCLUDE_DIRECTORIES
    ${TubeTK_SOURCE_DIR}/Base/Filtering )

#-----------------------------------------------------------------------------
if( TubeTK_USE_JSONCPP )

  set( _toolname CompareJsonFiles )

  message(STATUS "Configuring executable: ${_toolname}")
  add_executable( ${_toolname} ${_toolname}.cxx )
  target_link_libraries( ${_toolname}
    ${ITK_LIBRARIES}
    ITKIOMeta
    ITKMetaIO
    ${JsonCpp_LIBRARIES} )
  install( TARGETS
    ${_toolname}
    DESTINATION ${TubeTK_INSTALL_BIN_DIR}
    COMPONENT Runtime )

endif( TubeTK_USE_JSONCPP )

#-----------------------------------------------------------------------------
set( _toolname CompareTextFiles )
message(STATUS "Configuring executable: ${_toolname}")

add_executable( ${_toolname} ${_toolname}.cxx )
target_link_libraries( ${_toolname} ${ITK_LIBRARIES} ITKIOMeta ITKMetaIO)
install( TARGETS
  ${_toolname}
  DESTINATION ${TubeTK_INSTALL_BIN_DIR}
  COMPONENT Runtime )

#-----------------------------------------------------------------------------
if( TubeTK_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTK_BUILD_TESTING )
