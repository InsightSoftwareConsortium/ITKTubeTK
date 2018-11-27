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

project(TukeTKCompareTools)

add_executable( CompareImages
  CompareImages.cxx )
target_link_libraries( CompareImages
  ${ITK_LIBRARIES} TubeTKLibCommon )
target_include_directories( CompareImages PRIVATE
  ${TubeTKLib_SOURCE_DIR}/Filtering )

install( TARGETS
  CompareImages
  DESTINATION bin
  COMPONENT Runtime )

#-----------------------------------------------------------------------------
set( _toolname CompareTextFiles )
message(STATUS "Configuring executable: ${_toolname}")

add_executable( CompareTextFiles CompareTextFiles.cxx )
target_link_libraries( CompareTextFiles
  ${ITK_LIBRARIES}
  ITKIOMeta
  ITKMetaIO )

install( TARGETS
  CompareTextFiles
  DESTINATION bin
  COMPONENT Runtime )

#-----------------------------------------------------------------------------
if( TubeTKLib_BUILD_TESTING )
  add_subdirectory( Testing )
endif( TubeTKLib_BUILD_TESTING )
