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

add_executable( CompareImages
  CompareImages.cxx )
target_link_libraries( CompareImages
  ${ITK_LIBRARIES}
  TubeTKLibCommon )
target_include_directories( CompareImages PRIVATE
  ${ITK_INCLUDE_DIRS}
  ${TubeTK_SOURCE_DIR}/TubeTKLib/Filtering )

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
target_include_directories( CompareTextFiles PRIVATE
  ${ITK_INCLUDE_DIRS} )

install( TARGETS
  CompareTextFiles
  DESTINATION bin
  COMPONENT Runtime )

#-----------------------------------------------------------------------------
if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )
