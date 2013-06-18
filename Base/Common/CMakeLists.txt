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

project( TubeCommon )

find_package( ITK REQUIRED )
include( ${ITK_USE_FILE} )

set( TubeTK_Base_Common_H_Files
  itktubeDifferenceImageFilter2.h
  itktubeFilterWatcher.h
  tubeMacro.h
  tubeMessage.h
  tubeStringUtilities.h
  tubeTestMain.h )

set( TubeTK_Base_Common_HXX_Files
  itktubeDifferenceImageFilter2.hxx )

add_custom_target( TubeCommon SOURCES
  ${TubeTK_Base_Common_H_Files}
  ${TubeTK_Base_Common_HXX_Files} )

add_executable( ImageCompareCommand ImageCompareCommand.cxx )
target_link_libraries( ImageCompareCommand ${ITK_LIBRARIES} )

add_executable( JsonCompareCommand JsonCompareCommand.cxx )
target_link_libraries( JsonCompareCommand
  ${ITK_LIBRARIES}
  ${JsonCpp_LIBRARIES} )

add_executable( TextCompareCommand TextCompareCommand.cxx )
target_link_libraries( TextCompareCommand ${ITK_LIBRARIES} )

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
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  COMPONENT RUNTIME )
