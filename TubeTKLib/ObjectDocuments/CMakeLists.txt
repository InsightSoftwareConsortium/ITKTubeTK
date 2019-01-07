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

set( TubeTKLib_Base_ObjectDocuments_H_Files
  itktubeBlobSpatialObjectDocument.h
  itktubeDocument.h
  itktubeImageDocument.h
  itktubeObjectDocument.h
  itktubeObjectDocumentToImageFilter.h
  itktubeObjectDocumentToObjectSource.h
  itktubeSpatialObjectDocument.h
  tubeMetaDocument.h
  tubeMetaObjectDocument.h
  tubeOptionList.h )

set( TubeTk_Base_ObjectDocuments_HXX_Files
  itktubeObjectDocumentToImageFilter.hxx
  itktubeObjectDocumentToObjectSource.hxx )

set( TubeTk_Base_ObjectDocuments_SRCS
  tubeMetaDocument.cxx
  tubeMetaObjectDocument.cxx
  tubeOptionList.cxx )

add_library( TubeTKLibObjectDocuments STATIC
  ${TubeTKLib_Base_ObjectDocuments_H_Files}
  ${TubeTk_Base_ObjectDocuments_HXX_Files}
  ${TubeTk_Base_ObjectDocuments_SRCS} )

target_link_libraries( TubeTKLibObjectDocuments PUBLIC
  TubeTKLibCommon
  ${ITK_LIBRARIES}
  ITKIOMeta
  ITKIOSpatialObjects )

target_include_directories( TubeTKLibObjectDocuments PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTKLib_Base_ObjectDocuments_H_Files}
  ${TubeTKLib_Base_ObjectDocuments_HXX_Files}
  DESTINATION include
  COMPONENT Development )

set_property( GLOBAL APPEND PROPERTY TubeTKLib_TARGETS 
  TubeTKLibObjectDocuments )
