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

project( TubeTKObjectDocuments )

set( TubeTK_Base_ObjectDocuments_H_Files
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

add_library( TubeTKObjectDocuments STATIC
  ${TubeTK_Base_ObjectDocuments_H_Files}
  ${TubeTk_Base_ObjectDocuments_HXX_Files}
  ${TubeTk_Base_ObjectDocuments_SRCS} )

target_link_libraries( TubeTKObjectDocuments
  TubeTKCommon ${ITK_LIBRARIES} ITKIOMeta ITKIOSpatialObjects )

target_include_directories( TubeTKObjectDocuments PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_ObjectDocuments_H_Files}
  ${TubeTK_Base_ObjectDocuments_HXX_Files}
  DESTINATION include )

TubeTKMacroInstallPlugins( TubeTKObjectDocuments )

set_property( GLOBAL APPEND PROPERTY TubeTK_TARGETS TubeTKObjectDocuments )
