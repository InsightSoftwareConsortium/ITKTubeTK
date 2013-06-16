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

project( TubeObjectDocuments )

# Find ITK
find_package( ITK REQUIRED )
include( ${USE_ITK_FILE} )

set( TubeTK_Base_ObjectDocuments_H_Files
  itkBlobSpatialObjectDocument.h
  itkDocument.h
  itkImageDocument.h
  itkObjectDocument.h
  itkObjectDocumentToImageFilter.h
  itkObjectDocumentToObjectSource.h
  itkSpatialObjectDocument.h
  itkSpatialObjectToSpatialObjectFilter.h
  tubeMetaDocument.h
  tubeMetaObjectDocument.h
  tubeOptionList.h )

set( TubeTk_Base_ObjectDocuments_TXX_Files
  itkObjectDocumentToImageFilter.txx
  itkObjectDocumentToObjectSource.txx
  itkSpatialObjectToSpatialObjectFilter.txx )

set( TubeTk_Base_ObjectDocuments_SRCS
  tubeMetaDocument.cxx
  tubeMetaObjectDocument.cxx
  tubeOptionList.cxx )

add_library( TubeObjectDocuments STATIC
  ${TubeTK_Base_ObjectDocuments_H_Files}
  ${TubeTk_Base_ObjectDocuments_TXX_Files}
  ${TubeTk_Base_ObjectDocuments_SRCS} )

target_link_libraries( TubeObjectDocuments ${ITK_LIBRARIES} )

slicer4_set_plugins_output_path( TubeObjectDocuments )
set( TARGETS TubeObjectDocuments )
slicer4_install_plugins( ${TARGETS} )

if( BUILD_TESTING )
  add_subdirectory( Testing )
endif( BUILD_TESTING )

install( FILES
  ${TubeTK_Base_ObjectDocuments_H_Files}
  ${TubeTK_Base_ObjectDocuments_TXX_Files}
  DESTINATION include )
