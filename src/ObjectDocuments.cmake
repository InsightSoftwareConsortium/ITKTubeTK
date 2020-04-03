##############################################################################
#
# Library:   TubeTK
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

set( TubeTK_ObjectDocuments_H_Files
  ObjectDocuments/itktubeBlobSpatialObjectDocument.h
  ObjectDocuments/itktubeDocument.h
  ObjectDocuments/itktubeImageDocument.h
  ObjectDocuments/itktubeObjectDocument.h
  ObjectDocuments/itktubeObjectDocumentToImageFilter.h
  ObjectDocuments/itktubeObjectDocumentToObjectSource.h
  ObjectDocuments/itktubeSpatialObjectDocument.h
  ObjectDocuments/tubeMetaDocument.h
  ObjectDocuments/tubeMetaObjectDocument.h
  ObjectDocuments/tubeOptionList.h )

set( TubeTk_ObjectDocuments_HXX_Files
  ObjectDocuments/itktubeObjectDocumentToImageFilter.hxx
  ObjectDocuments/itktubeObjectDocumentToObjectSource.hxx )

set( TubeTk_ObjectDocuments_CXX_Files
  ObjectDocuments/tubeMetaDocument.cxx
  ObjectDocuments/tubeMetaObjectDocument.cxx
  ObjectDocuments/tubeOptionList.cxx )

list( APPEND TubeTK_SRCS
  ${TubeTK_Objectdocuments_H_Files}
  ${TubeTK_Objectdocuments_HXX_Files}
  ${TubeTK_Objectdocuments_CXX_Files} )
