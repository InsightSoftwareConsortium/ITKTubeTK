/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "itktubeBlobSpatialObjectDocument.h"
#include "itktubeDocument.h"
#include "itktubeImageDocument.h"
#include "itktubeObjectDocument.h"
#include "itktubeObjectDocumentToImageFilter.h"
#include "itktubeObjectDocumentToObjectSource.h"
#include "itktubeSpatialObjectDocument.h"
#include "tubeMetaDocument.h"
#include "tubeMetaObjectDocument.h"
#include "tubeOptionList.h"

int tubeObjectDocumentsPrintTest( int tubeNotUsed( argc ),
                                      char * tubeNotUsed( argv )[] )
{
  using BlobSpatialObjectDocumentType = itk::tube::BlobSpatialObjectDocument;

  BlobSpatialObjectDocumentType::Pointer blobSpatialObjectDocument
    = BlobSpatialObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::BlobSpatialObjectDocument"
                           << blobSpatialObjectDocument );

  using DocumentType = itk::tube::Document;

  DocumentType::Pointer document = DocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::Document" << document );

  using ImageDocumentType = itk::tube::ImageDocument;

  ImageDocumentType::Pointer imageDocument = ImageDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ImageDocument"
                           << imageDocument );

  using MetaDocumentType = tube::MetaDocument;

  MetaDocumentType::Pointer metaDocument = new MetaDocumentType();
  tubeStandardOutputMacro( << "-------------tube::MetaDocument"
                           << metaDocument );
  delete metaDocument;

  using MetaObjectDocumentType = tube::MetaObjectDocument;

  MetaObjectDocumentType::Pointer metaObjectDocument
    = new MetaObjectDocumentType();
  tubeStandardOutputMacro( << "-------------tube::MetaObjectDocument"
                           << metaObjectDocument );
  delete metaObjectDocument;

  using ObjectDocumentType = itk::tube::ObjectDocument;

  ObjectDocumentType::Pointer objectDocument = ObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocument"
                           << objectDocument );

  using ImageType = itk::Image< float, 3 >;
  using ObjectDocumentToImageFilterType = itk::tube::ObjectDocumentToImageFilter< ObjectDocumentType, ImageType >;

  ObjectDocumentToImageFilterType::Pointer objectDocumentToImageFilter
    = ObjectDocumentToImageFilterType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocumentToImageFilter"
                           << objectDocumentToImageFilter );

  using ObjectDocumentToObjectSourceType = itk::tube::ObjectDocumentToObjectSource< ObjectDocumentType, 3 >;

  ObjectDocumentToObjectSourceType::Pointer objectDocumentToObjectSource
    = ObjectDocumentToObjectSourceType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocumentToObjectSource"
                           << objectDocumentToObjectSource );

  using OptionListType = tube::OptionList;

  OptionListType::Pointer optionList = new OptionListType();
  tubeStandardOutputMacro( << "-------------tube::OptionList" << optionList );
  delete optionList;

  using SpatialObjectDocumentType = itk::tube::SpatialObjectDocument;

  SpatialObjectDocumentType::Pointer spatialObjectDocument
    = SpatialObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::SpatialObjectDocument"
                           << spatialObjectDocument );

  return EXIT_SUCCESS;
}
