/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

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

int tubeBaseObjectDocumentsPrintTest( int tubeNotUsed( argc ),
                                      char * tubeNotUsed( argv )[] )
{
  typedef itk::tube::BlobSpatialObjectDocument BlobSpatialObjectDocumentType;

  BlobSpatialObjectDocumentType::Pointer blobSpatialObjectDocument
    = BlobSpatialObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::BlobSpatialObjectDocument"
                           << blobSpatialObjectDocument );

  typedef itk::tube::Document DocumentType;

  DocumentType::Pointer document = DocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::Document" << document );

  typedef itk::tube::ImageDocument ImageDocumentType;

  ImageDocumentType::Pointer imageDocument = ImageDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ImageDocument"
                           << imageDocument );

  typedef tube::MetaDocument MetaDocumentType;

  MetaDocumentType::Pointer metaDocument = new MetaDocumentType();
  tubeStandardOutputMacro( << "-------------tube::MetaDocument"
                           << metaDocument );
  delete metaDocument;

  typedef tube::MetaObjectDocument MetaObjectDocumentType;

  MetaObjectDocumentType::Pointer metaObjectDocument
    = new MetaObjectDocumentType();
  tubeStandardOutputMacro( << "-------------tube::MetaObjectDocument"
                           << metaObjectDocument );
  delete metaObjectDocument;

  typedef itk::tube::ObjectDocument ObjectDocumentType;

  ObjectDocumentType::Pointer objectDocument = ObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocument"
                           << objectDocument );

  typedef itk::Image< float, 3 > ImageType;
  typedef itk::tube::ObjectDocumentToImageFilter< ObjectDocumentType, ImageType >
    ObjectDocumentToImageFilterType;

  ObjectDocumentToImageFilterType::Pointer objectDocumentToImageFilter
    = ObjectDocumentToImageFilterType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocumentToImageFilter"
                           << objectDocumentToImageFilter );

  typedef itk::tube::ObjectDocumentToObjectSource< ObjectDocumentType, 3 >
    ObjectDocumentToObjectSourceType;

  ObjectDocumentToObjectSourceType::Pointer objectDocumentToObjectSource
    = ObjectDocumentToObjectSourceType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::ObjectDocumentToObjectSource"
                           << objectDocumentToObjectSource );

  typedef tube::OptionList OptionListType;

  OptionListType::Pointer optionList = new OptionListType();
  tubeStandardOutputMacro( << "-------------tube::OptionList" << optionList );
  delete optionList;

  typedef itk::tube::SpatialObjectDocument SpatialObjectDocumentType;

  SpatialObjectDocumentType::Pointer spatialObjectDocument
    = SpatialObjectDocumentType::New();
  tubeStandardOutputMacro( << "-------------itk::tube::SpatialObjectDocument"
                           << spatialObjectDocument );

  return EXIT_SUCCESS;
}
