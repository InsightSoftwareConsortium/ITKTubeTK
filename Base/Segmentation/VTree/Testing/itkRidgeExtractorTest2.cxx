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
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSpatialObjectReader.h"
#include "itkGroupSpatialObject.h"
#include "itkImageRegionIteratorWithIndex.h"

#include "../itkRidgeExtractor.h"

int itkRidgeExtractorTest2( int argc, char * argv[] )
  {
  if( argc != 3 )
    {
    std::cout 
      << "itkRidgeExtractorTest <inputImage> <vessel.tre>" 
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::RidgeExtractor<ImageType> RidgeOpType;
  RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  ridgeOp->SetInputImage( im );
  ridgeOp->SetStepX( 0.75 );
  ridgeOp->SetScale( 2.0 );
  ridgeOp->SetExtent( 3.0 );
  ridgeOp->SetDynamicScale( true );

  typedef itk::SpatialObjectReader<> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  typedef itk::GroupSpatialObject<> GroupType;
  GroupType::Pointer group = reader->GetGroup();

  return EXIT_SUCCESS;
  }

