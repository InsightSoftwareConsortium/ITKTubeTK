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
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "itkTubeTubeExtractor.h"

int itkTubeTubeExtractorTest( int argc, char * argv[] )
  {
  if( argc != 3 )
    {
    std::cout << "itkTubeTubeExtractorTest <inputImage> <vessel.tre>"
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::tube::TubeExtractor<ImageType> TubeOpType;
  TubeOpType::Pointer tubeOp = TubeOpType::New();

  tubeOp->SetInputImage( im );
  tubeOp->SetRadius( 2.0 );

  typedef itk::SpatialObjectReader<>                   ReaderType;
  typedef itk::SpatialObject<>::ChildrenListType       ObjectListType;
  typedef itk::GroupSpatialObject<>                    GroupType;
  typedef itk::VesselTubeSpatialObject<>               TubeType;
  typedef TubeType::PointListType                      PointListType;
  typedef TubeType::PointType                          PointType;
  typedef TubeType::TubePointType                      TubePointType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  GroupType::Pointer group = reader->GetGroup();

  std::cout << "Number of children = " << group->GetNumberOfChildren()
    << std::endl;

  char tubeName[17];
  strcpy( tubeName, "Tube" );
  ObjectListType * tubeList = group->GetChildren( -1, tubeName );

  unsigned int numTubes = tubeList->size();
  std::cout << "Number of tubes = " << numTubes << std::endl;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
    RandGenType;
  RandGenType::Pointer rndGen = RandGenType::New();
  rndGen->Initialize(); // set seed here

  int failures = 0;
  for( unsigned int mcRun=0; mcRun<2; mcRun++ )
    {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube test *****" << std::endl;
    std::cout << "***** Beginning tube test *****" << std::endl;
    unsigned int rndTubeNum = rndGen->GetUniformVariate( 0, 1 ) * numTubes;
    if( rndTubeNum > numTubes-1 )
      {
      rndTubeNum = numTubes-1;
      }
    ObjectListType::iterator tubeIter = tubeList->begin();
    for( unsigned int i=0; i<rndTubeNum; i++ )
      {
      ++tubeIter;
      }
    TubeType::Pointer tube = static_cast< TubeType * >(
      tubeIter->GetPointer() );
    std::cout << "Test tube = " << rndTubeNum << std::endl;

    PointListType tubePointList = tube->GetPoints();
    unsigned int numPoints = tubePointList.size();
    unsigned int rndPointNum = rndGen->GetUniformVariate( 0, 1 )
      * numPoints;
    if( rndPointNum > numPoints-1 )
      {
      rndPointNum = numPoints-1;
      }
    PointListType::iterator pntIter = tubePointList.begin();
    for( unsigned int i=0; i<rndPointNum; i++ )
      {
      ++pntIter;
      }
    TubePointType * pnt = static_cast< TubePointType * >(&(*pntIter));
    std::cout << "Test point = " << rndPointNum << std::endl;

    TubeOpType::ContinuousIndexType x0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++)
      {
      x0[i] = pnt->GetPosition()[i];
      }
    std::cout << "Test index = " << x0 << std::endl;

    TubeOpType::ContinuousIndexType x1 = x0;
    tubeOp->SetDebug( true );
    //tubeOp->GetRidgeOp()->SetDebug( true );
    //tubeOp->GetRadiusOp()->SetDebug( true );
    if( !tubeOp->LocalTube( x1 ) )
      {
      std::cerr << "Local tube test failed.  No tube found." << std::endl;
      std::cerr << "   Source = " << x0 << std::endl;
      std::cerr << "   Result = " << x1 << std::endl;
      ++failures;
      continue;
      }

    double diff = 0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++)
      {
      double tf = x0[i]-x1[i];
      diff += tf * tf;
      }
    diff = vcl_sqrt( diff );
    if( diff > 2 )
      {
      std::cerr << "Local tube test failed.  Local tube too far."
        << std::endl;
      std::cerr << "   Source = " << x0 << std::endl;
      std::cerr << "   Result = " << x1 << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Local tube discovery a success!" << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube extraction *****" << std::endl;
    TubeType::Pointer xTube = tubeOp->ExtractTube( x1, mcRun );
    std::cout << "***** Ending tube extraction *****" << std::endl;

    if( xTube.IsNull() )
      {
      std::cerr << "Tube extraction failed" << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting duplicate extraction - should fail *****"
      << std::endl;
    TubeType::Pointer xTube2;
    xTube2 = tubeOp->ExtractTube( x1, 101 );
    if( xTube2.IsNotNull() )
      {
      std::cerr << "Tube extracted twice - test failed" << std::endl;
      std::cerr << "  Tube size = " << xTube2->GetPoints().size()
        << std::endl;
      ++failures;
      tubeOp->DeleteTube( xTube2 );
      }

    std::cout << "***** Attempting delete *****"
      << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cerr << "Delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting extract after delete *****"
      << std::endl;
    xTube = tubeOp->ExtractTube( x1, 101 );
    if( xTube.IsNull() )
      {
      std::cerr << "Tube extraction after delete failed." << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting second delete *****"
      << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cerr << "Second delete tube failed" << std::endl;
      ++failures;
      }

    if( xTube.IsNull() )
      {
      continue;
      }

    std::cout << "***** Attempting smooth tube *****" << std::endl;
    tubeOp->SmoothTube( xTube, 5 );

    std::cout << "***** Attempting add tube *****" << std::endl;
    if( !tubeOp->AddTube( xTube ) )
      {
      std::cerr << "Add tube failed" << std::endl;
      ++failures;
      }

    std::cout << "***** Attempting delete tube *****" << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cerr << "Third delete tube failed" << std::endl;
      ++failures;
      }
    }

  std::cout << "***** Verifying empty mask *****" << std::endl;
  TubeOpType::TubeMaskImageType::Pointer mask = tubeOp->GetTubeMaskImage();
  itk::ImageRegionIterator< TubeOpType::TubeMaskImageType > maskIt( mask,
    mask->GetLargestPossibleRegion() );
  unsigned int count = 0;
  while( !maskIt.IsAtEnd() )
    {
    if( maskIt.Get() != 0 )
      {
      ++count;
      }
    ++maskIt;
    }
  if( count > 0 )
    {
    std::cout << "Final mask not blank." << std::endl;
    std::cout << "  Number of non-blank = " << count << std::endl;
    ++failures;
    }

  std::cout << "Number of failures = " << failures << std::endl;
  if( failures > 1 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }
