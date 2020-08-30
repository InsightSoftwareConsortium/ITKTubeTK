/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "itktubeTubeExtractor.h"

#include "tubeTubeMathFilters.h"

#include <itkImageFileReader.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkSpatialObjectReader.h>

int itktubeTubeExtractorTest( int argc, char * argv[] )
{
  if( argc != 3 )
    {
    std::cout << "itktubeTubeExtractorTest <inputImage> <vessel.tre>"
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
  tubeOp->SetRadiusInObjectSpace( 2.0 );

  typedef itk::SpatialObjectReader<>                   ReaderType;
  typedef itk::SpatialObject<>::ChildrenListType       ObjectListType;
  typedef itk::GroupSpatialObject<>                    GroupType;
  typedef itk::TubeSpatialObject<>                     TubeType;
  typedef TubeType::TubePointListType                  TubePointListType;
  typedef TubeType::TubePointType                      TubePointType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  reader->Update();
  GroupType::Pointer group = reader->GetGroup();

  std::cout << "Number of children = " << group->GetNumberOfChildren()
    << std::endl;

  char tubeName[17];
  std::strcpy( tubeName, "Tube" );
  ObjectListType * tubeList = group->GetChildren( -1, tubeName );

  unsigned int numTubes = tubeList->size();
  std::cout << "Number of tubes = " << numTubes << std::endl;

  typedef itk::Statistics::MersenneTwisterRandomVariateGenerator
    RandGenType;
  RandGenType::Pointer rndGen = RandGenType::New();
  rndGen->Initialize(); // set seed here

  ImageType::IndexType imMinX = tubeOp->GetExtractBoundMinInIndexSpace();
  ImageType::IndexType imMaxX = tubeOp->GetExtractBoundMaxInIndexSpace();
  int margin = 10;

  int failures = 0;
  for( int mcRun=0; mcRun<2; mcRun++ )
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

    TubePointListType tubePointList = tube->GetPoints();
    unsigned int numPoints = tubePointList.size();
    unsigned int rndPointNum = rndGen->GetUniformVariate( 0, 1 )
      * numPoints;
    if( rndPointNum > numPoints-1 )
      {
      rndPointNum = numPoints-1;
      }
    TubePointListType::iterator pntIter = tubePointList.begin();
    for( unsigned int i=0; i<rndPointNum; i++ )
      {
      ++pntIter;
      }
    TubePointType * pnt = static_cast< TubePointType * >( &( *pntIter ) );
    std::cout << "Test point = " << rndPointNum << std::endl;

    tube->Update();
    ImageType::PointType pntX = pnt->GetPositionInWorldSpace();

    TubeOpType::ContinuousIndexType x0;
    bool inMargin = im->TransformPhysicalPointToContinuousIndex( pntX, x0 );
    if( inMargin )
      {
      for( unsigned int i=0; i<ImageType::ImageDimension; i++ )
        {
        if( x0[i] < imMinX[i]+margin || x0[i] > imMaxX[i]-margin )
          {
          inMargin = false;
          break;
          }
        }
      }
    if( !inMargin )
      {
      --mcRun;
      std::cout << "Warning: Tube point outside of image, repicking..."
        << std::endl;
      continue;
      }

    std::cout << "Test index = " << x0 << std::endl;

    std::cout << "Setting min and max z to save time." << std::endl;
    ImageType::IndexType minX = imMinX;
    ImageType::IndexType maxX = imMaxX;
    if( x0[2] - 5 > minX[2] )
      {
      minX[2] = x0[2] - 5;
      }
    if( x0[2] + 5 < maxX[2] )
      {
      maxX[2] = x0[2] + 5;
      }
    tubeOp->SetExtractBoundMinInIndexSpace( minX );
    tubeOp->SetExtractBoundMaxInIndexSpace( maxX );

    if( pnt->GetRadiusInObjectSpace() > 1 )
      {
      tubeOp->SetRadiusInObjectSpace( 0.8 * pnt->GetRadiusInObjectSpace() );
      }
    else
      {
      tubeOp->SetRadiusInObjectSpace( 0.5 );
      }


    ImageType::PointType pntX1 = pntX;
    tubeOp->SetDebug( true );
    tubeOp->GetRidgeExtractor()->SetDebug( true );
    tubeOp->GetRadiusExtractor()->SetDebug( true );
    if( tubeOp->FindLocalTubeInObjectSpace( pntX1 )
      != TubeOpType::RidgeExtractorType::SUCCESS )
      {
      TubeOpType::ContinuousIndexType x1;
      im->TransformPhysicalPointToContinuousIndex( pntX1, x1 );
      std::cout << "Local tube test failed.  No tube found." << std::endl;
      std::cout << "   Source = " << pntX << " (" << x0 << ")" << std::endl;
      std::cout << "   Result = " << pntX1 << " (" << x1 << ")" << std::endl;
      ++failures;
      continue;
      }
    TubeOpType::ContinuousIndexType x1;
    im->TransformPhysicalPointToContinuousIndex( pntX1, x1 );

    double diff = 0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++ )
      {
      double tf = x0[i]-x1[i];
      diff += tf * tf;
      }
    diff = std::sqrt( diff );
    if( diff > 2 )
      {
      std::cout << "Local tube test failed.  Local tube too far."
        << std::endl;
      std::cout << "   Source = " << x0 << std::endl;
      std::cout << "   Result = " << x1 << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Local tube discovery a success!" << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube extraction *****" << std::endl;
    TubeType::Pointer xTube = tubeOp->ExtractTubeInObjectSpace( pntX1, mcRun );
    std::cout << "***** Ending tube extraction *****" << std::endl;

    if( xTube.IsNull() )
      {
      std::cout << "Tube extraction failed" << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting duplicate extraction - should fail *****"
      << std::endl;
    TubeType::Pointer xTube2;
    xTube2 = tubeOp->ExtractTubeInObjectSpace( pntX1, 101 );
    if( xTube2.IsNotNull() )
      {
      std::cout << "Tube extracted twice - test failed" << std::endl;
      std::cout << "  Tube size = " << xTube2->GetPoints().size()
        << std::endl;
      ++failures;
      tubeOp->DeleteTube( xTube2 );
      }

    std::cout << "***** Attempting delete *****"
      << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cout << "Delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting extract after delete *****"
      << std::endl;
    xTube = tubeOp->ExtractTubeInObjectSpace( pntX1, 101 );
    if( xTube.IsNull() )
      {
      std::cout << "Tube extraction after delete failed." << std::endl;
      ++failures;
      continue;
      }

    std::cout << "***** Attempting second delete *****"
      << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cout << "Second delete tube failed" << std::endl;
      ++failures;
      }

    if( xTube.IsNull() )
      {
      continue;
      }

    std::cout << "***** Attempting smooth tube *****" << std::endl;
    ::tube::TubeMathFilters<3> filter;
    filter.SetInputTube( xTube );
    filter.SmoothTube( 5 );

    std::cout << "***** Attempting add tube *****" << std::endl;
    if( !tubeOp->AddTube( xTube ) )
      {
      std::cout << "Add tube failed" << std::endl;
      ++failures;
      }

    std::cout << "***** Attempting delete tube *****" << std::endl;
    if( !tubeOp->DeleteTube( xTube ) )
      {
      std::cout << "Third delete tube failed" << std::endl;
      ++failures;
      }
    }
  delete tubeList;

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
