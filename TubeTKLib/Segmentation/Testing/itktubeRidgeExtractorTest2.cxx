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

#include "itktubeRidgeExtractor.h"

#include <itkImageFileReader.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkSpatialObjectReader.h>

int itktubeRidgeExtractorTest2( int argc, char * argv[] )
{
  if( argc != 3 )
    {
    std::cout
      << "itktubeRidgeExtractorTest <inputImage> <vessel.tre>"
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::tube::RidgeExtractor<ImageType> RidgeOpType;
  RidgeOpType::Pointer ridgeOp = RidgeOpType::New();

  //ridgeOp->SetDebug( true );

  ridgeOp->SetInputImage( im );
  ridgeOp->SetStepX( 0.5 );
  ridgeOp->SetDynamicScale( true );

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

  RidgeOpType::IndexType imMinX = ridgeOp->GetExtractBoundMinInIndexSpace();
  RidgeOpType::IndexType imMaxX = ridgeOp->GetExtractBoundMaxInIndexSpace();
  int margin = 3;
  std::cout << "Bound min = " << imMinX << std::endl;
  std::cout << "Bound max = " << imMaxX << std::endl;

  int failures = 0;
  for( int mcRun=0; mcRun<2; mcRun++ )
    {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube test ***** " << std::endl;
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
    std::cout << "Test position = " << pntX << std::endl;

    RidgeOpType::ContinuousIndexType xContI;
    bool inMargin = im->TransformPhysicalPointToContinuousIndex( pntX, xContI );
    if( inMargin )
      {
      for( unsigned int i=0; i<ImageType::ImageDimension; i++ )
        {
        if( xContI[i] < imMinX[i]+margin || xContI[i] > imMaxX[i]-margin )
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
    std::cout << "Test index = " << xContI << std::endl;

    std::cout << "Setting min and max z to save time." << std::endl;
    RidgeOpType::IndexType minX = imMinX;
    RidgeOpType::IndexType maxX = imMaxX;
    if( xContI[2] - margin > minX[2] )
      {
      minX[2] = xContI[2] - margin;
      }
    if( xContI[2] + margin < maxX[2] )
      {
      maxX[2] = xContI[2] + margin;
      }
    ridgeOp->SetExtractBoundMinInIndexSpace( minX );
    ridgeOp->SetExtractBoundMaxInIndexSpace( maxX );

    if( pnt->GetRadiusInObjectSpace() > 1 )
      {
      ridgeOp->SetScale( 0.8 * pnt->GetRadiusInObjectSpace() );
      }
    else
      {
      ridgeOp->SetScale( 0.5 );
      }


    RidgeOpType::PointType xRidgePnt;
    im->TransformContinuousIndexToPhysicalPoint( xContI, xRidgePnt );
    if( ridgeOp->LocalRidge( xRidgePnt ) != RidgeOpType::SUCCESS )
      {
      RidgeOpType::ContinuousIndexType xRidgeContI;
      im->TransformPhysicalPointToContinuousIndex( xRidgePnt, xRidgeContI );
      std::cout << "*** FAILURE: Local ridge test failed.  No ridge found."
        << std::endl;
      std::cout << "   Source = " << xContI << std::endl;
      std::cout << "   Result = " << xRidgeContI << std::endl;
      ++failures;
      continue;
      }
    RidgeOpType::ContinuousIndexType xRidgeContI;
    im->TransformPhysicalPointToContinuousIndex( xRidgePnt, xRidgeContI );

    double diff = 0;
    for( unsigned int i=0; i<ImageType::ImageDimension; i++ )
      {
      double tf = xContI[i]-xRidgeContI[i];
      diff += tf * tf;
      }
    diff = std::sqrt( diff );
    if( diff > 2*pnt->GetRadiusInObjectSpace() && diff > 4 )
      {
      std::cout << "*** FAILURE: Local ridge test failed.  Local ridge too far."
        << std::endl;
      std::cout << "   Source = " << xContI << std::endl;
      std::cout << "   Result = " << xRidgeContI << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Local ridge discovery a success!" << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube extraction ***** " << std::endl;
    TubeType::Pointer xTube = ridgeOp->ExtractRidge( xRidgePnt, mcRun );
    std::cout << "***** Ending tube extraction ***** " << std::endl;
    std::cout << "   # of points = " << xTube->GetPoints().size() << std::endl;

    if( xTube.IsNull() )
      {
      std::cout << "*** FAILURE: Ridge extraction failed" << std::endl;
      ++failures;
      continue;
      }

    TubeType::Pointer xTube2;
    xTube2 = ridgeOp->ExtractRidge( xRidgePnt, 101 );
    if( xTube2.IsNotNull() )
      {
      std::cout << "*** FAILURE: Ridge extracted twice - test failed" << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "*** FAILURE: Delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Repeat extraction" << std::endl;
    xTube = ridgeOp->ExtractRidge( xRidgePnt, 101 );
    if( xTube.IsNull() )
      {
      std::cout << "*** FAILURE: Ridge extraction after delete failed." << std::endl;
      ++failures;
      continue;
      }
    std::cout << "   # of points = " << xTube->GetPoints().size() << std::endl;

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "*** FAILURE: Second delete tube failed" << std::endl;
      ++failures;
      continue;
      }

    if( xTube.IsNull() )
      {
      std::cout << "*** FAILURE: Ridge extraction failed" << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->AddTube( xTube ) )
      {
      std::cout << "*** FAILURE: Add tube failed" << std::endl;
      ++failures;
      continue;
      }

    if( !ridgeOp->DeleteTube( xTube ) )
      {
      std::cout << "*** FAILURE: Third delete tube failed" << std::endl;
      ++failures;
      continue;
      }
    }
  delete tubeList;

  RidgeOpType::TubeMaskImageType::Pointer mask =
    ridgeOp->GetTubeMaskImage();
  itk::ImageRegionIterator< RidgeOpType::TubeMaskImageType > maskIt( mask,
    mask->GetLargestPossibleRegion() );
  while( !maskIt.IsAtEnd() )
    {
    if( maskIt.Get() != 0 )
      {
      std::cout << "Final mask not blank." << std::endl;
      std::cout << "Number of failures = " << failures << std::endl;
      return EXIT_FAILURE;
      }
    ++maskIt;
    }

  std::cout << "Ridge termination code counts:" << std::endl;
  for( unsigned int code = 0; code < ridgeOp->GetNumberOfFailureCodes();
    ++code )
    {
    std::cout << "   " << ridgeOp->GetFailureCodeName(
      RidgeOpType::FailureCodeEnum( code ) ) << " : "
      << ridgeOp->GetFailureCodeCount( RidgeOpType::FailureCodeEnum(
      code ) ) << std::endl;
    }

  std::cout << "Number of failures = " << failures << std::endl;
  if( failures > 0 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
