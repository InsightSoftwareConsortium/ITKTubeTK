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

#include "../itkRadiusExtractor.h"

#include "tubeMessage.h"

int itkRadiusExtractorTest( int argc, char * argv[] )
  {
  if( argc != 3 )
    {
    std::cout 
      << "itkRadiusExtractorTest <inputImage> <vessel.tre>" 
      << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;

  typedef itk::ImageFileReader< ImageType > ImageReaderType;
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName( argv[1] );
  imReader->Update();

  ImageType::Pointer im = imReader->GetOutput();

  typedef itk::RadiusExtractor<ImageType> RadiusOpType;
  RadiusOpType::Pointer radiusOp = RadiusOpType::New();

  radiusOp->SetInputImage( im );
  radiusOp->SetScale( 2.0 );
  radiusOp->SetExtent( 3.0 );

  bool returnStatus = EXIT_SUCCESS;

  double dataMin = radiusOp->GetDataMin();
  if( dataMin != 54 )
    {
    tube::ErrorMessage( "Data min != 54" );
    returnStatus = EXIT_FAILURE;
    }

  double dataMax = radiusOp->GetDataMax();
  if( dataMax != 230 )
    {
    tube::ErrorMessage( "Data max != 230" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetRadiusMin( 0.5 );
  if( radiusOp->GetRadiusMin() != 0.5 )
    {
    tube::ErrorMessage( "Radius min != 0.5" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetRadiusMax( 10.0 );
  if( radiusOp->GetRadiusMax() != 10.0 )
    {
    tube::ErrorMessage( "Radius max != 10.0" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetRadius0( 2.0 );
  if( radiusOp->GetRadius0() != 2.0 )
    {
    tube::ErrorMessage( "Radius 0 != 2.0" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetThreshMedialness( 0.05 );
  if( radiusOp->GetThreshMedialness() != 0.05 )
    {
    tube::ErrorMessage( "ThreshMedialness != 0.05" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetThreshMedialnessStart( 0.02 );
  if( radiusOp->GetThreshMedialnessStart() != 0.02 )
    {
    tube::ErrorMessage( "ThreshMedialnessStart != 0.02" );
    returnStatus = EXIT_FAILURE;
    }

  radiusOp->SetExtractRidge( true );
  if( radiusOp->GetExtractRidge() != true )
    {
    tube::ErrorMessage( "ExtractRidge != true" );
    returnStatus = EXIT_FAILURE;
    }

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
  for( unsigned int mcRun=0; mcRun<10; mcRun++ )
    {
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "***** Beginning tube test ***** " << std::endl;
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
  
    double r0 = pnt->GetRadius();
    double r1 = r0;

    radiusOp->SetDebug( true );
    if( !radiusOp->ComputeOptimalRadiusAtPoint( *pnt, r1, r1/4, r1*2,
        r1/8, r1/16 ) )
      {
      std::cerr << "ComputeOptimalRadius test failed." << std::endl;
      std::cerr << "   Source = " << pnt->GetPosition() 
        << " at r = " << r0 << std::endl;
      std::cerr << "      Result r = " << r1 << std::endl;
      ++failures;
      continue;
      }
  
    double diff = vnl_math_abs( r1 - r0 );
    if( diff > 0.5 * r0 )
      {
      std::cerr << "Radius estimate inaccurate." << std::endl;
      std::cerr << "   Source = " << r0 << std::endl;
      std::cerr << "      Result = " << r1 << std::endl;
      ++failures;
      continue;
      }

    std::cout << "Radius estimate a success!" << std::endl;
    std::cout << std::endl;

    std::cout << "***** Beginning radius estimation for entire tube  ***** "
      << std::endl;
    if( ! radiusOp->ComputeTubeRadii( *tube ) )
      {
      std::cerr << "Tube radius estimation failed." << std::endl;
      ++failures;
      continue;
      }

    // Check results HERE 
    //

    std::cout << "***** Ending tube extraction ***** " << std::endl;
    }

  std::cout << "Number of failures = " << failures << std::endl;
  if( failures > 10 )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
  }

