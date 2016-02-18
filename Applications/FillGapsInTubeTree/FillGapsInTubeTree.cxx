/*=========================================================================
   Library:   TubeTK

   Copyright 2010 Kitware Inc. 28 Corporate Drive,
   Clifton Park, NY, 12065, USA.

   All rights reserved.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
=========================================================================*/
#include <iostream>
#include <sstream>

#include "metaScene.h"

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "FillGapsInTubeTreeCLP.h"

#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkNumericTraits.h"

#include <itkTimeProbesCollectorBase.h>

template< unsigned int VDimension >
int DoIt( int argc, char * argv[] );

template< unsigned int VDimension >
int FillGap( typename itk::GroupSpatialObject< VDimension >::Pointer &pTubeGroup,
             std::string InterpolationMethod )
{
  //typedefs
  typedef itk::GroupSpatialObject< VDimension >         TubeGroupType;
  typedef typename TubeGroupType::ChildrenListPointer   TubeListPointerType;
  typedef itk::VesselTubeSpatialObject< VDimension >    TubeType;
  typedef typename TubeType::Pointer                    TubePointerType;
  typedef typename TubeType::TubePointType              TubePointType;
  typedef typename TubeType::PointType                  PositionType;
  typedef itk::IndexValueType                           TubeIdType;
  typedef typename TubeType::PointListType              TubePointListType;

  char tubeName[] = "Tube";
  TubeListPointerType pTubeList
    = pTubeGroup->GetChildren(
    pTubeGroup->GetMaximumDepth(), tubeName );

  for( typename TubeGroupType::ChildrenListType::iterator
    itSourceTubes = pTubeList->begin();
    itSourceTubes != pTubeList->end(); ++itSourceTubes )
    {
    TubePointerType pCurTube
      = dynamic_cast< TubeType * >( itSourceTubes->GetPointer() );
    TubeIdType curTubeId = pCurTube->GetId();
    TubeIdType curParentTubeId = pCurTube->GetParentId();
    TubePointType* parentNearestPoint;
    if( pCurTube->GetRoot() == false &&  curParentTubeId != pTubeGroup->GetId() )
      {
      //find parent target tube
      for( typename TubeGroupType::ChildrenListType::iterator
        itTubes = pTubeList->begin();
        itTubes != pTubeList->end(); ++itTubes )
        {
        TubePointerType pTube
          = dynamic_cast< TubeType * >( itTubes->GetPointer() );
        if( pTube->GetId() == curParentTubeId )
          {
          double minDistance = itk::NumericTraits<double>::max();
          int flag =-1;
          for ( int index = 0; index < pTube->GetNumberOfPoints(); index++ )
            {
            TubePointType* tubePoint =
            dynamic_cast< TubePointType* >( pTube->GetPoint( index ) );
            PositionType tubePointPosition = tubePoint->GetPosition();
            double distance = tubePointPosition.SquaredEuclideanDistanceTo
              ( pCurTube->GetPoint( 0 )->GetPosition() );
            if ( minDistance > distance )
              {
              minDistance = distance;
              parentNearestPoint = tubePoint;
              flag = 1;
              }
            distance = tubePointPosition.SquaredEuclideanDistanceTo
              ( pCurTube->GetPoint( pCurTube->GetNumberOfPoints() - 1 )->GetPosition() );
            if ( minDistance > distance )
              {
              minDistance = distance;
              parentNearestPoint = tubePoint;
              flag = 2;
              }
            }

          if( InterpolationMethod == "Straight_Line" )
            {
            if( flag == 1 )//add as the starting point
              {
              TubePointListType targetTubePoints = pCurTube->GetPoints();
              pCurTube->Clear();
              pCurTube->GetPoints().push_back(*parentNearestPoint);
              for( int i=0; i < targetTubePoints.size(); i++ )
                {
                pCurTube->GetPoints().push_back( targetTubePoints[i] );
                }
              }
            if( flag == 2 )// add as an end point
              {
              pCurTube->GetPoints().push_back( *parentNearestPoint );
              }
            }
          break;
          }
        }
      }
    }
  return EXIT_SUCCESS;
}

template< unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // Ensure that the input image dimension is valid
  // We only support 2D and 3D Images due to the
  // limitation of itkTubeSpatialObject
  if ( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }
  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "FillGapsInTubeTree",
    CLPProcessInformation );
  progressReporter.Start();
  float progress = 0;

  timeCollector.Start( "Loading Input TRE File" );
  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();

  try
    {
    tubeFileReader->SetFileName( inputTree.c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Loading Input TRE File" );
  progress = 0.25;

  timeCollector.Start( "Filling gaps in input tree" );
  tubeStandardOutputMacro( << "\n>> Filling gaps in input tree" );

  typename itk::GroupSpatialObject< VDimension >::Pointer inputTubes;
  inputTubes = tubeFileReader->GetGroup();
  FillGap< VDimension >( inputTubes, InterpolationMethod );

  timeCollector.Stop( "Filling gaps in input tree" );
  progress = 0.75;

  timeCollector.Start( "Writing the updates TRE file." );

  typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();

  try
    {
    tubeWriter->SetFileName( outputTree.c_str() );
    tubeWriter->SetInput( tubeFileReader->GetGroup() );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Writing the updates TRE file." );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch ( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

  MetaScene *mScene = new MetaScene;
  mScene->Read( inputTree.c_str() );
  if( mScene->GetObjectList()->empty() )
    {
    tubeWarningMacro( << "Input TRE file has no spatial objects" );
    return EXIT_SUCCESS;
    }

  switch( mScene->GetObjectList()->front()->NDims() )
    {
    case 2:
      return DoIt<2>( argc, argv );
      break;

    case 3:
      return DoIt<3>( argc, argv );
      break;

    default:
      tubeErrorMacro(<< "Error: Only 2D and 3D data is currently supported.");
      return EXIT_FAILURE;
    }
}
