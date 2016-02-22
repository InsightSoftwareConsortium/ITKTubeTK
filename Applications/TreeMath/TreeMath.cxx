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

#include "tubeMessage.h"
#include "tubeMacro.h"
#include "TreeMathCLP.h"

#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkNumericTraits.h"

#include <metaCommand.h>


template< unsigned int VDimension >
int DoIt( int argc, char * argv[] );

template< unsigned int VDimension >
void InterpolatePath(
  typename itk::VesselTubeSpatialObject< VDimension >::TubePointType *parentNearestPoint,
  typename itk::VesselTubeSpatialObject< VDimension >::TubePointType *childEndPoint,
  typename itk::VesselTubeSpatialObject< VDimension >::PointListType &newTubePoints,
  char InterpolationMethod )
{
  if( InterpolationMethod == 'S' )
    {
    newTubePoints.push_back( *parentNearestPoint );
    }
  return;
}

template< unsigned int VDimension >
void FillGap( typename itk::GroupSpatialObject< VDimension >::Pointer &pTubeGroup,
             char InterpolationMethod )
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

          TubePointListType newTubePoints;
          if( flag == 1 )
            {
            TubePointType* childTubeStartPoint =
              dynamic_cast< TubePointType* >( pCurTube->GetPoint( 0 ) );
            InterpolatePath< VDimension >(parentNearestPoint,
              childTubeStartPoint, newTubePoints, InterpolationMethod );
            TubePointListType targetTubePoints = pCurTube->GetPoints();
            pCurTube->Clear();
            for( int index = 0; index < newTubePoints.size(); index++ )
              {
              pCurTube->GetPoints().push_back( newTubePoints[ index ] );
              }
            for( int i = 0; i < targetTubePoints.size(); i++ )
              {
              pCurTube->GetPoints().push_back( targetTubePoints[ i ] );
              }
            }
          if( flag == 2 )
            {
            TubePointType* childTubeEndPoint =
              dynamic_cast< TubePointType* >
              ( pCurTube->GetPoint( pCurTube->GetNumberOfPoints() - 1 ) );
            InterpolatePath< VDimension >(
              parentNearestPoint, childTubeEndPoint, newTubePoints, InterpolationMethod );
            for( int index = newTubePoints.size() - 1; index >= 0; index-- )
              {
              pCurTube->GetPoints().push_back( newTubePoints[ index ] );
              }
            }
          break;
          }
        }
      }
    }
}

template< unsigned int VDimension >
int DoIt( MetaCommand & command )
{
  if ( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  typename itk::GroupSpatialObject< VDimension >::Pointer inputTubes;

  MetaCommand::OptionVector parsed = command.GetParsedOptions();

  // Load TRE File
  tubeStandardOutputMacro( << "\n>> Loading TRE File" );

  typedef itk::SpatialObjectReader< VDimension > TubesReaderType;
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();
  try
    {
    tubeFileReader->SetFileName( command.GetValueAsString(
    "infile" ).c_str() );
    tubeFileReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error loading TRE File: "
      + std::string( err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  inputTubes = tubeFileReader->GetGroup();
  inputTubes->ComputeObjectToWorldTransform();

  MetaCommand::OptionVector::const_iterator it = parsed.begin();
  while( it != parsed.end() )
    {
    if( it->name == "Write" )
      {
      tubeStandardOutputMacro( << "\n>> Writing TRE File" );

      typedef itk::SpatialObjectWriter< VDimension > TubeWriterType;
      typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
      try
        {
        tubeWriter->SetFileName( command.GetValueAsString(
        *it, "filename" ).c_str() );
        tubeWriter->SetInput( inputTubes );
        tubeWriter->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        tube::ErrorMessage( "Error writing TRE file: "
          + std::string( err.GetDescription() ) );
        return EXIT_FAILURE;
        }
      }
    else if( it->name == "FillGapsInTubeTree" )
      {
      tubeStandardOutputMacro( << "\n>> Filling gaps in input tree" );
      FillGap< VDimension >( inputTubes, command.GetValueAsString
        ( *it, "InterpolationMethod" ).c_str()[0] );
      }
    ++it;
    }
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  //PARSE_ARGS;
  MetaCommand command;

  command.SetName( "TreeMath" );
  command.SetVersion( "1.0" );
  command.SetAuthor( "Sumedha Singla" );
  command.SetDescription( "Perform several filters on a tube tree." );

  command.AddField( "infile", "infile filename",
  MetaCommand::STRING, MetaCommand::DATA_IN );

  command.SetOption( "Write", "w", false,
    "Writes current tubes to the designated file." );
  command.AddOptionField( "Write", "filename", MetaCommand::STRING, true,
    "", "Output filename", MetaCommand::DATA_OUT );

  command.SetOption( "FillGapsInTubeTree", "f", false,
    "Connects the parent and child tube if they have a gap inbetween,"
    " by interpolating the path inbetween." );
  command.AddOptionField( "FillGapsInTubeTree", "InterpolationMethod",
    MetaCommand::STRING, true, "",
    "[S]traight Line, [L]Linear Interpolation, [C]urve Fitting, [M]inimal Path",
    MetaCommand::DATA_IN );

  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  return DoIt< 3 >( command );
}
