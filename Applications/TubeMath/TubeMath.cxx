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
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include "tubeTubeMath.h"

#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkVesselTubeSpatialObject.h>

#include <itkImageFileReader.h>
#include <itkDisplacementFieldTransform.h>
#include <itkTransformFactory.h>
#include <itkTransformFileReader.h>

#include <metaCommand.h>

#include "TubeMathCLP.h"

#define ALL_TUBES_CURRENT -98989

template< unsigned int DimensionT >
typename itk::GroupSpatialObject< DimensionT >::Pointer
ReadTubeFile( const char * fileName )
{
  typedef itk::SpatialObjectReader< DimensionT > SpatialObjectReaderType;

  typename SpatialObjectReaderType::Pointer reader =
    SpatialObjectReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  typename SpatialObjectReaderType::GroupType::Pointer group =
    reader->GetGroup();
  group->ComputeObjectToWorldTransform();
  return group;
}

template< unsigned int DimensionT >
void WriteTubeFile( typename itk::GroupSpatialObject< DimensionT >::Pointer
  object, const char * fileName )
{
  typedef itk::SpatialObjectWriter< DimensionT > SpatialObjectWriterType;

  typename SpatialObjectWriterType::Pointer writer =
    SpatialObjectWriterType::New();
  writer->SetInput( object );
  writer->SetFileName( fileName );
  writer->Update();
}

template< class PixelT, unsigned int DimensionT >
typename itk::Image< PixelT, DimensionT >::Pointer
ReadImageFile( const char * fileName )
{
  typedef itk::Image< PixelT, DimensionT >  ImageType;
  typedef itk::ImageFileReader< ImageType > ImageReaderType;

  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( fileName );
  reader->Update();

  return reader->GetOutput();
}

template< class PixelT, unsigned int DimensionT >
void
SetPropertyFromImage( typename itk::GroupSpatialObject< DimensionT >::
  Pointer & inputTubes, int currentTube, typename itk::Image< PixelT,
  DimensionT>::Pointer & inputImage, char propertyId )
{
  typedef itk::VesselTubeSpatialObject< DimensionT >  TubeType;
  typedef typename TubeType::PointListType            TubePointListType;
  typedef typename TubeType::TubePointType            TubePointType;
  typedef itk::Image< PixelT, DimensionT >            ImageType;

  typename TubeType::ChildrenListType::iterator tubeIterator;

  char soTypeName[80];
  strcpy( soTypeName, "VesselTubeSpatialObject" );
  typename TubeType::ChildrenListPointer inputTubeList =
    inputTubes->GetChildren( inputTubes->GetMaximumDepth(), soTypeName );
  for( tubeIterator = inputTubeList->begin(); tubeIterator !=
    inputTubeList->end(); tubeIterator++ )
    {
    typename TubeType::Pointer inputTube = ((TubeType *)( tubeIterator->
      GetPointer()));

    if( currentTube == ALL_TUBES_CURRENT || inputTube->GetId() ==
      currentTube )
      {
      inputTube->ComputeObjectToWorldTransform();

      unsigned int pointListSize = inputTube->GetNumberOfPoints();
      for( unsigned int pointNum = 0; pointNum < pointListSize; ++pointNum)
        {
        TubePointType * currentPoint = static_cast< TubePointType * >(
          inputTube->GetPoint( pointNum ) );
        typename TubeType::PointType pointIndex;
        pointIndex = currentPoint->GetPosition();

        typename TubeType::PointType pointWorld;
        pointWorld = inputTube->GetIndexToWorldTransform()->TransformPoint(
          pointIndex );

        typename ImageType::IndexType imageIndex;
        double val = 0;
        if( inputImage->TransformPhysicalPointToIndex( pointWorld,
          imageIndex ) )
          {
          val = inputImage->GetPixel( imageIndex );
          }
        switch( propertyId )
          {
          default:
          case 'R':
            currentPoint->SetRidgeness( val );
            break;
          case 'M':
            currentPoint->SetMedialness( val );
            break;
          case 'B':
            currentPoint->SetBranchness( val );
            break;
          case 'r':
            currentPoint->SetRadius( val );
            break;
          }
        }
      }
    }
}

template< unsigned int DimensionT >
int DoIt( MetaCommand & command )
{
  typename itk::GroupSpatialObject< DimensionT >::Pointer inputTubes;

  MetaCommand::OptionVector parsed = command.GetParsedOptions();

  inputTubes = ReadTubeFile< DimensionT >( command.GetValueAsString(
    "infile" ).c_str() );

  int currentTube = ALL_TUBES_CURRENT;
  MetaCommand::OptionVector::const_iterator it = parsed.begin();
  while( it != parsed.end() )
    {
    if( it->name == "Write" )
      {
      WriteTubeFile< DimensionT >( inputTubes, command.GetValueAsString(
        *it, "filename" ).c_str() );
      }
    else if( it->name == "SelectTube" )
      {
      ::tube::Message( "Select Tube" );
      currentTube = command.GetValueAsInt( *it, "tubeNumber" );
      }
    else if( it->name == "SelectAllTubes" )
      {
      ::tube::Message( "Select All Tubes" );
      currentTube = ALL_TUBES_CURRENT;
      }
    else if( it->name == "LoadRidgeness" )
      {
      ::tube::Message( "Load Ridgness" );
      typedef itk::Image< double, DimensionT > ImageType;
      typename ImageType::Pointer ridgeImage = ReadImageFile< double,
        DimensionT >( command.GetValueAsString( *it, "filename" ).c_str() );
      SetPropertyFromImage< double, DimensionT >( inputTubes, currentTube,
        ridgeImage, 'R' );
      }
    else if( it->name == "LoadMedialness" )
      {
      ::tube::Message( "Load Medialness" );
      typedef itk::Image< double, DimensionT > ImageType;
      typename ImageType::Pointer ridgeImage = ReadImageFile< double,
        DimensionT >( command.GetValueAsString( *it, "filename" ).c_str() );
      SetPropertyFromImage< double, DimensionT >( inputTubes, currentTube,
        ridgeImage, 'M' );
      }
    else if( it->name == "LoadBranchness" )
      {
      ::tube::Message( "Load Branchness" );
      typedef itk::Image< double, DimensionT > ImageType;
      typename ImageType::Pointer ridgeImage = ReadImageFile< double,
        DimensionT >( command.GetValueAsString( *it, "filename" ).c_str() );
      SetPropertyFromImage< double, DimensionT >( inputTubes, currentTube,
        ridgeImage, 'B' );
      }
    else if( it->name == "LoadRadius" )
      {
      ::tube::Message( "Load Radius" );
      typedef itk::Image< double, DimensionT > ImageType;
      typename ImageType::Pointer ridgeImage = ReadImageFile< double,
        DimensionT >( command.GetValueAsString( *it, "filename" ).c_str() );
      SetPropertyFromImage< double, DimensionT >( inputTubes, currentTube,
        ridgeImage, 'r' );
      }
    else if( it->name == "Delete" )
      {
      ::tube::Message( "Delete" );
      }
    else if( it->name == "Reverse" )
      {
      ::tube::Message( "Reverse" );
      }
    else if( it->name == "ComputeTangentsAndNormals" )
      {
      ::tube::Message( "Compute Tangents and Normals" );
      typedef itk::VesselTubeSpatialObject< DimensionT >  TubeType;
      typedef typename TubeType::PointListType            TubePointListType;
      typedef typename TubeType::TubePointType            TubePointType;

      typename TubeType::ChildrenListType::iterator tubeIterator;

      char soTypeName[80];
      strcpy( soTypeName, "VesselTubeSpatialObject" );
      typename TubeType::ChildrenListPointer inputTubeList =
        inputTubes->GetChildren( inputTubes->GetMaximumDepth(),
        soTypeName );
      for( tubeIterator = inputTubeList->begin(); tubeIterator !=
        inputTubeList->end(); tubeIterator++ )
        {
        typename TubeType::Pointer inputTube = ((TubeType *)( tubeIterator->
          GetPointer()));

        if( currentTube == ALL_TUBES_CURRENT || inputTube->GetId() ==
          currentTube )
          {
          ::tube::ComputeTubeTangentsAndNormals< TubeType >( inputTube );
          }
        }
      }
    else if( it->name == "SplitAtPoint" )
      {
      ::tube::Message( "Split at Point" );
      }
    else if( it->name == "MarkAsArtery" )
      {
      ::tube::Message( "Mark as Artery" );
      typedef itk::VesselTubeSpatialObject< DimensionT >  TubeType;
      typedef typename TubeType::PointListType            TubePointListType;
      typedef typename TubeType::TubePointType            TubePointType;

      typename TubeType::ChildrenListType::iterator tubeIterator;

      char soTypeName[80];
      strcpy( soTypeName, "VesselTubeSpatialObject" );
      typename TubeType::ChildrenListPointer inputTubeList =
        inputTubes->GetChildren( inputTubes->GetMaximumDepth(),
        soTypeName );
      for( tubeIterator = inputTubeList->begin(); tubeIterator !=
        inputTubeList->end(); tubeIterator++ )
        {
        typename TubeType::Pointer inputTube = ((TubeType *)( tubeIterator->
          GetPointer()));

        if( currentTube == ALL_TUBES_CURRENT || inputTube->GetId() ==
          currentTube )
          {
          inputTube->SetArtery( true );
          }
        }
      }
    else if( it->name == "MarkAsRoot" )
      {
      ::tube::Message( "Mark as Root" );
      typedef itk::VesselTubeSpatialObject< DimensionT >  TubeType;
      typedef typename TubeType::PointListType            TubePointListType;
      typedef typename TubeType::TubePointType            TubePointType;

      typename TubeType::ChildrenListType::iterator tubeIterator;

      char soTypeName[80];
      strcpy( soTypeName, "VesselTubeSpatialObject" );
      typename TubeType::ChildrenListPointer inputTubeList =
        inputTubes->GetChildren( inputTubes->GetMaximumDepth(),
        soTypeName );
      for( tubeIterator = inputTubeList->begin(); tubeIterator !=
        inputTubeList->end(); tubeIterator++ )
        {
        typename TubeType::Pointer inputTube = ((TubeType *)( tubeIterator->
          GetPointer()));

        if( currentTube == ALL_TUBES_CURRENT || inputTube->GetId() ==
          currentTube )
          {
          inputTube->SetRoot( true );
          }
        }
      }
    else if( it->name == "UniqueIDs" )
      {
      ::tube::Message( "Assign Unique IDs" );
      typedef itk::VesselTubeSpatialObject< DimensionT >  TubeType;
      typedef typename TubeType::PointListType            TubePointListType;
      typedef typename TubeType::TubePointType            TubePointType;

      typename TubeType::ChildrenListType::iterator tubeIterator;

      char soTypeName[80];
      strcpy( soTypeName, "VesselTubeSpatialObject" );
      typename TubeType::ChildrenListPointer inputTubeList =
        inputTubes->GetChildren( inputTubes->GetMaximumDepth(),
        soTypeName );
      unsigned int count = 100;
      for( tubeIterator = inputTubeList->begin(); tubeIterator !=
        inputTubeList->end(); tubeIterator++ )
        {
        typename TubeType::Pointer inputTube = ((TubeType *)( tubeIterator->
          GetPointer()));

        inputTube->SetId( count++ );
        }
      }
    else if( it->name == "MergeWithTube" )
      {
      ::tube::Message( "Merge with Tube" );
      }
    ++it;
    }

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
   //PARSE_ARGS;
  MetaCommand command;

  command.SetName( "ImageMath" );
  command.SetVersion( "1.0" );
  command.SetAuthor( "CADDLab @ UNC" );
  command.SetDescription( "Perform several filters on an set of tubes" );

  command.AddField( "infile", "infile filename",
    MetaCommand::STRING, MetaCommand::DATA_IN );

  command.SetOption( "Write", "w", false,
    "writes current tubes to the designated file" );
  command.AddOptionField( "Write", "filename", MetaCommand::STRING, true,
    "", "output filename", MetaCommand::DATA_OUT );

  command.SetOption( "SelectTube", "t", false,
    "Sets the target tube by ID." );
  command.AddOptionField( "SelectTube", "tubeNumber",
    MetaCommand::INT, true, "", "Tube Number", MetaCommand::DATA_IN );

  command.SetOption( "SelectAllTubes", "T", false,
    "Sets the target tube to be all tubes." );

  command.SetOption( "LoadRidgeness", "R", false,
    "Sets tube points' ridgeness values to the values in the image." );
  command.AddOptionField( "LoadRidgeness", "filename",
    MetaCommand::STRING, true, "", "Image filename", MetaCommand::DATA_IN );

  command.SetOption( "LoadMedialness", "M", false,
    "Sets tube points' medialness values to the values in the image." );
  command.AddOptionField( "LoadMedialness", "filename",
    MetaCommand::STRING, true, "", "Image filename", MetaCommand::DATA_IN );

  command.SetOption( "LoadBranchness", "B", false,
    "Sets tube points' branchness values to the values in the image." );
  command.AddOptionField( "LoadBranchness", "filename",
    MetaCommand::STRING, true, "", "Image filename", MetaCommand::DATA_IN );

  command.SetOption( "LoadRadius", "r", false,
    "Sets tube points' radius values to the values in the image." );
  command.AddOptionField( "LoadRadius", "filename",
    MetaCommand::STRING, true, "", "Image filename", MetaCommand::DATA_IN );

  command.SetOption( "ComputeTangentsAndNormals", "n", false,
    "Compute the current tube's tangents and normals from point sequence."
    );

  command.SetOption( "MarkAsArtery", "a", false,
    "Flag the current tube as an Artery." );

  command.SetOption( "MarkAsRoot", "o", false,
    "Flag the current tube as an Root in a tree." );

  command.SetOption( "UniqueIDs", "u", false,
    "Assign current tube a unique ID." );

  if( !command.Parse( argc, argv ) )
    {
    return EXIT_FAILURE;
    }

  return DoIt< 3 >( command );
}
