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

#include "tubeSegmentTubes.h"

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "SegmentTubesCLP.h"

#include <sstream>

#define PARSE_ARGS_FLOAT_ONLY

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and
//   forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "SegmentTubes",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                     PixelType;
  typedef itk::Image< PixelType, VDimension >        ImageType;
  typedef itk::ImageFileReader< ImageType >          ReaderType;

  typedef tube::SegmentTubes< ImageType >            SegmentTubesFilterType;

  typedef typename SegmentTubesFilterType::TubeMaskImageType
    MaskImageType;

  typedef itk::ImageFileReader< MaskImageType >      MaskReaderType;
  typedef itk::ImageFileWriter< MaskImageType >      MaskWriterType;

  typedef typename ImageType::PointType              PointType;

  typedef typename SegmentTubesFilterType::ContinuousIndexType
    IndexType;

  typedef std::vector< IndexType >                   IndexListType;
  typedef std::vector< double >                      RadiusListType;
  typedef std::vector< PointType >                   PointListType;

  typedef itk::SpatialObjectWriter< VDimension >     TubesWriterType;
  typedef itk::SpatialObjectReader< VDimension >     TubesReaderType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();
  typename SegmentTubesFilterType::Pointer segmentTubesFilter
    = SegmentTubesFilterType::New();
  segmentTubesFilter->SetInputImage( inputImage );

  if( !radiusInputVolume.empty() )
    {
    typename ReaderType::Pointer radiusReader = ReaderType::New();
    radiusReader->SetFileName( radiusInputVolume.c_str() );
    try
      {
      radiusReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading radius volume: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    segmentTubesFilter->SetRadiusInputImage( radiusReader->GetOutput() );
    }
  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  segmentTubesFilter->SetRadiusInObjectSpace( radiusInObjectSpace );

  IndexType seedIndex;
  IndexListType seedIndexList;
  if( !seedI.empty() )
    {
    seedIndexList.clear();
    for( unsigned int seedINum = 0; seedINum < seedI.size(); ++seedINum )
      {
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        seedIndex[i] = seedI[seedINum][i];
        }
      seedIndexList.push_back( seedIndex );
      }
    segmentTubesFilter->SetSeedsInIndexSpaceList( seedIndexList );
    }

  PointType point;
  PointListType pointList;
  RadiusListType radiusList;
  if( !seedP.empty() )
    {
    pointList.clear();
    radiusList.clear();
    for( size_t seedNum = 0; seedNum < seedP.size(); ++seedNum )
      {
      for( unsigned int i = 0; i < VDimension; i++ )
        {
        point[i] = seedP[seedNum][i];
        }

      pointList.push_back( point );
      radiusList.push_back( radiusInObjectSpace );
      }
    segmentTubesFilter->SetSeedsInObjectSpaceList( pointList );
    segmentTubesFilter->SetSeedRadiiInObjectSpaceList( radiusList );
    }

  if( !seedsInIndexSpaceListFile.empty() )
    {
    seedIndexList.clear();
    std::ifstream readStream;
    readStream.open( seedsInIndexSpaceListFile.c_str(), std::ios::binary |
      std::ios::in );
    std::string line;
    double seedRadius;
    while( std::getline( readStream, line ) )
      {
      std::istringstream iss( line );
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        iss >> seedIndex[i];
        }
      iss >> seedRadius;
      seedIndexList.push_back( seedIndex );
      radiusList.push_back( seedRadius );
      }
    segmentTubesFilter->SetSeedsInIndexSpaceList( seedIndexList );
    segmentTubesFilter->SetSeedRadiiInObjectSpaceList( radiusList );
    }

  if( !seedsInObjectSpaceListFile.empty() )
    {
    double seedRadius;
    pointList.clear();
    radiusList.clear();
    std::ifstream readStream;
    readStream.open( seedsInObjectSpaceListFile.c_str(), std::ios::binary |
      std::ios::in );
    std::string line;
    while( std::getline( readStream, line ) )
      {
      std::istringstream iss( line );
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        iss >> point[i];
        }
      iss >> seedRadius;
      pointList.push_back( point );
      radiusList.push_back( seedRadius );
      }
    segmentTubesFilter->SetSeedsInObjectSpaceList( pointList );
    segmentTubesFilter->SetSeedRadiiInObjectSpaceList( radiusList );
    }

  if( !seedMask.empty() )
    {
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( seedMask.c_str() );
    try
      {
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading mask: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    segmentTubesFilter->SetSeedMask( maskReader->GetOutput() );

    if( !radiusMask.empty() )
      {
      typename ReaderType::Pointer radiusReader = ReaderType::New();
      radiusReader->SetFileName( radiusMask.c_str() );
      try
        {
        radiusReader->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        tube::ErrorMessage( "Reading scale: Exception caught: "
                            + std::string( err.GetDescription() ) );
        timeCollector.Report();
        return EXIT_FAILURE;
        }
      segmentTubesFilter->SetSeedRadiusMask( radiusReader->GetOutput() );
      }
    segmentTubesFilter->SetSeedMaskStride( seedMaskStride );
    }

  if( !existingVesselsMask.empty() )
    {
    typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( existingVesselsMask.c_str() );
    try
      {
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading vessels mask: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    segmentTubesFilter->SetTubeMaskImage( maskReader->GetOutput() );
    }

  unsigned int numberOfPriorChildren = 0;
  if( !existingVessels.empty() )
    {
    typename TubesReaderType::Pointer tubeReader = TubesReaderType::New();

    try
      {
      tubeReader->SetFileName( existingVessels.c_str() );
      tubeReader->Update();
      }
    catch( ... )
      {
      std::cerr << "Error:: No readable Tubes found " << std::endl;
      }
    segmentTubesFilter->SetTubeGroup( tubeReader->GetGroup() );
    numberOfPriorChildren = tubeReader->GetGroup()->GetNumberOfChildren();
    }

  if( !parametersFile.empty() )
    {
    segmentTubesFilter->LoadParameterFile( parametersFile.c_str() );
    }

  segmentTubesFilter->SetBorderInIndexSpace( border );

  timeCollector.Start( "Ridge Extractor" );

  segmentTubesFilter->ProcessSeeds();

  timeCollector.Stop( "Ridge Extractor" );

  if( segmentTubesFilter->GetTubeGroup()->GetNumberOfChildren() ==
    numberOfPriorChildren )
    {
    std::cerr << "Failed to extract any tubes." << std::endl;
    progressReporter.Report( 1.0 );
    progressReporter.End();
 
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Start( "Save tubes and mask" );
  // Save Tubes
  typename TubesWriterType::Pointer soWriter = TubesWriterType::New();
  soWriter->SetFileName( outputTubeFile );
  soWriter->SetInput( segmentTubesFilter->GetTubeGroup() );
  soWriter->Update();

  // Save Tube Mask Image
  if( !outputTubeImage.empty() )
    {
    typename MaskWriterType::Pointer writer = MaskWriterType::New();
    writer->SetFileName( outputTubeImage );
    writer->SetInput( segmentTubesFilter->GetTubeMaskImage() );
    writer->SetUseCompression( true );
    writer->Update();
    }
  timeCollector.Stop( "Save tubes and mask" );

  progressReporter.Report( 1.0 );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  try
    {
    PARSE_ARGS;
    }
  catch( const std::exception & err )
    {
    tube::ErrorMessage( err.what() );
    return EXIT_FAILURE;
    }
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
