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

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "metaScene.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkImageFileReader.h"
#include "itkGroupSpatialObject.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"

#include "itkSpeedFunctionToPathFilter.h"
#include "itkSpeedFunctionPathInformation.h"
#include "itkIterateNeighborhoodOptimizer.h"
#include "itkSingleImageCostFunction.h"
#include "itkPathIterator.h"

#include "itktubeRadiusExtractor2.h"

#include "tubeSegmentTubesUsingMinimalPath.h"
#include <sstream>

#include "SegmentTubeUsingMinimalPathCLP.h"

template< class TPixel, unsigned int DimensionT >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int DimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;
  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "Extract Minimal Path",
    CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                                PixelType;
  typedef itk::Image< PixelType, DimensionT >                   ImageType;
  typedef itk::ImageFileReader< ImageType >                     ReaderType;
  typedef itk::SpatialObjectReader< DimensionT >                TubesReaderType;
  typedef itk::GroupSpatialObject< DimensionT >                 TubeGroupType;
  typedef itk::Point< double, DimensionT >                      PointType;

  timeCollector.Start( "Load data" );

  typedef tube::SegmentTubesUsingMinimalPath< DimensionT, PixelType >
                                                         FilterType;
  typename FilterType::Pointer filter = FilterType::New();

  //Read input Image
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( InputImage.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    tube::ErrorMessage( out.str() );
    timeCollector.Stop( "Load data" );
    return EXIT_FAILURE;
    }

  typename ImageType::Pointer speed = reader->GetOutput();
  speed->DisconnectPipeline();
  filter->SetSpeedImage( speed );

  //Read radius extraction Image
  if( !RadiusImage.empty() )
    {
    reader->SetFileName( RadiusImage.c_str() );
    try
      {
      reader->Update();
      filter->SetRadiusImage( reader->GetOutput() );
      filter->SetExtractRadius( true );
      filter->SetStartRadius( StartRadius );
      filter->SetMaxRadius( MaxRadius );
      filter->SetStepSizeForRadiusEstimation( StepRadius );
      }
    catch( itk::ExceptionObject & err )
      {
      std::stringstream out;
      out << "ExceptionObject caught !" << std::endl;
      out << err << std::endl;
      tube::ErrorMessage( out.str() );
      timeCollector.Stop( "Load data" );
      return EXIT_FAILURE;
      }
    }

  //Read target tube
  typename TubesReaderType::Pointer tubeFileReader = TubesReaderType::New();
  if( !TargetTubeFileName.empty() )
    {
    try
      {
      tubeFileReader->SetFileName( TargetTubeFileName.c_str() );
      tubeFileReader->Update();
      filter->SetInput( tubeFileReader->GetGroup() );
      filter->SetExtractEndPointFromTargetTube( true );
      if( ConnectionOption == "Connect_To_Target_Tube_Surface" )
        {
        filter->SetConnectToTargetTubeSurface( true );
        }
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Error loading TRE File: "
        + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }
  else
    {
    typename TubeGroupType::Pointer inputTubeGroup =TubeGroupType::New();
    filter->SetInput( inputTubeGroup );
    }

  timeCollector.Stop( "Load data" );
  progressReporter.Report( 0.1 );

  timeCollector.Start( "Set parameters" );

  PointType startPathPoint;
  if( StartPoint.size() == 1 )
    {
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      startPathPoint[i]=StartPoint[0][i];
      }
   filter->SetStartPoint( startPathPoint );
    }
  else
    {
    tubeErrorMacro(
      << "Error: Path must contain at only one Start Point" );
    timeCollector.Stop( "Set parameters" );
    return EXIT_FAILURE;
    }

  if( IntermediatePoints.size() >= 1 )
  {
  std::vector< PointType > intermediatePathPoints;
  for( unsigned int k = 0; k < IntermediatePoints.size(); k++ )
    {
    PointType pathPoint;
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
       pathPoint[i]=IntermediatePoints[k][i];
      }
    intermediatePathPoints.push_back( pathPoint );
    }
  filter->SetIntermediatePoints( intermediatePathPoints );
  }

  if( EndPoint.size() == 1 )
    {
    PointType targetPathPoint;
    for( unsigned int i = 0; i < DimensionT; i++ )
      {
      targetPathPoint[i]=EndPoint[0][i];
      }
    filter->SetEndPoint( targetPathPoint );
    }
  else if ( TargetTubeFileName.empty() && EndPoint.size() == 0 )
    {
    tubeErrorMacro(
      << "Error: Atmost one End/Target Point or Target Tube should be provided. " );
    timeCollector.Stop( "Set parameters" );
    return EXIT_FAILURE;
    }

 filter->SetOptimizationMethod( Optimizer );
 filter->SetOptimizerTerminationValue( TerminationValue );
 filter->SetOptimizerNumberOfIterations( NumberOfIterations );
 filter->SetOptimizerStepLengthFactor( StepLengthFactor );
 filter->SetOptimizerStepLengthRelax( StepLengthRelax );

  timeCollector.Stop( "Set parameters" );
  progressReporter.Report( 0.2 );

  timeCollector.Start( "Extract minimal path" );

  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::stringstream out;
    out << "ExceptionObject caught !" << std::endl;
    out << err << std::endl;
    tube::ErrorMessage( out.str() );
    timeCollector.Stop( "Extract minimal path" );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Extract minimal path" );
  progressReporter.Report( 0.4 );

  timeCollector.Start( "Write output data" );

  // Write output TRE file
  typedef itk::SpatialObjectWriter< DimensionT > TubeWriterType;
  typename TubeWriterType::Pointer tubeWriter = TubeWriterType::New();
  try
    {
    tubeWriter->SetFileName( OutputTREFile.c_str() );
    tubeWriter->SetInput( filter->GetOutput() );
    tubeWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Error writing TRE file: "
      + std::string( err.GetDescription() ) );
    timeCollector.Stop( "Write output data" );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Write output data" );
  progressReporter.Report( 1.0 );

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;
  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImage.
  return tube::ParseArgsAndCallDoIt( InputImage, argc, argv );
}
