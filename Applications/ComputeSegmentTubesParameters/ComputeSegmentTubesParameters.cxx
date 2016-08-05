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
#include <itkTimeProbesCollectorBase.h>

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itktubeComputeSegmentTubesParameters.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ComputeSegmentTubesParametersCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< int VDimension >
int WriteOutputData( std::ofstream & fileStream,
  itk::ContinuousIndex< double, VDimension > & cIndx, double intensity,
  double ridgeness, double roundness, double curvature, double levelness )
{
  for( unsigned int i = 0; i < VDimension; ++i )
    {
    fileStream << cIndx[i] << " ";
    }
  fileStream << intensity << " ";
  fileStream << ridgeness << " ";
  fileStream << roundness << " ";
  fileStream << curvature << " ";
  fileStream << levelness << std::endl;

  return EXIT_SUCCESS;
}

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "ComputeSegmentTubesParameters", CLPProcessInformation );
  progressReporter.Start();
  typedef itk::tube::ComputeSegmentTubesParameters< TPixel, VDimension >
    FilterType;
  FilterType::Pointer filter = FilterType::New();

  typedef typename FilterType::InputImageType       InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;

  typedef typename FilterType::MaskImageType        MaskImageType;
  typedef itk::ImageFileReader< MaskImageType >     MaskReaderType;

  typedef typename FilterType::ScaleImageType       ScaleImageType;
  typedef itk::ImageFileReader< ScaleImageType >    ScaleReaderType;

  timeCollector.Start("Load data");

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName.c_str() );
  try
    {
    reader->Update();
    filter->SetInputImage(  reader->GetOutput() );
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
  maskReader->SetFileName( maskImageFileName.c_str() );
  try
    {
    maskReader->Update();
    filter->SetMaskInputImage( maskReader->GetOutput() );
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading mask: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typename ScaleReaderType::Pointer scaleReader = ScaleReaderType::New();
  scaleReader->SetFileName( scaleImageFileName.c_str() );
  try
    {
    scaleReader->Update();
    filter->SetScaleInputImage( scaleReader->GetOutput() );
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading scale: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  filter->SetMaskBackGroundId( maskBackgroundId );
  filter->SetMaskTubeId( maskTubeId );
  filter->SetParameterFile( outputParametersFile );

  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Compute ridgeness images");
  try
    {
    filter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "ComputeSegmentTubesParameters Update error: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Compute ridgeness images");

  timeCollector.Start("Save Data");

  std::string fileName = outputParametersFile + ".init.txt";
  std::ofstream outputDataStreamInit;
  outputDataStreamInit.open( fileName.c_str(), std::ios::binary |
    std::ios::out );
  outputDataStreamInit.precision( 6 );

  fileName = outputParametersFile + ".tube.txt";
  std::ofstream outputDataStreamTube;
  outputDataStreamTube.open( fileName.c_str(), std::ios::binary |
    std::ios::out );
  outputDataStreamTube.precision( 6 );

  fileName = outputParametersFile + ".bkg.txt";
  std::ofstream outputDataStreamBkg;
  outputDataStreamBkg.open( fileName.c_str(), std::ios::binary |
    std::ios::out );
  outputDataStreamBkg.precision( 6 );

  std::vector< vnl_vector< double > > seedData = filter->GetSeedData();
  std::vector< vnl_vector< double > > tubeData = filter->GetTubeData();
  std::vector< vnl_vector< double > > bkgData = filter->GetBkgData();

  std::vector< itk::ContinuousIndex< double, VDimension > > seedIndex =
    filter->GetSeedDataIndexList();
  std::vector< itk::ContinuousIndex< double, VDimension > > tubeIndex =
    filter->GetTubeDataIndexList();
  std::vector< itk::ContinuousIndex< double, VDimension > > bkgIndex =
    filter->GetBkgDataIndexList();

  for( int i = 0; i < seedData.size(); i++)
    {
    vnl_vector< double > instance = seedData[i];
    WriteOutputData< VDimension >( outputDataStreamInit, seedIndex[i],
            instance[0], instance[1], instance[2], instance[3], instance[4] );
    }
  for( int i = 0; i < tubeData.size(); i++)
    {
    vnl_vector< double > instance = tubeData[i];
    WriteOutputData< VDimension >( outputDataStreamTube, tubeIndex[i],
            instance[0], instance[1], instance[2], instance[3], instance[4] );
    }
  for( int i = 0; i < bkgData.size(); i++)
    {
    vnl_vector< double > instance = bkgData[i];
    WriteOutputData< VDimension >( outputDataStreamBkg, bkgIndex[i],
            instance[0], instance[1], instance[2], instance[3], instance[4] );
    }

  outputDataStreamBkg.close();
  outputDataStreamTube.close();
  outputDataStreamInit.close();

  timeCollector.Stop("Save Data");

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );
}
