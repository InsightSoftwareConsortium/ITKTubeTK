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

#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeMacro.h"
#include "ComputeTrainingMaskCLP.h"

#include <tubeComputeTrainingMask.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkBinaryThinningImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkCastImageFilter.h>

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  if ( VDimension != 2 && VDimension != 3 )
    {
    tube::ErrorMessage(
      "Error: Only 2D and 3D data is currently supported.");
    return EXIT_FAILURE;
    }

  typedef itk::Image< TPixel, VDimension >            ImageType;
  typedef itk::ImageFileReader< ImageType >           ImageReaderType;

  // The timeCollector to perform basic profiling of algorithmic components
  itk::TimeProbesCollectorBase timeCollector;

  timeCollector.Start( "Loading Input Volume Mask File" );
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "ComputeTrainingMask",
    CLPProcessInformation );
  progressReporter.Start();
  float progress = 0;

  typename ImageReaderType::Pointer imReader;
  imReader = ImageReaderType::New();
  typename ImageType::Pointer image;
  tube::InfoMessage( "Reading volume mask..." );
  imReader->SetFileName(inputVolume.c_str());
  try
    {
    imReader->Update();
    image = imReader->GetOutput();
    }
  catch ( itk::ExceptionObject & err )
    {
    tube::FmtErrorMessage( "Cannot read volume mask file: %s",
      err.what() );
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Loading Input Volume Mask File" );
  progress = 0.35;
  progressReporter.Report( progress );
  timeCollector.Start( "Compute training mask" );
  tube::InfoMessage( "Compute training mask..." );
  typedef tube::ComputeTrainingMask<ImageType> ComputeTrainingMaskType;
  typename ComputeTrainingMaskType::Pointer filter =
    ComputeTrainingMaskType::New();
  filter->SetInput(imReader->GetOutput());
  filter->SetGap(gap);
  filter->SetNotVesselWidth(notVesselWidth);
  filter->Update();
  progress = 0.65;
  progressReporter.Report( progress );
  timeCollector.Stop( "Compute training mask" );
  typedef typename ComputeTrainingMaskType::ImageTypeShort ImageTypeShort;
  typedef itk::ImageFileWriter<ImageTypeShort>             VolumeWriterType;

  typename VolumeWriterType::Pointer writer = VolumeWriterType::New();
  if ( !notVesselMask.empty() )
    {
    timeCollector.Start( "Creating not-Vessel Mask" );
    tube::InfoMessage( "Creating not-Vessel Mask..." );
    writer->SetFileName(notVesselMask);
    writer->SetInput(filter->GetNotVesselMask());
    writer->Update();
    timeCollector.Stop( "Creating not-Vessel Mask" );
    }
  timeCollector.Start( "Creating Vessel Mask" );
  tube::InfoMessage( "Creating Vessel Mask..." );
  writer->SetFileName(outputVolume.c_str());
  writer->SetInput( filter->GetOutput() );
  writer->Update();
  timeCollector.Stop("Creating Vessel Mask" );
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

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
