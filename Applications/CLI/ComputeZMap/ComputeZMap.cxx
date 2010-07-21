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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkJointHistogramImageFunction.h"
#include "itkImageRandomConstIteratorWithIndex.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "ComputeZMapCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "ComputeZMap",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                                 PixelType;
  typedef itk::OrientedImage<PixelType, dimensionT>             ImageType;
  typedef itk::ImageFileReader<ImageType>                       ReaderType;
  typedef itk::ImageFileWriter<ImageType>                       WriterType;

  typedef itk::JointHistogramImageFunction<ImageType>           HistCalcType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>     IterType;
  typedef itk::ImageRandomConstIteratorWithIndex<ImageType>     RandIterType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume1.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer image1 = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( inputVolume2.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer image2 = reader->GetOutput();

  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename HistCalcType::Pointer calc = HistCalcType::New();
  calc->SetInputImage( image1 );
  calc->SetInputMask( image2 );

  RandIterType itr( image1, image1->GetLargestPossibleRegion() );
  itr.SetNumberOfSamples( samples );
  itr.GoToBegin();
  while( !itr.IsAtEnd() )
    {
    typename HistCalcType::PointType curPoint;
    image1->TransformIndexToPhysicalPoint( itr.GetIndex(), curPoint );
    calc->Precompute( curPoint );
    ++itr;
    }

  typename ImageType::Pointer curImage = ImageType::New();
  typename ImageType::RegionType region = image1->GetLargestPossibleRegion();  
  curImage->SetRegions( region );
  curImage->Allocate();
  curImage->FillBuffer( 0 );

  itr.GoToBegin();
  PixelType sum = 0;
  while( !itr.IsAtEnd() )
    {
    typename HistCalcType::PointType curPoint;
    image1->TransformIndexToPhysicalPoint( itr.GetIndex(), curPoint );
    PixelType cur = calc->Evaluate( curPoint );
    sum += cur;
    curImage->SetPixel( itr.GetIndex(), cur );
    ++itr;
    }

  std::cout << sum / static_cast<float>( samples ) << std::endl;


  timeCollector.Start("Save data");
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( curImage );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End( );
  
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume1, argc, argv );
}
