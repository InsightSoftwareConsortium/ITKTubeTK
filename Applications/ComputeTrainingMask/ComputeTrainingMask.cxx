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

#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
typename itk::Image< TPixel, VDimension >::Pointer
FindCenterLines( typename itk::Image< TPixel, VDimension >::Pointer input )
{
  typedef itk::Image< TPixel, VDimension >                        ImageType;
  typedef itk::BinaryThinningImageFilter< ImageType, ImageType >  FilterType;

  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput( input );
  filter->Update();
  return filter->GetOutput();
}

template< class TPixel, unsigned int VDimension >
void
ThresholdVolume( typename itk::Image< TPixel, VDimension >::Pointer &input,
                 float threshLow,
                 float threshHigh,
                 float valTrue,
                 float valFalse )
{
  typedef itk::Image< TPixel, VDimension >                        ImageType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdType;

  typename ThresholdType::Pointer threshold = ThresholdType::New();

  threshold->SetInput( input );
  threshold->SetLowerThreshold( threshLow);
  threshold->SetUpperThreshold( threshHigh );
  threshold->SetInsideValue( valTrue );
  threshold->SetOutsideValue( valFalse );
  threshold->Update();
  input = threshold->GetOutput();
  return;
}

template< class TPixel, unsigned int VDimension >
void
ApplyDilateMorphologyFilter( typename itk::Image< TPixel, VDimension >::Pointer &input,
                  float radius,
                  float foregroundValue )
{
  typedef itk::Image< TPixel, VDimension >                          ImageType;
  typedef itk::BinaryBallStructuringElement< TPixel, VDimension >   BallType;
  BallType ball;
  ball.SetRadius( 1 );
  ball.CreateStructuringElement();

  typedef itk::ErodeObjectMorphologyImageFilter
    < ImageType, ImageType, BallType >       ErodeFilterType;
  typedef itk::DilateObjectMorphologyImageFilter
    < ImageType, ImageType, BallType >       DilateFilterType;

  for ( int r = 0; r<radius; r++ )
    {
    typename DilateFilterType::Pointer filter =
      DilateFilterType::New();
    filter->SetKernel( ball );
    filter->SetObjectValue( foregroundValue );
    filter->SetInput( input );
    filter->Update();
    input = filter->GetOutput();
    }
  return;
}

template< class TPixel, unsigned int VDimension >
void
AddVolume( typename itk::Image< TPixel, VDimension >::Pointer &input1,
           typename itk::Image< TPixel, VDimension >::Pointer input2,
           float weight1,
           float weight2 )
{
  typedef itk::Image< TPixel, VDimension >   ImageType;

  itk::ImageRegionIterator< ImageType > it1(input1,
    input1->GetLargestPossibleRegion());
  itk::ImageRegionIterator< ImageType > it2(input2,
    input2->GetLargestPossibleRegion());
  it1.GoToBegin();
  it2.GoToBegin();
  while ( !it1.IsAtEnd() )
    {
    double tf1 = it1.Get();
    double tf2 = it2.Get();
    double tf = weight1*tf1 + weight2*tf2;
    it1.Set( ( TPixel )tf );
    ++it1;
    ++it2;
    }
  return;
}

template< class TPixel, unsigned int VDimension >
void
SaveVolumeAsShort( typename itk::Image< TPixel, VDimension >::Pointer input,
                   const char* fileName )
{
  typedef itk::Image< short, VDimension >                      ImageTypeShort;
  typedef itk::Image< TPixel, VDimension >                     ImageType;
  typedef itk::CastImageFilter< ImageType, ImageTypeShort >    CastFilterType;
  typedef itk::ImageFileWriter< ImageTypeShort >               VolumeWriterType;

  typename CastFilterType::Pointer castFilter =
    CastFilterType::New();
  typename VolumeWriterType::Pointer writer =
    VolumeWriterType::New();

  castFilter->SetInput( input );

  writer->SetFileName( fileName );
  writer->SetInput( castFilter->GetOutput() );
  writer->SetUseCompression( true );
  writer->Update();
  writer->Write();
}

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

  timeCollector.Start( "Find Center Lines" );
  tube::InfoMessage( "Finding Center Lines..." );
  typename ImageType::Pointer centerLines =
    FindCenterLines< TPixel, VDimension >( image );
  timeCollector.Stop( "Find Center Lines" );
  progress = 0.65;
  progressReporter.Report( progress );
  timeCollector.Start( "Threshold and Mathematical Morphology" );
  tube::InfoMessage( "Thresholding..." );
  ThresholdVolume< TPixel, VDimension >( image, 0, gap, 0, 255 );
  ApplyDilateMorphologyFilter< TPixel, VDimension >( image, notVesselWidth, 255 );
  typename ImageType::Pointer dialatedImage = image;
  ApplyDilateMorphologyFilter< TPixel, VDimension >( image, notVesselWidth, 255 );
  tube::InfoMessage( "Creating Not-Vessel Mask..." );
  AddVolume< TPixel, VDimension >( image, dialatedImage, 1, -1 );
  if ( !notVesselMask.empty() )
    {
    SaveVolumeAsShort< TPixel, VDimension >( image, notVesselMask.c_str() );
    }
  tube::InfoMessage( "Creating Vessel Mask..." );
  AddVolume< TPixel, VDimension >( image, centerLines, 0.5, 255 );
  SaveVolumeAsShort< TPixel, VDimension >( image, outputVolume.c_str() );
  timeCollector.Stop( "Threshold and Mathematical Morphology" );
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
