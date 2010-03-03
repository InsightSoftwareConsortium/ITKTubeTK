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

#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Application-specific includes
#include "itkCropImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"

#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"

#include "itkRecursiveGaussianImageFilter.h"

#include "itkRigidImageToImageRegistrationMethod.h"
#include "itkResampleImageFilter.h"

#include "itkShiftScaleImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "MatchROIsCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "SampleCLIApplication",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef float                                                   PixelType;
  typedef itk::Image< PixelType,  dimensionT >                    ImageType;
  
  /** Read input images */
  typename ImageType::Pointer curVolume;
  typename ImageType::Pointer curMask;

  timeCollector.Start("Read");
  {
  typedef itk::ImageFileReader< ImageType >                       ReaderType;

  typename ReaderType::Pointer readerVolume = ReaderType::New();
  typename ReaderType::Pointer readerMask = ReaderType::New();

  //read input image  
  readerVolume->SetFileName( inputVolume.c_str() );
  readerMask->SetFileName( inputMask.c_str() );

  try
    {
    readerVolume->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume. Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  try
    {
    readerMask->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading mask. Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  curVolume = readerVolume->GetOutput();
  curMask = readerMask->GetOutput();
  }
  timeCollector.Stop("Read");

  typename ImageType::SizeType inputSize = curVolume->GetLargestPossibleRegion().GetSize();
  typename ImageType::SizeType lowerCropSize;
  typename ImageType::SizeType upperCropSize;

  /** Crop input images to ROI */
  timeCollector.Start("Crop");
  {
  for( unsigned int i=0; i<dimensionT; i++ )
    {
    int ti = roiCenter[i] - (roiSize[i]-1)/2;
    if( ti < 0 )
      {
      lowerCropSize[i] = 0;
      }
    else if( ti >= (int)(inputSize[i]) )
      {
      lowerCropSize[i] = inputSize[i]-1;
      }
    else
      {
      lowerCropSize[i] = ti;
      }

    ti = inputSize[i] - ( lowerCropSize[i] + roiSize[i] );
    if( ti < 0 )
      {
      upperCropSize[i] = 0;
      }
    else if( ti >= (int)(inputSize[i]) )
      {
      ti = inputSize[i]-1;
      }
    upperCropSize[i] = ti;
    }

  typedef itk::CropImageFilter< ImageType, ImageType > CropFilterType;
  typename CropFilterType::Pointer cropVolumeFilter = CropFilterType::New();
  typename CropFilterType::Pointer cropMaskFilter = CropFilterType::New();

  cropVolumeFilter->SetLowerBoundaryCropSize( lowerCropSize );
  cropVolumeFilter->SetUpperBoundaryCropSize( upperCropSize );
  cropVolumeFilter->SetInput( curVolume );
  cropVolumeFilter->Update();
  curVolume = cropVolumeFilter->GetOutput();

  cropMaskFilter->SetLowerBoundaryCropSize( lowerCropSize );
  cropMaskFilter->SetUpperBoundaryCropSize( upperCropSize );
  cropMaskFilter->SetInput( curMask );
  cropMaskFilter->Update();
  curMask = cropMaskFilter->GetOutput();
  }
  timeCollector.Stop("Crop");

  if( foreground != 1 || background != 0 )
    {
    timeCollector.Start("Fg/Bg");

    typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( curMask );
    if( foreground != 1 )
      {
      filter->SetLowerThreshold( foreground );
      filter->SetUpperThreshold( foreground );
      filter->SetInsideValue( 1 );
      filter->SetOutsideValue( 0 );
      }
    else
      {
      filter->SetLowerThreshold( background );
      filter->SetUpperThreshold( background );
      filter->SetInsideValue( 0 );
      filter->SetOutsideValue( 1 );
      }
    filter->Update();
    curMask = filter->GetOutput();

    timeCollector.Stop("Fg/Bg");
    }
  
  if( erode > 0 )
    {
    timeCollector.Start("Erode");

    typedef itk::BinaryBallStructuringElement<PixelType, dimensionT >  BallType;
    BallType ball;
    ball.SetRadius( 1 );
    ball.CreateStructuringElement();

    typedef itk::BinaryErodeImageFilter
                 <ImageType, ImageType, BallType>       ErodeFilterType;

    for(int r=0; r<erode; r++)
      {
      typename ErodeFilterType::Pointer filter = ErodeFilterType::New();
      filter->SetBackgroundValue( 0 );
      filter->SetErodeValue( 1 );
      filter->SetKernel( ball );
      filter->SetInput( curMask );
      filter->Update();
      curMask = filter->GetOutput();
      }

    timeCollector.Stop("Erode");
    }

  if( gaussianBlur > 0 )
    {
    timeCollector.Start("Blur");

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    for(unsigned int i=0; i<dimensionT; i++)
      {
      filter = FilterType::New();
      filter->SetInput( curMask );
      // filter->SetNormalizeAcrossScale( true );
      filter->SetSigma( gaussianBlur );

      filter->SetOrder( 
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );

      filter->Update();
      curMask = filter->GetOutput();
      }

    timeCollector.Stop("Blur");
    }

  if( registerROIs )
    {
    timeCollector.Start("RegisterROIs");
    typedef itk::RigidImageToImageRegistrationMethod< ImageType >
                                                    RegistrationMethodType;
    typename RegistrationMethodType::Pointer reg = RegistrationMethodType::New();
    reg->SetFixedImage( curVolume );
    reg->SetMovingImage( curMask );
    reg->SetMaxIterations( 100 );
    reg->SetUseEvolutionaryOptimization( false );
    int numSamples = 1;
    for( unsigned int i=0; i<dimensionT; i++)
      {
      numSamples *= inputSize[i];
      }
    reg->SetNumberOfSamples( numSamples * 0.2 );
    reg->Update();

    typedef itk::ResampleImageFilter< ImageType, ImageType, double > ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();

    typedef itk::LinearInterpolateImageFunction< ImageType, double > InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage( curMask );

    resampler->SetInput( curMask );
    resampler->SetInterpolator( interpolator.GetPointer() );
    resampler->SetOutputParametersFromImage( curMask );
    resampler->SetTransform( reg->GetTypedTransform() );
    resampler->Update();
    curMask = resampler->GetOutput();

    timeCollector.Stop("RegisterROIs");
    }

  /*
  if( outOffset != 0 || outScale != 1 )
    {
    typedef itk::ShiftScaleImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter = FilterType::New();
    filter->SetInput( curMask );
    filter->SetShift( outOffset );
    filter->SetScale( outScale );

    filter->Update();
    curMask = filter->GetOutput();
    }
    */

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;
  typename ImageWriterType::Pointer writerVolume = ImageWriterType::New();
  typename ImageWriterType::Pointer writerMask = ImageWriterType::New();

  writerVolume->SetFileName( outputVolumeROI.c_str() );
  writerVolume->SetInput ( curVolume );
  writerMask->SetFileName( outputMaskROI.c_str() );
  writerMask->SetInput ( curMask );

  try
    {
    writerVolume->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume. Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  
  try
    {
    writerMask->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing mask. Exception caught: " 
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  
  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char **argv )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
