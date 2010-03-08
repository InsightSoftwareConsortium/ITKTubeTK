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

#include "itkNormalizeImageFilter.h"
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
  tube::CLIProgressReporter    progressReporter( "MatchROIs",
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
  
      ti = inputSize[i] - (int)( lowerCropSize[i] + roiSize[i] );
      if( ti < 0 )
        {
        upperCropSize[i] = 0;
        }
      else if( ti >= (int)(inputSize[i]) - (int)(lowerCropSize[i]) )
        {
        ti = (int)(inputSize[i]) - (int)(lowerCropSize[i]);
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

  typename ImageType::Pointer orgMask = curMask;

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

  if( !disableRegisterROIs )
    {
    timeCollector.Start("RegisterROIs");
    typedef itk::NormalizeImageFilter< ImageType, ImageType > NormFilterType;

    typename NormFilterType::Pointer norm1 = NormFilterType::New();
    norm1->SetInput( curVolume );
    norm1->Update();

    typename NormFilterType::Pointer norm2 = NormFilterType::New();
    norm2->SetInput( curMask );
    norm2->Update();

    typedef itk::RigidImageToImageRegistrationMethod< ImageType >
                                                    RegistrationMethodType;
    typename RegistrationMethodType::Pointer reg = RegistrationMethodType::New();
    reg->SetFixedImage( norm1->GetOutput() );
    reg->SetMovingImage( norm2->GetOutput() );
    reg->SetMaxIterations( 100 );
    typename RegistrationMethodType::TransformParametersScalesType scales;
    scales = reg->GetTransformParametersScales();
    scales[0] = 1.0/0.02;
    scales[1] = 1.0/5.0;
    scales[2] = 1.0/5.0;
    reg->SetTransformParametersScales( scales );
    reg->SetUseEvolutionaryOptimization( true );
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

    typename ResamplerType::Pointer orgResampler = ResamplerType::New();
    typename InterpolatorType::Pointer orgInterpolator = InterpolatorType::New();
    orgInterpolator->SetInputImage( orgMask );
    orgResampler->SetInput( orgMask );
    orgResampler->SetInterpolator( orgInterpolator.GetPointer() );
    orgResampler->SetOutputParametersFromImage( orgMask );
    orgResampler->SetTransform( reg->GetTypedTransform() );
    orgResampler->Update();
    orgMask = orgResampler->GetOutput();

    timeCollector.Stop("RegisterROIs");
    }

  if( outputSize.size() > 0 )
    {
    timeCollector.Start("Crop2");
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      int ti = (roiSize[i])/2 - (outputSize[i]-1)/2;
      if( ti < 0 )
        {
        lowerCropSize[i] = 0;
        }
      else if( ti >= (int)(roiSize[i]) )
        {
        lowerCropSize[i] = roiSize[i]-1;
        }
      else
        {
        lowerCropSize[i] = ti;
        }
  
      ti = roiSize[i] - (int)( lowerCropSize[i] + outputSize[i] );
      if( ti < 0 )
        {
        upperCropSize[i] = 0;
        }
      else if( ti >= (int)(roiSize[i]) - (int)(lowerCropSize[i]) )
        {
        ti = (int)(roiSize[i]) - (int)(lowerCropSize[i]);
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
    timeCollector.Stop("Crop2");
    }
  
  if( !disableNormalize )
    {
    timeCollector.Start("Normalize");
  
    typedef itk::ImageRegionIterator< ImageType > IterType;
    IterType volIter( curVolume, curVolume->GetLargestPossibleRegion() );
    IterType maskIter( curMask, curMask->GetLargestPossibleRegion() );
    IterType orgMaskIter( orgMask, orgMask->GetLargestPossibleRegion() );
  
    int numBins = 50;
    float volFgHisto[50];
    float volBgHisto[50];
    float maskFgHisto[50];
    float maskBgHisto[50];
    for( int i=0; i<numBins; i++ )
      {
      volFgHisto[i] = 0;
      volBgHisto[i] = 0;
      maskFgHisto[i] = 0;
      maskBgHisto[i] = 0;
      }
    volIter.GoToBegin();
    maskIter.GoToBegin();
    orgMaskIter.GoToBegin();
    double volMax = volIter.Get();
    double maskMax = maskIter.Get();
    double volMin = volMax;
    double maskMin = maskMax;
    while( !volIter.IsAtEnd() )
      {
      if( volIter.Get() > volMax )
        {
        volMax = volIter.Get();
        }
      if( volIter.Get() < volMin )
        {
        volMin = volIter.Get();
        }
      if( maskIter.Get() > maskMax )
        {
        maskMax = maskIter.Get();
        }
      if( maskIter.Get() < maskMin )
        {
        maskMin = maskIter.Get();
        }
  
      ++volIter;
      ++maskIter;
      }
  
    int bin;
    volIter.GoToBegin();
    maskIter.GoToBegin();
    orgMaskIter.GoToBegin();
    while( !volIter.IsAtEnd() )
      {
      if( orgMaskIter.Get() == foreground )
        {
        bin = (int)( (((double)(volIter.Get())-volMin)/(volMax-volMin)) * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++volFgHisto[ bin ];
        if( bin > 0 )
          {
          volFgHisto[ bin-1 ] += 0.5;
          }
        if( bin < numBins - 1 )
          {
          volFgHisto[ bin+1 ] += 0.5;
          }
  
        bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin)) * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++maskFgHisto[ bin ];
        if( bin > 0 )
          {
          maskFgHisto[ bin-1 ] += 0.5;
          }
        if( bin < numBins - 1 )
          {
          maskFgHisto[ bin+1 ] += 0.5;
          }
        }
      else
        {
        bin = (int)( (((double)(volIter.Get())-volMin)/(volMax-volMin)) * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++volBgHisto[ bin ];
        if( bin > 0 )
          {
          volBgHisto[ bin-1 ] += 0.5;
          }
        if( bin < numBins - 1 )
          {
          volBgHisto[ bin+1 ] += 0.5;
          }
  
        bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin)) * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++maskBgHisto[ bin ];
        if( bin > 0 )
          {
          maskBgHisto[ bin-1 ] += 0.5;
          }
        if( bin < numBins - 1 )
          {
          maskBgHisto[ bin+1 ] += 0.5;
          }
        }
  
      ++volIter;
      ++maskIter;
      ++orgMaskIter;
      }
    /*
    for( bin=0; bin<numBins; bin++ )
      {
      std::cout << volBgHisto[bin] << " : " << volFgHisto[bin] << " - "
                << maskBgHisto[bin] << " : " << maskFgHisto[bin]
                << std::endl;
      }
      */
  
    float maxVFg = 0;
    int maxVFgI = -1;
    float maxVBg = 0;
    int maxVBgI = -1;
    float maxMFg = 0;
    int maxMFgI = -1;
    float maxMBg = 0;
    int maxMBgI = -1;
    for( int i=0; i<numBins; i++ )
      {
      if( volFgHisto[i] > maxVFg )
        {
        maxVFg = volFgHisto[i];
        maxVFgI = i;
        }
      if( volBgHisto[i] > maxVBg )
        {
        maxVBg = volBgHisto[i];
        maxVBgI = i;
        }
      if( maskFgHisto[i] > maxMFg )
        {
        maxMFg = maskFgHisto[i];
        maxMFgI = i;
        }
      if( maskBgHisto[i] > maxMBg )
        {
        maxMBg = maskBgHisto[i];
        maxMBgI = i;
        }
      }
  
    maxVFg = (double)(maxVFgI)/(numBins-1) * (volMax - volMin) + volMin;
    maxVBg = (double)(maxVBgI)/(numBins-1) * (volMax - volMin) + volMin;
    maxMFg = (double)(maxMFgI)/(numBins-1) * (maskMax - maskMin) + maskMin;
    maxMBg = (double)(maxMBgI)/(numBins-1) * (maskMax - maskMin) + maskMin;
    /*
    std::cout << "VFg = " << maxVFg << std::endl;
    std::cout << "VBg = " << maxVBg << std::endl;
    std::cout << "MFg = " << maxMFg << std::endl;
    std::cout << "MBg = " << maxMBg << std::endl;
    */
   
    if( maxVFgI != -1 &&
        maxVBgI != -1 &&
        maxMFgI != -1 &&
        maxMBgI != -1 )
      {
      typedef itk::ShiftScaleImageFilter< ImageType, ImageType > FilterType;
      typename FilterType::Pointer filter = FilterType::New();
   
      double scale = (maxVFg - maxVBg) / (maxMFg - maxMBg);
      double shift = (maxVBg - (maxMBg*scale));
      filter = FilterType::New();
      filter->SetInput( curMask );
      filter->SetShift( shift );
      filter->SetScale( scale );
   
      filter->Update();
      curMask = filter->GetOutput();
      }
  
    timeCollector.Stop("Normalize");
    }

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
