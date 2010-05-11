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
#ifndef __tubeCompareCroppedROIs_txx
#define __tubeCompareCroppedROIs_txx

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

#include "itkLinearInterpolateImageFunction.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizeImageFilter.h"

namespace tube
{

template< class pixelT, unsigned int dimensionT >
CompareCroppedROIs< pixelT, dimensionT>::
CompareCroppedROIs( void )
{
  m_VolImage = NULL;
  m_MaskImage = NULL;
  m_OrgMaskImage = NULL;
  m_Foreground = 0;
  m_Erode = 4;
  m_Dilate = 2;
  m_GaussianBlur = 10;
  m_UseRegistration = true;
  m_UseRegistrationTransform = false;
  m_RegistrationTransform = RegistrationMethodType::TransformType::New();
  m_RegistrationTransform->SetIdentity();
  m_Normalize = true;
  m_OutputSize.resize(0);
  m_UseMeanSquaresMetric = true;
  m_UseCorrelationMetric = false;
  m_SamplingRate = 0.2;

  m_TimeCollector = NULL;
  m_ProgressReporter = NULL;
  m_ProgressStart = 0;
  m_ProgressRange = 1.0;

  m_GoF = 0;
}

template< class pixelT, unsigned int dimensionT >
CompareCroppedROIs< pixelT, dimensionT>::
~CompareCroppedROIs( void )
{
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetVolumeImage( typename ImageType::Pointer volImage )
{
  m_VolImage = volImage;
}

template< class pixelT, unsigned int dimensionT >
typename itk::OrientedImage< float, dimensionT >::Pointer 
CompareCroppedROIs< pixelT, dimensionT>::
GetVolumeImage( void )
{
  return m_VolImage;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetMaskImage( typename ImageType::Pointer maskImage )
{
  m_MaskImage = maskImage;
}

template< class pixelT, unsigned int dimensionT >
typename itk::OrientedImage< float, dimensionT >::Pointer 
CompareCroppedROIs< pixelT, dimensionT>::
GetMaskImage( void )
{
  return m_MaskImage;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetOriginalMaskImage( typename ImageType::Pointer orgMaskImage )
{
  m_OrgMaskImage = orgMaskImage;
}

template< class pixelT, unsigned int dimensionT >
typename itk::OrientedImage< float, dimensionT >::Pointer 
CompareCroppedROIs< pixelT, dimensionT>::
GetOriginalMaskImage( void )
{
  return m_OrgMaskImage;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetForeground( float foreground )
{
  m_Foreground = foreground;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetErode( int erode )
{
  m_Erode = erode;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetDilate( int dilate )
{
  m_Dilate = dilate;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetGaussianBlur( float gaussianBlur )
{
  m_GaussianBlur = gaussianBlur;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetUseRegistration( bool reg )
{
  m_UseRegistration = reg;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetUseRegistrationTransform( bool reg )
{
  m_UseRegistrationTransform = reg;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetRegistrationTransform( 
  typename itk::RigidImageToImageRegistrationMethod< 
    itk::OrientedImage< float, dimensionT > >::TransformType::Pointer tfm )
{
  m_RegistrationTransform = tfm;
}

template< class pixelT, unsigned int dimensionT >
typename itk::RigidImageToImageRegistrationMethod< 
  itk::OrientedImage< float, dimensionT > >::TransformType::Pointer 
CompareCroppedROIs< pixelT, dimensionT>::
GetRegistrationTransform( void )
{
  return m_RegistrationTransform;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetOutputSize( std::vector< int > & outputSize )
{
  m_OutputSize.resize( outputSize.size() );
  for(unsigned int i=0; i<outputSize.size(); i++ )
    {
    m_OutputSize[i] = outputSize[i];
    }
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetUseMeanSquaresMetric( bool useMeanSquaresMetric )
{
  m_UseMeanSquaresMetric = useMeanSquaresMetric;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetTimeCollector( itk::TimeProbesCollectorBase * timeCollector )
{
  m_TimeCollector = timeCollector;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
SetProgressReporter( tube::CLIProgressReporter * progressReporter,
                     float progressStart, float progressRange )
{
  m_ProgressReporter = progressReporter;
  m_ProgressStart = progressStart;
  m_ProgressRange = progressRange;
}

template< class pixelT, unsigned int dimensionT >
void CompareCroppedROIs< pixelT, dimensionT>::
Update( void )
{

  if( m_Dilate > 0 )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Dilate");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart );
      }

    typedef itk::BinaryBallStructuringElement< PixelType, dimensionT >  
      BallType;
    BallType ball;
    ball.SetRadius( 1 );
    ball.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter< ImageType, ImageType, BallType >
      DilateFilterType;

    for(int r=0; r<m_Dilate; r++)
      {
      typename DilateFilterType::Pointer filter = DilateFilterType::New();
      filter->SetBackgroundValue( 0 );
      filter->SetDilateValue( 1 );
      filter->SetKernel( ball );
      filter->SetInput( m_MaskImage );
      filter->Update();
      m_MaskImage = filter->GetOutput();
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Dilate");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.2*m_ProgressRange );
      }
    }

  if( m_Erode > 0 )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Erode");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart );
      }

    typedef itk::BinaryBallStructuringElement< PixelType, dimensionT >  
      BallType;
    BallType ball;
    ball.SetRadius( 1 );
    ball.CreateStructuringElement();

    typedef itk::BinaryErodeImageFilter< ImageType, ImageType, BallType >
      ErodeFilterType;

    for(int r=0; r<m_Erode; r++)
      {
      typename ErodeFilterType::Pointer filter = ErodeFilterType::New();
      filter->SetBackgroundValue( 0 );
      filter->SetErodeValue( 1 );
      filter->SetKernel( ball );
      filter->SetInput( m_MaskImage );
      filter->Update();
      m_MaskImage = filter->GetOutput();
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Erode");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.1*m_ProgressRange );
      }
    }

  if( m_GaussianBlur > 0 )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Blur");
      }

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > 
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    for(unsigned int i=0; i<dimensionT; i++)
      {
      filter = FilterType::New();
      filter->SetInput( m_MaskImage );
      // filter->SetNormalizeAcrossScale( true );
      filter->SetSigma( m_GaussianBlur );

      filter->SetOrder( 
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );

      filter->Update();
      m_MaskImage = filter->GetOutput();
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Blur");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.3*m_ProgressRange );
      }
    }

  if( m_UseRegistration )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("RegisterROIs");
      }

    if( !m_UseRegistrationTransform )
      {
      typedef itk::NormalizeImageFilter< ImageType, ImageType > 
        NormFilterType;
  
      typename NormFilterType::Pointer norm1 = NormFilterType::New();
      norm1->SetInput( m_VolImage );
      norm1->Update();
  
      typename NormFilterType::Pointer norm2 = NormFilterType::New();
      norm2->SetInput( m_MaskImage );
      norm2->Update();
  
      typename RegistrationMethodType::Pointer reg = 
        RegistrationMethodType::New();
      reg->SetFixedImage( norm1->GetOutput() );
      reg->SetMovingImage( norm2->GetOutput() );
      reg->SetMaxIterations( 500 );
      typename RegistrationMethodType::TransformParametersScalesType scales;
      scales = reg->GetTransformParametersScales();
      for( unsigned int i=0; i<dimensionT-1; i++)
        {
        scales[i] = 1.0/0.02;
        }
      for( unsigned int i=0; i<dimensionT; i++)
        {
        scales[i + dimensionT-1] = 1.0/5.0;
        }
      reg->SetTransformParametersScales( scales );
      reg->SetUseEvolutionaryOptimization( true );
      int numSamples = 1;
      typename ImageType::SizeType imageSize = m_VolImage
                                               ->GetLargestPossibleRegion()
                                               .GetSize();
      for( unsigned int i=0; i<dimensionT; i++)
        {
        numSamples *= imageSize[i];
        }
      reg->SetNumberOfSamples( numSamples * 0.2 );
  
      if( m_ProgressReporter )
        {
        tube::CLIFilterWatcher watcher( reg, "Registration",
                                m_ProgressReporter->GetProcessInformation(),
                                0.4 * m_ProgressRange,
                                m_ProgressStart + 0.2*m_ProgressRange,
                                true );
        }

      try
        {
        reg->Update();
        m_RegistrationTransform = reg->GetTypedTransform();
        }
      catch( ... )
        {
        m_RegistrationTransform = 
          RegistrationMethodType::TransformType::New();
        m_RegistrationTransform->SetIdentity();
        }
  
      if( m_ProgressReporter )
        {
        m_ProgressReporter->Report( m_ProgressStart + 0.6*m_ProgressRange );
        }
      }
  
    typedef itk::ResampleImageFilter< ImageType, ImageType, double > 
      ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();

    typedef itk::LinearInterpolateImageFunction< ImageType, double > 
      InterpolatorType;
    typename InterpolatorType::Pointer interpolator =
      InterpolatorType::New();
    interpolator->SetInputImage( m_MaskImage );

    resampler->SetInput( m_MaskImage );
    resampler->SetInterpolator( interpolator.GetPointer() );
    resampler->SetOutputParametersFromImage( m_VolImage );
    resampler->SetTransform( m_RegistrationTransform );
    resampler->Update();
    m_MaskImage = resampler->GetOutput();

    typename ResamplerType::Pointer orgResampler =
      ResamplerType::New();
    typename InterpolatorType::Pointer orgInterpolator =
      InterpolatorType::New();
    orgInterpolator->SetInputImage( m_OrgMaskImage );
    orgResampler->SetInput( m_OrgMaskImage );
    orgResampler->SetInterpolator( orgInterpolator.GetPointer() );
    orgResampler->SetOutputParametersFromImage( m_OrgMaskImage );
    orgResampler->SetTransform( m_RegistrationTransform );
    orgResampler->Update();
    m_OrgMaskImage = orgResampler->GetOutput();

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("RegisterROIs");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.7*m_ProgressRange );
      }
    }

  if( m_OutputSize.size() > 0 )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Crop2");
      }
    typename ImageType::SizeType roiSize = m_VolImage
                                             ->GetLargestPossibleRegion()
                                             .GetSize();
    typename ImageType::SizeType lowerCropSize;
    typename ImageType::SizeType upperCropSize;

    for( unsigned int i=0; i<dimensionT; i++ )
      {
      int ti = (roiSize[i])/2 - (m_OutputSize[i]-1)/2;
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
  
      ti = roiSize[i] - (int)( lowerCropSize[i] + m_OutputSize[i] );
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
    typename CropFilterType::Pointer cropVolumeFilter =
      CropFilterType::New();
    typename CropFilterType::Pointer cropMaskFilter =
      CropFilterType::New();
  
    cropVolumeFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropVolumeFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropVolumeFilter->SetInput( m_VolImage );
    cropVolumeFilter->Update();
    m_VolImage = cropVolumeFilter->GetOutput();
  
    cropMaskFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropMaskFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropMaskFilter->SetInput( m_MaskImage );
    cropMaskFilter->Update();
    m_MaskImage = cropMaskFilter->GetOutput();
    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Crop2");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.8*m_ProgressRange );
      }
    }
  
  if( m_Normalize )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Normalize");
      }
  
    typedef itk::ImageRegionIterator< ImageType > IterType;
    IterType volIter( m_VolImage, m_VolImage->GetLargestPossibleRegion() );
    IterType maskIter( m_MaskImage,
                       m_MaskImage->GetLargestPossibleRegion() );
    IterType orgMaskIter( m_OrgMaskImage,
                          m_OrgMaskImage->GetLargestPossibleRegion() );
  
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
  
    double binCount = 0;
    int bin;
    volIter.GoToBegin();
    maskIter.GoToBegin();
    orgMaskIter.GoToBegin();
    while( !volIter.IsAtEnd() )
      {
      if( orgMaskIter.Get() == m_Foreground )
        {
        bin = (int)( (((double)(volIter.Get())-volMin)/(volMax-volMin)) 
                               * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++volFgHisto[ bin ];
        ++binCount;
        if( bin > 0 )
          {
          volFgHisto[ bin-1 ] += 0.5;
          binCount += 0.5;
          }
        if( bin < numBins - 1 )
          {
          volFgHisto[ bin+1 ] += 0.5;
          binCount += 0.5;
          }
  
        bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin))
                               * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++maskFgHisto[ bin ];
        ++binCount;
        if( bin > 0 )
          {
          maskFgHisto[ bin-1 ] += 0.5;
          binCount += 0.5;
          }
        if( bin < numBins - 1 )
          {
          maskFgHisto[ bin+1 ] += 0.5;
          binCount += 0.5;
          }
        }
      else
        {
        bin = (int)( (((double)(volIter.Get())-volMin)/(volMax-volMin)) 
                               * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++volBgHisto[ bin ];
        ++binCount;
        if( bin > 0 )
          {
          volBgHisto[ bin-1 ] += 0.5;
          binCount += 0.5;
          }
        if( bin < numBins - 1 )
          {
          volBgHisto[ bin+1 ] += 0.5;
          binCount += 0.5;
          }
  
        bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin))
                               * (numBins-1) );
        if( bin < 0 )
          {
          bin = 0;
          }
        else if( bin >= numBins )
          {
          bin = numBins - 1;
          }
        ++maskBgHisto[ bin ];
        ++binCount;
        if( bin > 0 )
          {
          maskBgHisto[ bin-1 ] += 0.5;
          binCount += 0.5;
          }
        if( bin < numBins - 1 )
          {
          maskBgHisto[ bin+1 ] += 0.5;
          binCount += 0.5;
          }
        }
  
      ++volIter;
      ++maskIter;
      ++orgMaskIter;
      }
    binCount /= 2;

    /*  Debugguing - printout bins
    for( bin=0; bin<numBins; bin++ )
      {
      std::cout << volBgHisto[bin] << " : " << volFgHisto[bin] << " - "
                << maskBgHisto[bin] << " : " << maskFgHisto[bin]
                << std::endl;
      } */
  
    float maxVFg = 0;
    int maxVFgI = -1;
    float maxVBg = 0;
    int maxVBgI = -1;
    float maxMFg = 0;
    int maxMFgI = -1;
    float maxMBg = 0;
    int maxMBgI = -1;
    /*
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
    */
    double clip = 0.05;
    for( int i=0; i<numBins; i++ )
      {
      maxVBg += volBgHisto[i];
      if( maxVBg > clip * binCount )
        {
        maxVBgI = i;
        break;
        }
      }
    for( int i=numBins-1; i>=0; i-- )
      {
      maxVFg += volFgHisto[i];
      if( maxVFg > clip * binCount )
        {
        maxVFgI = i;
        break;
        }
      }
    clip = 0.02;
    for( int i=0; i<numBins; i++ )
      {
      maxMBg += maskBgHisto[i];
      if( maxMBg > clip * binCount )
        {
        maxMBgI = i;
        break;
        }
      }
    for( int i=numBins-1; i>=0; i-- )
      {
      maxMFg += maskFgHisto[i];
      if( maxMFg > clip * binCount )
        {
        maxMFgI = i;
        break;
        }
      }
  
    maxVFg = (double)(maxVFgI)/(numBins-1) * (volMax - volMin) + volMin;
    maxVBg = (double)(maxVBgI)/(numBins-1) * (volMax - volMin) + volMin;
    maxMFg = (double)(maxMFgI)/(numBins-1) * (maskMax - maskMin) + maskMin;
    maxMBg = (double)(maxMBgI)/(numBins-1) * (maskMax - maskMin) + maskMin;

    /*
    std::cout << "VBg = " << maxVBg <<  " - "
              << "VFg = " << maxVFg << std::endl;
    std::cout << "MBg = " << maxMBg <<  " - "
              << "MFg = " << maxMFg << std::endl;
              */
   
    if( maxVFgI != -1 &&
        maxVBgI != -1 &&
        maxMFgI != -1 &&
        maxMBgI != -1 )
      {
      typedef itk::ShiftScaleImageFilter< ImageType, ImageType > 
        FilterType;
      typename FilterType::Pointer filter = FilterType::New();
   
      double scale = (maxVFg - maxVBg) / (maxMFg - maxMBg);
      double shift = maxVBg/scale - maxMBg;
      filter = FilterType::New();
      filter->SetInput( m_MaskImage );
      filter->SetShift( shift );
      filter->SetScale( scale );
   
      filter->Update();
      m_MaskImage = filter->GetOutput();
      }
  
    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Normalize");
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.9*m_ProgressRange );
      }
    }

  /* Compute similarity */
  if( true )
    {
    if( m_TimeCollector )
      {
      m_TimeCollector->Start("Match metric");
      }

    typedef itk::NormalizeImageFilter< ImageType, ImageType > 
      NormFilterType;
  
    typename NormFilterType::Pointer norm1 = NormFilterType::New();
    norm1->SetInput( m_VolImage );
    norm1->Update();
    typename ImageType::Pointer image1 = norm1->GetOutput();
  
    typename NormFilterType::Pointer norm2 = NormFilterType::New();
    norm2->SetInput( m_MaskImage );
    norm2->Update();
    typename ImageType::Pointer image2 = norm2->GetOutput();
  
    typedef itk::AffineTransform< double, dimensionT > TransformType;
    typename TransformType::Pointer transform = TransformType::New();
    transform->SetIdentity();
  
    typedef itk::LinearInterpolateImageFunction< ImageType, double > 
      InterpolatorType;
    typename InterpolatorType::Pointer interpolator =
      InterpolatorType::New();
    interpolator->SetInputImage( image2 );
  
    typedef itk::ImageToImageMetric< ImageType, ImageType >   MetricType;
    typename MetricType::Pointer metric;
  
    if( m_UseCorrelationMetric )
      {
      typedef itk::NormalizedCorrelationImageToImageMetric< ImageType,
        ImageType > CorMetricType;
      metric = CorMetricType::New();
      }
    else if( m_UseMeanSquaresMetric )
      {
      typedef itk::MeanSquaresImageToImageMetric< ImageType,
        ImageType > MSMetricType;
      metric = MSMetricType::New();
      }
    else
      {
      typedef itk::MutualInformationImageToImageMetric< ImageType,
        ImageType > MIMetricType;
      metric = MIMetricType::New();
      }
  
    typename ImageType::SizeType size = image1
                                        ->GetLargestPossibleRegion()
                                        .GetSize();
  
    metric->SetFixedImage( image1 );
    metric->SetMovingImage( image2 );
    metric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
    metric->SetTransform( transform );
    metric->SetInterpolator( interpolator );
    int numSamples = 1;
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      numSamples *= size[i];
      }
    metric->SetNumberOfSpatialSamples( numSamples*m_SamplingRate );
    metric->Initialize();
    metric->MultiThreadingInitialize();
  
    if( !m_UseCorrelationMetric && !m_UseMeanSquaresMetric )
      {
      m_GoF = metric->GetValue( transform->GetParameters() );
      }
    else
      {
      m_GoF = -metric->GetValue( transform->GetParameters() );
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop("Match metric");
      }
    }

  if( m_ProgressReporter )
    {
    m_ProgressReporter->Report( m_ProgressStart + 1.0*m_ProgressRange );
    }
}

template< class pixelT, unsigned int dimensionT >
float CompareCroppedROIs< pixelT, dimensionT>::
GetGoodnessOfFit( void )
{
  return m_GoF;
}

}

#endif
