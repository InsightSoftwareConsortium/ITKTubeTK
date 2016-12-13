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

#ifndef __tubeCompareImageWithPrior_hxx
#define __tubeCompareImageWithPrior_hxx

#include "tubeCompareImageWithPrior.h"

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCropImageFilter.h>
#include <itkImage.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkNormalizeImageFilter.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkRigidImageToImageRegistrationMethod.h>
#include <itkShiftScaleImageFilter.h>
#include <itkTimeProbesCollectorBase.h>

#include <sstream>

namespace tube
{

template< class TPixel, unsigned int TDimension >
CompareImageWithPrior< TPixel, TDimension>
::CompareImageWithPrior( void )
{
  m_VolImage = NULL;
  m_MaskImage = NULL;
  m_OutputVolImage = NULL;
  m_OutputMaskImage = NULL;
  m_MetricMask = NULL;
  m_Foreground = 0;
  m_Background = 0;
  m_Erode = 4;
  m_Dilate = 2;
  m_GaussianBlur = 10;
  m_UseRegistration = true;
  m_UseRegistrationOptimization = true;
  m_UseRegistrationTransform = false;
  m_RegistrationTransform = RegistrationMethodType::TransformType::New();
  m_RegistrationTransform->SetIdentity();
  m_Normalize = true;
  m_BoundarySize.resize( 0 );
  m_SamplingRate = 0.2;
  m_Seed = 0;

  m_TimeCollector = NULL;
  m_ProgressReporter = NULL;
  m_ProgressStart = 0;
  m_ProgressRange = 1.0;

  m_GoF = 0;
}

template< class TPixel, unsigned int TDimension >
CompareImageWithPrior< TPixel, TDimension>
::~CompareImageWithPrior( void )
{
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetInput( typename ImageType::Pointer volImage )
{
  m_VolImage = volImage;
}

template< class TPixel, unsigned int TDimension >
typename itk::Image< float, TDimension >::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetInput( void )
{
  return m_VolImage;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetMaskImage( typename ImageType::Pointer maskImage )
{
  m_MaskImage = maskImage;
}

template< class TPixel, unsigned int TDimension >
typename itk::Image< float, TDimension >::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetMaskImage( void )
{
  return m_MaskImage;
}

template< class TPixel, unsigned int TDimension >
typename itk::Image< float, TDimension >::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetOutput( void )
{
  return m_OutputVolImage;
}

template< class TPixel, unsigned int TDimension >
typename itk::Image< float, TDimension >::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetOutputMaskImage( void )
{
  return m_OutputMaskImage;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetMetricMask( typename ImageType::Pointer metricMask )
{
  m_MetricMask = metricMask;
}

template< class TPixel, unsigned int TDimension >
typename itk::Image< float, TDimension >::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetMetricMask( void )
{
  return m_MetricMask;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetForeground( float foreground )
{
  m_Foreground = foreground;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetBackground( float background )
{
  m_Background = background;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetErode( int erode )
{
  m_Erode = erode;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetDilate( int dilate )
{
  m_Dilate = dilate;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetGaussianBlur( float gaussianBlur )
{
  m_GaussianBlur = gaussianBlur;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetSamplingRate( float samplingRate )
{
  m_SamplingRate = samplingRate;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetSeed( unsigned int seed )
{
  m_Seed = seed;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetNormalize( bool normalize )
{
  m_Normalize = normalize;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetUseRegistration( bool reg )
{
  m_UseRegistration = reg;
}

template< class TPixel, unsigned int TDimension >
bool
CompareImageWithPrior< TPixel, TDimension>
::GetUseRegistration( void )
{
  return m_UseRegistration;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetUseRegistrationTransform( bool reg )
{
  m_UseRegistrationTransform = reg;
}

template< class TPixel, unsigned int TDimension >
bool
CompareImageWithPrior< TPixel, TDimension>
::GetUseRegistrationTransform( void )
{
  return m_UseRegistrationTransform;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetUseRegistrationOptimization( bool reg )
{
  m_UseRegistrationOptimization = reg;
}

template< class TPixel, unsigned int TDimension >
bool
CompareImageWithPrior< TPixel, TDimension>
::GetUseRegistrationOptimization( void )
{
  return m_UseRegistrationOptimization;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetRegistrationTransform( 
  typename itk::RigidImageToImageRegistrationMethod<
    itk::Image< float, TDimension > >::TransformType::Pointer tfm )
{
  m_RegistrationTransform = tfm;
}

template< class TPixel, unsigned int TDimension >
typename itk::RigidImageToImageRegistrationMethod<
  itk::Image< float, TDimension > >::TransformType::Pointer
CompareImageWithPrior< TPixel, TDimension>
::GetRegistrationTransform( void )
{
  return m_RegistrationTransform;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetBoundarySize( std::vector< int > & boundarySize )
{
  m_BoundarySize.resize( boundarySize.size() );
  for( unsigned int i=0; i<boundarySize.size(); i++ )
    {
    m_BoundarySize[i] = boundarySize[i];
    }
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetTimeCollector( itk::TimeProbesCollectorBase * timeCollector )
{
  m_TimeCollector = timeCollector;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::SetProgressReporter( tube::CLIProgressReporter * progressReporter,
  float progressStart,
  float progressRange )
{
  m_ProgressReporter = progressReporter;
  m_ProgressStart = progressStart;
  m_ProgressRange = progressRange;
}

template< class TPixel, unsigned int TDimension >
void
CompareImageWithPrior< TPixel, TDimension>
::Update( void )
{
  m_OutputVolImage = m_VolImage;
  m_OutputMaskImage = m_MaskImage;
  typename ImageType::Pointer regMaskImage = m_MaskImage;

  if( m_Erode > 0 )
    {
    std::cout << "Eroding = " << m_Erode << std::endl;

    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "Erode" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart );
      }

    typedef itk::BinaryBallStructuringElement< PixelType, TDimension >
      BallType;
    BallType ball;
    ball.SetRadius( 1 );
    ball.CreateStructuringElement();

    typedef itk::BinaryErodeImageFilter< ImageType, ImageType, BallType >
      ErodeFilterType;

    for( int r=0; r<m_Erode; r++ )
      {
      typename ErodeFilterType::Pointer filter = ErodeFilterType::New();
      filter->SetBackgroundValue( m_Background );
      filter->SetErodeValue( m_Foreground );
      filter->SetKernel( ball );
      filter->SetInput( m_OutputMaskImage );
      filter->Update();
      m_OutputMaskImage = filter->GetOutput();

      regMaskImage = m_OutputMaskImage;
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "Erode" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.1*m_ProgressRange );
      }
    }

  if( m_Dilate > 0 )
    {
    std::cout << "Dilating = " << m_Dilate << std::endl;

    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "Dilate" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart );
      }

    typedef itk::BinaryBallStructuringElement< PixelType, TDimension >
      BallType;
    BallType ball;
    ball.SetRadius( 1 );
    ball.CreateStructuringElement();

    typedef itk::BinaryDilateImageFilter< ImageType, ImageType, BallType >
      DilateFilterType;

    for( int r=0; r<m_Dilate; r++ )
      {
      typename DilateFilterType::Pointer filter = DilateFilterType::New();
      filter->SetBackgroundValue( m_Background );
      filter->SetDilateValue( m_Foreground );
      filter->SetKernel( ball );
      filter->SetInput( m_OutputMaskImage );
      filter->Update();
      m_OutputMaskImage = filter->GetOutput();

      regMaskImage = m_OutputMaskImage;
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "Dilate" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.2*m_ProgressRange );
      }
    }

  if( m_GaussianBlur > 0 )
    {
    std::cout << "Blurring = " << m_GaussianBlur << std::endl;
    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "Blur" );
      }

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    for( unsigned int i=0; i<TDimension; i++ )
      {
      filter = FilterType::New();
      filter->SetInput( m_OutputMaskImage );
      filter->SetSigma( m_GaussianBlur );

      filter->SetOrder( 
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );

      filter->Update();
      m_OutputMaskImage = filter->GetOutput();
      }

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "Blur" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.3*m_ProgressRange );
      }
    }

  if( m_UseRegistration )
    {
    std::cout << "Registering" << std::endl;
    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "RegisterROIs" );
      }

    if( m_UseRegistrationOptimization )
      {
      typedef itk::NormalizeImageFilter< ImageType, ImageType >
        NormFilterType;

      typename NormFilterType::Pointer norm1 = NormFilterType::New();
      norm1->SetInput( m_OutputVolImage );
      norm1->Update();

      typename NormFilterType::Pointer norm2 = NormFilterType::New();
      norm2->SetInput( m_OutputMaskImage );
      norm2->Update();

      typename RegistrationMethodType::Pointer reg =
        RegistrationMethodType::New();
      reg->SetFixedImage( norm1->GetOutput() );
      reg->SetMovingImage( norm2->GetOutput() );
      reg->SetMaxIterations( 500 );
      typename RegistrationMethodType::TransformParametersScalesType scales;
      scales = reg->GetTransformParametersScales();
      for( unsigned int i=0; i<TDimension-1; i++ )
        {
        scales[i] = 1.0/0.02;
        }
      for( unsigned int i=0; i<TDimension; i++ )
        {
        scales[i + TDimension-1] = 1.0/10.0;
        }
      reg->SetTransformParametersScales( scales );
      reg->SetUseEvolutionaryOptimization( true );
      int numSamples = 1;
      if( m_UseRegistrationTransform )
        {
        reg->SetInitialTransformParameters( 
          m_RegistrationTransform->GetParameters() );
        }
      typename ImageType::SizeType imageSize = m_OutputVolImage
                                               ->GetLargestPossibleRegion()
                                               .GetSize();
      for( unsigned int i=0; i<TDimension; i++ )
        {
        numSamples *= imageSize[i];
        }
      reg->SetNumberOfSamples( numSamples * m_SamplingRate );
      if( m_Seed > 0 )
        {
        reg->SetRandomNumberSeed( m_Seed );
        }

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
        tube::WarningMessage( 
          "Exception thrown during registration. Compensating..." );
        m_RegistrationTransform =
          RegistrationMethodType::TransformType::New();
        m_RegistrationTransform->SetIdentity();
        }

      std::stringstream str( std::stringstream::in | std::stringstream::out );
      str << "Registration params =";
      for( unsigned int i=0;
           i<m_RegistrationTransform->GetNumberOfParameters();
           i++ )
        {
        str << " " << m_RegistrationTransform->GetParameters()[i];
        }
      tube::InfoMessage( str.str() );
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
    interpolator->SetInputImage( m_OutputMaskImage );

    resampler->SetInput( m_OutputMaskImage );
    resampler->SetInterpolator( interpolator.GetPointer() );
    resampler->SetOutputParametersFromImage( m_OutputVolImage );
    resampler->SetTransform( m_RegistrationTransform );
    resampler->Update();
    m_OutputMaskImage = resampler->GetOutput();

    typedef itk::NearestNeighborInterpolateImageFunction< ImageType,
      double > NNInterpolatorType;
    typename ResamplerType::Pointer orgResampler =
      ResamplerType::New();
    typename NNInterpolatorType::Pointer orgInterpolator =
      NNInterpolatorType::New();
    orgInterpolator->SetInputImage( regMaskImage );
    orgResampler->SetInput( regMaskImage );
    orgResampler->SetInterpolator( orgInterpolator.GetPointer() );
    orgResampler->SetOutputParametersFromImage( m_OutputVolImage );
    orgResampler->SetTransform( m_RegistrationTransform );
    orgResampler->Update();
    regMaskImage = orgResampler->GetOutput();

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "RegisterROIs" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.7*m_ProgressRange );
      }
    }

  if( m_BoundarySize.size() > 0 )
    {
    std::cout << "Cropping = " << m_BoundarySize[0] << std::endl;

    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "Crop2" );
      }
    typename ImageType::SizeType roiSize = m_OutputVolImage
                                             ->GetLargestPossibleRegion()
                                             .GetSize();
    typename ImageType::SizeType lowerCropSize;
    typename ImageType::SizeType upperCropSize;

    for( unsigned int i=0; i<TDimension; i++ )
      {
      if( m_BoundarySize[i] < 0 )
        {
        lowerCropSize[i] = 0;
        upperCropSize[i] = 0;
        }
      else if( m_BoundarySize[i] >= ( int )( roiSize[i] )/2 )
        {
        lowerCropSize[i] = roiSize[i]/2-1;
        upperCropSize[i] = roiSize[i]/2-1;
        }
      else
        {
        lowerCropSize[i] = m_BoundarySize[i];
        upperCropSize[i] = m_BoundarySize[i];
        }
      }

    typedef itk::CropImageFilter< ImageType, ImageType > CropFilterType;
    typename CropFilterType::Pointer cropInputFilter =
      CropFilterType::New();
    typename CropFilterType::Pointer cropMaskFilter =
      CropFilterType::New();
    typename CropFilterType::Pointer cropRegMaskFilter =
      CropFilterType::New();

    cropInputFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropInputFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropInputFilter->SetInput( m_OutputVolImage );
    cropInputFilter->Update();
    m_OutputVolImage = cropInputFilter->GetOutput();

    cropMaskFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropMaskFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropMaskFilter->SetInput( m_OutputMaskImage );
    cropMaskFilter->Update();
    m_OutputMaskImage = cropMaskFilter->GetOutput();

    cropRegMaskFilter->SetLowerBoundaryCropSize( lowerCropSize );
    cropRegMaskFilter->SetUpperBoundaryCropSize( upperCropSize );
    cropRegMaskFilter->SetInput( regMaskImage );
    cropRegMaskFilter->Update();
    regMaskImage = cropRegMaskFilter->GetOutput();

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "Crop2" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.8*m_ProgressRange );
      }
    }

  if( m_Normalize )
    {
    std::cout << "Normalizing" << std::endl;

    if( m_TimeCollector )
      {
      m_TimeCollector->Start( "Normalize" );
      }

    typedef itk::ImageRegionIterator< ImageType > IterType;
    IterType volIter( m_OutputVolImage,
      m_OutputVolImage->GetLargestPossibleRegion() );
    IterType maskIter( m_OutputMaskImage,
      m_OutputMaskImage->GetLargestPossibleRegion() );
    IterType regMaskIter( regMaskImage,
      regMaskImage->GetLargestPossibleRegion() );

    int countVolFg = 0;
    double meanVolFg = 0;
    int countVolBg = 0;
    double meanVolBg = 0;
    int countMaskFg = 0;
    double meanMaskFg = 0;
    int countMaskBg = 0;
    double meanMaskBg = 0;

    volIter.GoToBegin();
    maskIter.GoToBegin();
    regMaskIter.GoToBegin();
    while( !volIter.IsAtEnd() )
      {
      if( regMaskIter.Get() == m_Foreground )
        {
        meanVolFg += volIter.Get();
        meanMaskFg += maskIter.Get();
        ++countVolFg;
        ++countMaskFg;
        }
      else if( regMaskIter.Get() == m_Background )
        {
        meanVolBg += volIter.Get();
        meanMaskBg += maskIter.Get();
        ++countVolBg;
        ++countMaskBg;
        }

      ++volIter;
      ++maskIter;
      ++regMaskIter;
      }

    if( countVolFg != 0 )
      {
      meanVolFg /= countVolFg;
      }
    if( countVolBg != 0 )
      {
      meanVolBg /= countVolBg;
      }
    if( countMaskFg != 0 )
      {
      meanMaskFg /= countMaskFg;
      }
    if( countMaskBg != 0 )
      {
      meanMaskBg /= countMaskBg;
      }

    std::cout << "VBg = " << meanVolBg <<  " - "
              << "VFg = " << meanVolFg << std::endl;
    std::cout << "MBg = " << meanMaskBg <<  " - "
              << "MFg = " << meanMaskFg << std::endl;

    double scale = ( meanVolFg - meanVolBg ) / ( meanMaskFg - meanMaskBg );
    double shift = meanVolBg/scale - meanMaskBg;

     typedef itk::ShiftScaleImageFilter< ImageType, ImageType >
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();
    filter->SetInput( m_OutputMaskImage );
    filter->SetShift( shift );
    filter->SetScale( scale );

    std::cout << "Scale = " << scale << std::endl;
    std::cout << "Shift = " << shift << std::endl;

    filter->Update();
    m_OutputMaskImage = filter->GetOutput();

    if( m_TimeCollector )
      {
      m_TimeCollector->Stop( "Normalize" );
      }
    if( m_ProgressReporter )
      {
      m_ProgressReporter->Report( m_ProgressStart + 0.9*m_ProgressRange );
      }
    }

  /* Compute similarity */
  std::cout << "Computing similarity" << std::endl;
  if( m_TimeCollector )
    {
    m_TimeCollector->Start( "Match metric" );
    }

  itk::ImageRegionIteratorWithIndex< ImageType > volIter( 
    m_OutputVolImage, m_OutputVolImage->GetLargestPossibleRegion() );
  itk::ImageRegionIteratorWithIndex< ImageType > maskIter( 
    m_OutputMaskImage, m_OutputMaskImage->GetLargestPossibleRegion() );

  m_GoF = 0;
  unsigned int count = 0;
  typename ImageType::PointType pnt;
  typename ImageType::IndexType indx;
  while( !volIter.IsAtEnd() )
    {
    bool valid = true;
    if( m_MetricMask.IsNotNull() )
      {
      valid = false;
      m_OutputVolImage->TransformIndexToPhysicalPoint( volIter.GetIndex(),
        pnt );
      if( m_MetricMask->TransformPhysicalPointToIndex( pnt, indx ) )
        {
        if( m_MetricMask->GetPixel( indx ) != 0 )
          {
          valid = true;
          }
        }
      }
    if( valid )
      {
      double tf = volIter.Get() - maskIter.Get();
      m_GoF += tf * tf;
      ++count;
      }
    ++volIter;
    ++maskIter;
    }

  if( m_GoF>0 && count>0 )
    {
    m_GoF = -std::sqrt( m_GoF/count );
    }

  if( m_TimeCollector )
    {
    m_TimeCollector->Stop( "Match metric" );
    }
  if( m_ProgressReporter )
    {
    m_ProgressReporter->Report( m_ProgressStart + 0.95*m_ProgressRange );
    }

  std::cout << "gof = " << m_GoF << std::endl;

  if( m_ProgressReporter )
    {
    m_ProgressReporter->Report( m_ProgressStart + 1.0*m_ProgressRange );
    }
}

template< class TPixel, unsigned int TDimension >
float
CompareImageWithPrior< TPixel, TDimension>
::GetGoodnessOfFit( void )
{
  return m_GoF;
}

} // End namespace tube

#endif // End !defined( __tubeCompareImageWithPrior_hxx )
