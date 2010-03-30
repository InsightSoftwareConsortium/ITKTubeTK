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

#include "tubeSubImageGenerator.h"

#include "itkImageRegionConstIterator.h"
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

namespace tube
{

template< class pixelT, unsigned int dimensionT>
SubImageGenerator<pixelT,dimensionT>
::SubImageGenerator()
  : m_RoiCenter(),
    m_RoiSize(),
    m_InputVolume(NULL),
    m_InputMask(NULL),
    m_OutputVolume(NULL),
    m_OutputMask(NULL)
{
}


  /// Full populated constructor
template< class pixelT, unsigned int dimensionT>
SubImageGenerator<pixelT,dimensionT>
::SubImageGenerator( std::vector<int> roiCenter,
                     std::vector<int> roiSize,
                     typename ImageType::Pointer inputVolume,
                     typename ImageType::Pointer inputMask )
  : m_RoiCenter(roiCenter),
    m_RoiSize(roiSize),
    m_InputVolume(inputVolume),
    m_InputMask(inputMask),
    m_OutputVolume(NULL),
    m_OutputMask(NULL)
{
}

  /// Default Destructor
template< class pixelT, unsigned int dimensionT>
SubImageGenerator<pixelT,dimensionT>
::~SubImageGenerator()
{
}

  /// Update function for doing the processing and producing the output 
  /// volumes.
template< class pixelT, unsigned int dimensionT>
void
SubImageGenerator<pixelT,dimensionT>
::Update()
{
  /**typedef pixelT PixelType;

  typename ImageType::Pointer curVolume = m_InputVolume;
  typename ImageType::Pointer curMask = m_InputMask;

  typename ImageType::SizeType inputSize = 
    curVolume->GetLargestPossibleRegion().GetSize();
  typename ImageType::SizeType lowerCropSize;
  typename ImageType::SizeType upperCropSize;

  std::vector<int> roiCenter = m_RoiCenter;
  std::vector<int> roiSize = m_RoiSize;
  std::vector<int> outputSize = m_RoiSize;
  std::vector<int>::iterator oSizeItr;
  for( oSizeItr = outputSize.begin(); oSizeItr != outputSize.end(); 
       ++oSizeItr )
    {
    *oSizeItr = *oSizeItr - 20;
    }

  // Crop input images to ROI 
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

  typename ImageType::Pointer orgMask = curMask;

  float foreground = 255;
  float background = 0;
  float erode = 7;
  float gaussianBlur = 10;

  if( foreground != 1 || background != 0 )
    {

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

    }
  
  if( erode > 0 )
    {
    typedef itk::BinaryBallStructuringElement<PixelType, dimensionT > BallType;
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
    }

  if( gaussianBlur > 0 )
    {

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > 
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    for(unsigned int i=0; i<dimensionT; i++)
      {
      filter = FilterType::New();
      filter->SetInput( curMask );
      filter->SetSigma( gaussianBlur );

      filter->SetOrder( 
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      filter->SetDirection( i );

      filter->Update();
      curMask = filter->GetOutput();
      }

    }

    {
    typedef itk::NormalizeImageFilter< ImageType, ImageType > NormFilterType;

    typename NormFilterType::Pointer norm1 = NormFilterType::New();
    norm1->SetInput( curVolume );
    norm1->Update();

    typename NormFilterType::Pointer norm2 = NormFilterType::New();
    norm2->SetInput( curMask );
    norm2->Update();

    typedef itk::RigidImageToImageRegistrationMethod< ImageType >
                                                    RegistrationMethodType;
    typename RegistrationMethodType::Pointer reg = 
      RegistrationMethodType::New();
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

    typedef itk::ResampleImageFilter< ImageType, ImageType, double > 
      ResamplerType;
    typename ResamplerType::Pointer resampler = ResamplerType::New();

    typedef itk::LinearInterpolateImageFunction< ImageType, double > 
      InterpolatorType;
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
    interpolator->SetInputImage( curMask );

    resampler->SetInput( curMask );
    resampler->SetInterpolator( interpolator.GetPointer() );
    resampler->SetOutputParametersFromImage( curMask );
    resampler->SetTransform( reg->GetTypedTransform() );
    resampler->Update();
    curMask = resampler->GetOutput();

    typename ResamplerType::Pointer orgResampler = ResamplerType::New();
    typename InterpolatorType::Pointer orgInterpolator = 
      InterpolatorType::New();
    orgInterpolator->SetInputImage( orgMask );
    orgResampler->SetInput( orgMask );
    orgResampler->SetInterpolator( orgInterpolator.GetPointer() );
    orgResampler->SetOutputParametersFromImage( orgMask );
    orgResampler->SetTransform( reg->GetTypedTransform() );
    orgResampler->Update();
    orgMask = orgResampler->GetOutput();
    }

  if( outputSize.size() > 0 )
    {
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
    }
  
  typedef itk::CropImageFilter< ImageType, ImageType > CropFilterType;
  cropVolumeFilter = CropFilterType::New();
  cropMaskFilter = CropFilterType::New();
  
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
      bin = (int)((((double)(volIter.Get())-volMin)/(volMax-volMin)) * 
                  (numBins-1) );
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
      
      bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin)) * 
                   (numBins-1) );
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
      bin = (int)( (((double)(volIter.Get())-volMin)/(volMax-volMin)) * 
                   (numBins-1) );
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
      
      bin = (int)( (((double)(maskIter.Get())-maskMin)/(maskMax-maskMin)) * 
                   (numBins-1) );
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

  m_OutputVolume = curVolume;
  m_OutputMask = curMask;**/

  /*******/
  // Old method 
  // Create blank output images of the appropriate size
  m_OutputVolume = ImageType::New();
  m_OutputMask = ImageType::New();
  typename ImageType::SizeType size;
  typename ImageType::IndexType start;
  typename ImageType::IndexType inputIndex;
  for( unsigned int i = 0; i < dimensionT; ++i )
    {
    size[i] = m_RoiSize[i];
    inputIndex[i] = m_RoiCenter[i] - (m_RoiSize[i]/2); // proper roi start
    start[i] = 0;
    }
  typename ImageType::RegionType region;
  region.SetSize( size );
  region.SetIndex( start );
  m_OutputVolume->CopyInformation( m_InputVolume );
  m_OutputVolume->SetRegions( region );
  m_OutputMask->CopyInformation( m_InputMask );
  m_OutputMask->SetRegions( region );
  m_OutputVolume->Allocate();
  m_OutputMask->Allocate();
  m_OutputVolume->FillBuffer(0);
  m_OutputMask->FillBuffer(0);

  // Iterate through the input and mask and populate the outputs
  typename ImageType::RegionType inputRegion;
  inputRegion.SetSize( size );
  inputRegion.SetIndex( inputIndex );
  typedef itk::ImageRegionConstIterator<ImageType> ConstIteratorType;
  typedef itk::ImageRegionIterator<ImageType> IteratorType;
  ConstIteratorType inputItr( m_InputVolume, inputRegion );
  ConstIteratorType maskItr( m_InputMask, inputRegion );
  IteratorType outputItr( m_OutputVolume, region );
  IteratorType outputMaskItr( m_OutputMask, region );
  inputItr.GoToBegin();
  maskItr.GoToBegin();
  outputItr.GoToBegin();
  outputMaskItr.GoToBegin();
  while( !inputItr.IsAtEnd() && !maskItr.IsAtEnd() &&
         !outputItr.IsAtEnd() && !outputMaskItr.IsAtEnd() )
    {
    outputItr.Set( inputItr.Get() );
    outputMaskItr.Set( maskItr.Get() );
    ++inputItr;
    ++maskItr;
    ++outputItr;
    ++outputMaskItr;
    }

}

template< class pixelT, unsigned int dimensionT>
void
SubImageGenerator<pixelT,dimensionT>
::SetRoiCenter( std::vector<int> roiCenter )
{
  m_RoiCenter = roiCenter;
}
 
template< class pixelT, unsigned int dimensionT>
void
SubImageGenerator<pixelT,dimensionT>
::SetRoiSize( std::vector<int> roiSize )
{
  m_RoiSize = roiSize;
}
 
template< class pixelT, unsigned int dimensionT>
void
SubImageGenerator<pixelT,dimensionT>
::SetInputVolume( typename ImageType::Pointer inputVolume )
{
  m_InputVolume = inputVolume;
}

template< class pixelT, unsigned int dimensionT>
void
SubImageGenerator<pixelT,dimensionT>
::SetInputMask( typename ImageType::Pointer inputMask )
{
  m_InputMask = inputMask;
}

template< class pixelT, unsigned int dimensionT>
typename SubImageGenerator<pixelT,dimensionT>::ImageType::Pointer
SubImageGenerator<pixelT,dimensionT>
::GetOutputVolume()
{
  return m_OutputVolume;
}

template< class pixelT, unsigned int dimensionT>
typename SubImageGenerator<pixelT,dimensionT>::ImageType::Pointer
SubImageGenerator<pixelT,dimensionT>
::GetOutputMask()
{
  return m_OutputMask;
}

} // end namespace tube
