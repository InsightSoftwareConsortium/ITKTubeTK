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
#ifndef __tubeSubImageGenerator_txx
#define __tubeSubImageGenerator_txx

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
  typedef itk::ImageRegionIterator<ImageType>      IteratorType;
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

#endif
