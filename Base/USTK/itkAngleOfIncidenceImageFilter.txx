/*=========================================================================

Library:   TubeTK

Copyright 2012 Kitware Inc. 28 Corporate Drive,
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
#ifndef __itkAngleOfIncidenceImageFilter_hxx
#define __itkAngleOfIncidenceImageFilter_hxx

#include "itkAngleOfIncidenceImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkProgressReporter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

namespace itk
{
/**
 * Constructor
 */
template< class TInputImage, class TOutputImage >
AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::AngleOfIncidenceImageFilter()
{
  m_UltrasoundOrigin.Fill(0);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::GenerateData(void)
{
  typedef typename TInputImage::SizeType          InputImageSizeType;
  typedef typename TInputImage::PointType         InputImagePointType;
  typedef typename TInputImage::SpacingType       InputImageSpacingType;
  typedef ImageRegionConstIterator< TInputImage > InputImageConstIteratorType;

  typedef typename TOutputImage::IndexType    ImageIndexType;
  typedef typename TOutputImage::PixelType    PixelType;
  typedef ImageRegionIterator< TOutputImage > OutputImageIteratorType;


  typename TInputImage::ConstPointer inputImagePtr(
    dynamic_cast< const TInputImage  * >(
      this->ProcessObject::GetInput(0) ) );
  typename TOutputImage::Pointer outputImagePtr(
    dynamic_cast< TOutputImage * >(
      this->ProcessObject::GetOutput(0) ) );

  VectorType originInput;
  originInput.Fill(0.0);
  outputImagePtr->SetOrigin(originInput);
  outputImagePtr->SetSpacing( inputImagePtr->GetSpacing() );
  outputImagePtr->SetDirection( inputImagePtr->GetDirection() );

  outputImagePtr->SetRequestedRegion( inputImagePtr->GetRequestedRegion() );
  outputImagePtr->SetBufferedRegion( inputImagePtr->GetBufferedRegion() );
  outputImagePtr->SetLargestPossibleRegion( inputImagePtr->GetLargestPossibleRegion() );
  outputImagePtr->Allocate();
  outputImagePtr->FillBuffer(0);

  InputImageConstIteratorType inputIt( inputImagePtr, inputImagePtr->GetLargestPossibleRegion() );
  OutputImageIteratorType     outputIt( outputImagePtr, outputImagePtr->GetLargestPossibleRegion() );

  outputIt.GoToBegin();
  inputIt.GoToBegin();

  InputImageSpacingType inputImageSpacing;
  InputImagePointType   inputImageOrigin;

  while ( !inputIt.IsAtEnd() )
    {
    outputIt.Set( inputIt.Get() );
    ++inputIt;
    ++outputIt;
    }
}

template< class TInputImage, class TOutputImage >
void AngleOfIncidenceImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Ultrasound origin vector : "
     << static_cast< typename NumericTraits< VectorType >::PrintType >( m_UltrasoundOrigin )
     << std::endl;
}
} // end namespace itk
#endif
