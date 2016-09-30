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
#ifndef __tubeShrinkWithBlendingImage_hxx
#define __tubeShrinkWithBlendingImage_hxx

#include "tubeShrinkWithBlendingImage.h"

namespace tube {

/**
 *
 */
template< class TInputImage, class TOutputImage >
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::ShrinkWithBlendingImage( void )
{
  m_Filter = FilterType::New();
}

template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::SetShrinkFactor(unsigned int i, unsigned int factor)
{
  if( m_Filter->GetShrinkFactor(i) != factor )
    {
    m_Filter->SetShrinkFactor(i,factor);
    this->Modified();
    }
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
unsigned int
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::GetShrinkFactor(unsigned int i)
{
  return m_Filter->GetShrinkFactor(i);
}

/**
 *
 */
template< class TInputImage, class TOutputImage >
void
ShrinkWithBlendingImage< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, itk::Indent  indent ) const
{
  os << indent << "ShrinkFactors:"
     << m_Filter->GetShrinkFactors() << std::endl;

  os << indent << "NewSize:"
     << m_Filter->GetNewSize() << std::endl;

  os << indent << "Overlap:" << m_Filter->GetOverlap() << std::endl;

  os << indent << "BlendWithMean:"
     << m_Filter->GetBlendWithMean() << std::endl;

  os << indent << "BlendWithMax:"
     << m_Filter->GetBlendWithMax() << std::endl;

  os << indent << "BlendWithGaussianWeighting:"
     << m_Filter->GetBlendWithGaussianWeighting() << std::endl;

  os << indent << "UseLog:"
     << m_Filter->GetUseLog() << std::endl;

  if( m_Filter->GetInputMipPointImage() != ITK_NULLPTR )
    {
    os << indent << "Input MIP Point Image: "
       << m_Filter->GetInputMipPointImage() << std::endl;
    }
  else
    {
    os << indent << "Input MIP Point Image: NULL" << std::endl;
    }

  if( m_Filter->GetOutputMipPointImage() != ITK_NULLPTR )
    {
    os << indent << "Output MIP Point Image: "
       << m_Filter->GetOutputMipPointImage() << std::endl;
    }
  else
    {
    os << indent << "Output MIP Point Image: NULL" << std::endl;
    }
}


} // end namespace tube

#endif
