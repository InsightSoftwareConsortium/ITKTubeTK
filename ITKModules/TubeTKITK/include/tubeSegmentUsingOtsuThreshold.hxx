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
#ifndef __tubeSegmentUsingOtsuThreshold_hxx
#define __tubeSegmentUsingOtsuThreshold_hxx

#include "tubeSegmentUsingOtsuThreshold.h"

namespace tube
{

template< class TInputPixel, unsigned int Dimension, class TMaskPixel >
SegmentUsingOtsuThreshold< TInputPixel, Dimension, TMaskPixel >
::SegmentUsingOtsuThreshold( void )
{
  m_Filter = FilterType::New();
  m_Filter->SetMaskOutput( false );
  m_Filter->SetOutsideValue( 0 );
  m_Filter->SetInsideValue( 1 );
}

template< class TInputPixel, unsigned int Dimension, class TMaskPixel >
void
SegmentUsingOtsuThreshold< TInputPixel, Dimension, TMaskPixel >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "Inside Value: " << m_Filter->GetInsideValue() << std::endl;
  os << "Outside Value: " << m_Filter->GetOutsideValue() << std::endl;
  os << "Mask Value: " << m_Filter->GetMaskValue() << std::endl;
}

}

#endif
