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
#ifndef __tubeEnhanceEdgesUsingDiffusion_hxx
#define __tubeEnhanceEdgesUsingDiffusion_hxx

#include "tubeEnhanceEdgesUsingDiffusion.h"

namespace tube
{

template< class TInputImage, class TOutputImage >
EnhanceEdgesUsingDiffusion< TInputImage, TOutputImage >
::EnhanceEdgesUsingDiffusion( void )
{
  m_Filter = FilterType::New();
}

template< class TInputImage, class TOutputImage >
void
EnhanceEdgesUsingDiffusion< TInputImage, TOutputImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Contrast parameter LambdaE: "
    << m_Filter->GetContrastParameterLambdaE() << std::endl;
  os << indent << "Sigma : " << m_Filter->GetSigma() << std::endl;
  os << indent << "SigmaOuter : " << m_Filter->GetSigmaOuter() << std::endl;
  os << indent << "Threshold parameter C "
    << m_Filter->GetThresholdParameterC() << std::endl;
}

}

#endif
