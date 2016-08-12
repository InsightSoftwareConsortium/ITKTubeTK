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

#ifndef __itktubeComputeTubeMeasuresFilter_hxx
#define __itktubeComputeTubeMeasuresFilter_hxx

#include "itktubeComputeTubeMeasuresFilter.h"

namespace itk
{

namespace tube
{

/**
 * Constructor
 */
template< class TPixel, unsigned int Dimension >
ComputeTubeMeasuresFilter< TPixel, Dimension >
::ComputeTubeMeasuresFilter( void )
{
  m_Scale = 0;
  m_Ridgeness = NULL;
  m_Roundness = NULL;
  m_Curvature = NULL;
  m_Levelness = NULL;
}

template< class TPixel, unsigned int Dimension >
void
ComputeTubeMeasuresFilter< TPixel, Dimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  SuperClass::PrintSelf( os, indent );
  os << "Scale: " << m_Scale << std::endl;
}

template< class TPixel, unsigned int Dimension >
void
ComputeTubeMeasuresFilter< TPixel, Dimension >
::GenerateData( void )
{
  itkDebugMacro( << "ComputeTubeMeasuresFilter::Update() called." );

  m_InputImage = this->GetInput();
  if ( m_Scale > 0 )
    {
    typename RescaleFilterType::Pointer rescaleFilter
      = RescaleFilterType::New();
    rescaleFilter->SetInput( m_InputImage );
    rescaleFilter->SetOutputMinimum( 0 );
    rescaleFilter->SetOutputMaximum( 1 );

    typename RidgeFilterType::Pointer filter;
    filter = RidgeFilterType::New();
    filter->SetInput( rescaleFilter->GetOutput() );
    filter->SetScale( m_Scale );
    filter->Update();

    m_Ridgeness = filter->GetRidgeness();
    m_Roundness = filter->GetRoundness();
    m_Curvature = filter->GetCurvature();
    m_Levelness = filter->GetLevelness();
    }

  itkDebugMacro( << "ComputeTubeMeasuresFilter::Update() finished." );
}

} // End namespace tube

} // End namespace itk

#endif
