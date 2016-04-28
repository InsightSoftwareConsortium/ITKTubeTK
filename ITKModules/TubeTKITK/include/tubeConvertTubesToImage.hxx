/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeConvertTubesToImage_hxx
#define __tubeConvertTubesToImage_hxx

#include "tubeConvertTubesToImage.h"

namespace tube
{

template< unsigned int Dimension, class TOutputPixel >
ConvertTubesToImage< Dimension, TOutputPixel >
::ConvertTubesToImage( void )
{
  m_TubesToImageFilter = TubesToImageFilterType::New();
  m_TubesToImageFilter->SetBuildRadiusImage( false );
  m_TubesToImageFilter->SetBuildTangentImage( false );
}

template< unsigned int Dimension, class TOutputPixel >
void
ConvertTubesToImage< Dimension, TOutputPixel >
::SetTemplateImage( typename ConvertTubesToImage< Dimension,
  TOutputPixel >::OutputImageType::Pointer pTemplateImage )
{
  typename TubesToImageFilterType::SizeType size;

  double spacing[Dimension];
  for(int i = 0; i < Dimension; i++ )
    {
    size[i] = pTemplateImage->GetLargestPossibleRegion().GetSize()[i];
    spacing[i] = pTemplateImage->GetSpacing()[i];
    }

  m_TubesToImageFilter->SetSize( size );
  m_TubesToImageFilter->SetSpacing( spacing );
}

template< unsigned int Dimension, class TOutputPixel >
void
ConvertTubesToImage< Dimension, TOutputPixel >
::SetUseRadius( bool useRadius )
{
  m_TubesToImageFilter->SetUseRadius( useRadius );
}

template< unsigned int Dimension, class TOutputPixel >
void
ConvertTubesToImage< Dimension, TOutputPixel >
::SetInput( typename ConvertTubesToImage< Dimension,
  TOutputPixel >::TubesType::Pointer pTubes )
{
  m_TubesToImageFilter->SetInput( pTubes );
}

template< unsigned int Dimension, class TOutputPixel >
void
ConvertTubesToImage< Dimension, TOutputPixel >
::Update()
{
  m_TubesToImageFilter->Update();
}

template< unsigned int Dimension, class TOutputPixel >
typename ConvertTubesToImage< Dimension,
  TOutputPixel >::OutputImageType::Pointer
ConvertTubesToImage< Dimension, TOutputPixel >
::GetOutput()
{
  return m_TubesToImageFilter->GetOutput();
}

template< unsigned int Dimension, class TOutputPixel >
void
ConvertTubesToImage< Dimension, TOutputPixel >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  m_TubesToImageFilter->PrintSelf( os, indent );
}

}

#endif
