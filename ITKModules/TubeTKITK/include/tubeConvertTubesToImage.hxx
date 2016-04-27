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

template< class TPixel, unsigned int Dimension >
ConvertTubesToImage< TPixel, Dimension >
::ConvertTubesToImage( void )
{
  m_TubesToImageFilter = TubesToImageFilterType::New();
}

template< class TPixel, unsigned int Dimension >
void
ConvertTubesToImage< TPixel, Dimension >
::SetTemplateImage( TemplateImageType::Pointer pTemplateImage )
{
  TubetoImageFilterType::SizeType size;

  double spacing[Dimension];
  for(int i = 0; i < Dimension; i++ )
    {
    size[i] = pTemplateImage->GetLargestPossibleRegion().GetSize()[i];
    spacing[i] = pTemplateImage->GetSpacing()[i];
    }

  m_TubesToImageFilter->SetSize( size );
  m_TubesToImageFilter->SetSpacing( spacing );
}

template< class TPixel, unsigned int Dimension >
void
ConvertTubesToImage< TPixel, Dimension >
::SetInput( TubesType::Pointer pTubes )
{
  m_TubesToImageFilter->SetInput( pTubes );
}

template< class TPixel, unsigned int Dimension >
void
ConvertTubesToImage< TPixel, Dimension >
::Update()
{
  m_TubesToImageFilter->Update();
}

template< class TPixel, unsigned int Dimension >
ConvertTubesToImage< TPixel, Dimension >::OutputImageType::Pointer
ConvertTubesToImage< TPixel, Dimension >
::GetOutput()
{
  return m_TubesToImageFilter->GetOutput();
}

template< class TPixel, unsigned int Dimension >
void
ConvertTubesToImage< TPixel, Dimension >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  m_TubesToImageFilter->PrintSelf( os, indent );
}

}

#endif
