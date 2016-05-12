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
#ifndef __itktubeMergeAdjacentImagesFilter_hxx
#define __itktubeMergeAdjacentImagesFilter_hxx

#include "itktubeMergeAdjacentImagesFilter.h"

namespace itk
{

namespace tube
{

template< class TImage>
MergeAdjacentImagesFilter<TImage>
::MergeAdjacentImagesFilter( void )
{
  this->SetNumberOfRequiredInputs(2);

  m_Background = 0;
  m_MaskZero = false;
  m_Iterations = 300;
}

template< class TImage>
void
MergeAdjacentImagesFilter<TImage>
::SetInput1(const TImage* image)
{
  this->SetInput(0, image);
}

template< class TImage>
void
MergeAdjacentImagesFilter<TImage>
::SetInput2(const TImage* image)
{
  this->SetInput(1, image);
}

template< class TImage>
void
MergeAdjacentImagesFilter<TImage>
::GenerateData()
{
  typename TImage::ConstPointer input1 = this->GetInput( 0 );
  typename TImage::ConstPointer input2 = this->GetInput( 1 );

  typename TImage::Pointer output = this->GetOutput();
  output->SetRegions(input->GetLargestPossibleRegion());
  output->Allocate();

}


} // End namespace tube

} // End namespace itk

#endif
