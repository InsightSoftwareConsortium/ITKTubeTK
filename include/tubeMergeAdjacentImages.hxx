/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeMergeAdjacentImages_hxx
#define __tubeMergeAdjacentImages_hxx


namespace tube
{

template< class TImage >
MergeAdjacentImages< TImage >
::MergeAdjacentImages()
{
  m_Filter = FilterType::New();
}

template< class TImage >
void
MergeAdjacentImages< TImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << "Background: " << this->GetBackground() << std::endl;
  os << "MaskZero: " << this->GetMaskZero() << std::endl;
  os << "MaxIterations: " << this->GetMaxIterations() << std::endl;
  os << "ExpectedOffset: " << this->GetExpectedOffset() << std::endl;
  os << "ExpectedRotation: " << this->GetExpectedRotation() << std::endl;
  os << "SamplingRatio: " << this->GetSamplingRatio() << std::endl;
  os << "BlendUsingAverage: " << this->GetBlendUsingAverage() << std::endl;
  os << "UseFastBlending: " << this->GetUseFastBlending() << std::endl;

}

} // End namespace tube

#endif
