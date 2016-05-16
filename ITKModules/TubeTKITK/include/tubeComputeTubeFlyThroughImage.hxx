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
#ifndef __tubeComputeTubeFlyThroughImage_hxx
#define __tubeComputeTubeFlyThroughImage_hxx

#include "tubeComputeTubeFlyThroughImage.h"

namespace tube
{

template< class TPixel, unsigned int Dimension >
ComputeTubeFlyThroughImage< TPixel, Dimension >
::ComputeTubeFlyThroughImage( void )
{
  m_Filter = FilterType::New();
}

template< class TPixel, unsigned int Dimension >
void
ComputeTubeFlyThroughImage< TPixel, Dimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << "TubeId: " << this->GetTubeId() << std::endl;
}

}

#endif
