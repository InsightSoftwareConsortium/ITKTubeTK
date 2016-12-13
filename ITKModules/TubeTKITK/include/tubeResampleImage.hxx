/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 ( the "License" );
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
#ifndef __tubeResampleImage_hxx
#define __tubeResampleImage_hxx

#include "tubeResampleImage.h"

namespace tube
{
template< class TPixel, unsigned int VDimension >
ResampleImage< TPixel, VDimension >
::ResampleImage( void )
{
  m_Filter = FilterType::New();
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::SetTransform( TransformType* t )
{
  m_Filter->SetTransform( t );
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::SetSpacing( std::vector<double> s )
{
  m_Filter->SetSpacing( s );
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::SetOrigin( std::vector<double> o )
{
  m_Filter->SetOrigin( o );
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::SetIndex( std::vector<int> i )
{
  m_Filter->SetIndex( i );
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::SetResampleFactor( std::vector<double> rf )
{
  m_Filter->SetResampleFactor( rf );
}

template< class TPixel, unsigned int VDimension >
void
ResampleImage< TPixel, VDimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << m_Filter << std::endl;
}

} // end namespace tube

#endif
