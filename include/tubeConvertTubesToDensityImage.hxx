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
#ifndef __tubeConvertTubesToDensityImage_hxx
#define __tubeConvertTubesToDensityImage_hxx

#include "tubeConvertTubesToDensityImage.h"

namespace tube
{

template< class TOutputPixel, unsigned int Dimension >
ConvertTubesToDensityImage< TOutputPixel, Dimension >
::ConvertTubesToDensityImage( void )
{
  m_Filter = FilterType::New();
}

template< class TOutputPixel, unsigned int Dimension >
void
ConvertTubesToDensityImage< TOutputPixel, Dimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_Spacing: " << m_Filter->GetSpacing() << std::endl;
  os << indent << "m_Size: " << m_Filter->GetSize() << std::endl;
  os << indent << "m_MaxDensityIntensity: " <<
    m_Filter->GetMaxDensityIntensity() << std::endl;
  os << indent << "m_UseSquaredDistance: " <<
    m_Filter->GetUseSquaredDistance() << std::endl;
}

}

#endif
