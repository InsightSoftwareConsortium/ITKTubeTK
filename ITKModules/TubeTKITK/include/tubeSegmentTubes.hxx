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
#ifndef __tubeSegmentTubes_hxx
#define __tubeSegmentTubes_hxx

#include "tubeSegmentTubes.h"


namespace tube
{
template< class TInputImage >
SegmentTubes< TInputImage >
::SegmentTubes( void )
{
  m_Filter = FilterType::New();
}

template< class TInputImage >
bool
SegmentTubes< TInputImage >
::AddTube( TubeType * tube )
{
  return m_Filter->AddTube( tube );
}

template< class TInputImage >
bool
SegmentTubes< TInputImage >
::DeleteTube( TubeType * tube )
{
  return m_Filter->DeleteTube( tube );
}

template< class TInputImage >
typename itk::tube::TubeExtractor< TInputImage >::TubeType *
SegmentTubes< TInputImage >
::ExtractTube( const ContinuousIndexType & x,
    unsigned int tubeID,
    bool verbose )
{
  return m_Filter->ExtractTube( x, tubeID, verbose );
}

template< class TInputImage >
void
SegmentTubes< TInputImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << m_Filter << std::endl;

}

} // end namespace tube


#endif
