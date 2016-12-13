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
#ifndef __tubeComputeSegmentTubesParameters_hxx
#define __tubeComputeSegmentTubesParameters_hxx

#include "tubeComputeSegmentTubesParameters.h"


namespace tube
{
template< class TPixel, unsigned int VDimension >
ComputeSegmentTubesParameters< TPixel, VDimension >
::ComputeSegmentTubesParameters( void )
{
  m_Filter = FilterType::New();
}

template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetSeedData( void )
{
  return m_Filter->GetSeedData();
}

template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetTubeData( void )
{
  return m_Filter->GetTubeData();
}

template< class TPixel, unsigned int VDimension >
std::vector< vnl_vector< double > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetBkgData( void )
{
  return m_Filter->GetBkgData();
}

template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetSeedDataIndexList( void )
{
  return m_Filter->GetSeedDataIndexList();
}

template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetTubeDataIndexList( void )
{
  return m_Filter->GetTubeDataIndexList();
}

template< class TPixel, unsigned int VDimension >
std::vector< itk::ContinuousIndex< double, VDimension > >
ComputeSegmentTubesParameters< TPixel, VDimension >
::GetBkgDataIndexList( void )
{
  return m_Filter->GetBkgDataIndexList();
}

template< class TPixel, unsigned int VDimension >
void
ComputeSegmentTubesParameters< TPixel, VDimension >
::SetParameterFile( std::string s )
{
  return m_Filter->SetParameterFile( s );
}

template< class TPixel, unsigned int VDimension >
void
ComputeSegmentTubesParameters< TPixel, VDimension >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << m_Filter << std::endl;
}

} // end namespace tube

#endif
