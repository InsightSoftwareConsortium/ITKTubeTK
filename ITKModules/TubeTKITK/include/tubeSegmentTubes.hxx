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
#ifndef __tubeSegmentTubes_hxx
#define __tubeSegmentTubes_hxx

#include "tubeSegmentTubes.h"


namespace tube
{
template< class TInputImage >
SegmentTubes< TInputImage >
::SegmentTubes( void )
{
  m_TubeExtractorFilter = TubeExtractorFilterType::New();
}

template< class TInputImage >
void
SegmentTubes< TInputImage >
::SetTubeMaskImage( typename
  TubeMaskImageType::Pointer & mask )
{
  return m_TubeExtractorFilter->SetTubeMaskImage( mask );
}

template< class TInputImage >
bool
SegmentTubes< TInputImage >
::AddTube( TubeType * tube )
{
  return m_TubeExtractorFilter->AddTube( tube );
}

template< class TInputImage >
void
SegmentTubes< TInputImage >
::SetExtractBoundMin( const IndexType & dataMin )
{
  return m_TubeExtractorFilter->SetExtractBoundMin( dataMin );
}

template< class TInputImage >
void
SegmentTubes< TInputImage >
::SetExtractBoundMax( const IndexType & dataMax )
{
  return m_TubeExtractorFilter->SetExtractBoundMax( dataMax );
}

template< class TInputImage >
typename itk::tube::TubeExtractor< TInputImage >::TubeType::Pointer
SegmentTubes< TInputImage >
::ExtractTube( const ContinuousIndexType & x,
    unsigned int tubeID,
    bool verbose )
{
  return m_TubeExtractorFilter->ExtractTube( x, tubeID, verbose );
}

template< class TInputImage >
typename itk::tube::TubeExtractor< TInputImage >::TubeGroupType::Pointer
SegmentTubes< TInputImage >
::GetTubeGroup( void )
{
  return m_TubeExtractorFilter->GetTubeGroup();
}

template< class TInputImage >
typename itk::tube::RidgeExtractor< TInputImage >::Pointer
SegmentTubes< TInputImage >
::GetRidgeOp( void )
{
  return m_TubeExtractorFilter->GetRidgeOp();
}

template< class TInputImage >
typename itk::tube::RadiusExtractor2< TInputImage >::Pointer
SegmentTubes< TInputImage >
::GetRadiusOp( void )
{
  return m_TubeExtractorFilter->GetRadiusOp();
}

template< class TInputImage >
void
SegmentTubes< TInputImage >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_TubeExtractorFilter->GetInputImage() )
    {
    os << indent << "Input Image = " <<
      m_TubeExtractorFilter->GetInputImage() << std::endl;
    }
  else
    {
    os << indent << "Input Image = NULL" << std::endl;
    }

  if( m_TubeExtractorFilter->GetRadiusInputImage() )
    {
    os << indent << "Radius Input Image = " <<
      m_TubeExtractorFilter->GetRadiusInputImage()
      << std::endl;
    }
  else
    {
    os << indent << "Radius Input Image = NULL" << std::endl;
    }

}

} // end namespace tube


#endif
