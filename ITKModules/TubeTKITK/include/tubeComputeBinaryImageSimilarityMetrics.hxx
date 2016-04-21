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
#ifndef __tubeComputeBinaryImageSimilarityMetrics_hxx
#define __tubeComputeBinaryImageSimilarityMetrics_hxx

#include "tubeComputeBinaryImageSimilarityMetrics.h"


namespace tube
{
template< typename TInputImage >
ComputeBinaryImageSimilarityMetrics< TInputImage >
::ComputeBinaryImageSimilarityMetrics( void )
{
  m_MetricFilter = MetricFilterType::New();
}

template< typename TInputImage >
void
ComputeBinaryImageSimilarityMetrics< TInputImage >
::SetSourceImage( const TInputImage *image )
{
  m_MetricFilter->SetSourceImage( image );
}

template< typename TInputImage >
void
ComputeBinaryImageSimilarityMetrics< TInputImage >
::SetTargetImage( const TInputImage *image )
{
  m_MetricFilter->SetTargetImage( image );
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetTotalOverlap( void ) const
{
  return m_MetricFilter->GetTotalOverlap();
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetUnionOverlap( void ) const
{
  return m_MetricFilter->GetUnionOverlap();
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetMeanOverlap( void ) const
{
  return m_MetricFilter->GetMeanOverlap();
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetSimilarity( void ) const
{
  return m_MetricFilter->GetSimilarity();
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetFalseNegativeError( void ) const
{
  return m_MetricFilter->GetFalseNegativeError();
}

template< typename TInputImage >
float
ComputeBinaryImageSimilarityMetrics< TInputImage >
::GetFalsePositiveError( void ) const
{
  return m_MetricFilter->GetFalsePositiveError();
}

template< typename TInputImage >
void
ComputeBinaryImageSimilarityMetrics< TInputImage >
::Update()
{
  m_MetricFilter->Update();
}

template< typename TInputImage >
void
ComputeBinaryImageSimilarityMetrics< TInputImage >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  m_MetricFilter->PrintSelf( os, indent);
}

} // end namespace tube


#endif
