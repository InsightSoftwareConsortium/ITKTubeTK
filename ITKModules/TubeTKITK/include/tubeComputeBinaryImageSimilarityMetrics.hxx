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
template< class TInputImage >
ComputeBinaryImageSimilarityMetrics< TInputImage >
::ComputeBinaryImageSimilarityMetrics( void )
{
  m_Filter = FilterType::New();
}

template< class TInputImage >
void
ComputeBinaryImageSimilarityMetrics< TInputImage >
::PrintSelf(std::ostream & os, itk::Indent indent) const
{
  Superclass::PrintSelf( os, indent );
  os << "TotalOverlap: " << this->GetTotalOverlap() << std::endl;
  os << "UnionOverlap: " << this->GetUnionOverlap() << std::endl;
  os << "MeanOverlap: " << this->GetMeanOverlap() << std::endl;
  os << "Similarity: " << this->GetSimilarity() << std::endl;
  os << "FalseNegativeError: " << this->GetFalseNegativeError() << std::endl;
  os << "FalsePositiveError: " << this->GetFalsePositiveError() << std::endl;
}

} // end namespace tube


#endif
