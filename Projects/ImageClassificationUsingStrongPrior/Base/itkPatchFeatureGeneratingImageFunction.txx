/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkPatchFeatureGeneratingImageFunction_txx
#define __itkPatchFeatureGeneratingImageFunction_txx

#include "itkPatchFeatureGeneratingImageFunction.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

namespace itk
{

/**
 * Set the input Image
 */
template < class TInputImage, class TCoordRep >
PatchFeatureGeneratingImageFunction< TInputImage, TCoordRep >
::PatchFeatureGeneratingImageFunction()
{
  this->m_Features.clear();
  this->m_Features.push_back( "PatchMean-PriorMean" );
  this->m_Features.push_back( "PatchStdDev-PriorStdDev" );
  this->m_Features.push_back( "PatchCrossCorrelation" );
  this->m_Features.push_back( "PatchDifferenceMin" );
  this->m_Features.push_back( "PatchDifferenceP25" );
  this->m_Features.push_back( "PatchDifferenceP50" );
  this->m_Features.push_back( "PatchDifferenceP75" );
  this->m_Features.push_back( "PatchDifferenceP95" );
  this->m_Features.push_back( "PatchDifferenceMax" );
  this->m_Features.push_back( "PatchDifferenceSum" );
  this->m_Features.push_back( "PatchDifferenceL2Norm" );
}

template < class TInputImage, class TCoordRep >
typename PatchFeatureGeneratingImageFunction< TInputImage,
  TCoordRep >::OutputType
PatchFeatureGeneratingImageFunction< TInputImage, TCoordRep >
::EvaluateAtIndex( const IndexType & index ) const
{
  typename NeighborIterType::SizeType patchSize  = {{m_PatchWidth,
                                                     m_PatchWidth}};
  NeighborIterType imageItr( patchSize, this->GetInputImage(),
                            this->GetInputImage()->GetBufferedRegion() );
  NeighborIterType priorItr( patchSize, m_Prior,
                            m_Prior->GetBufferedRegion() );

  imageItr.SetLocation( index );
  priorItr.SetLocation( index );

  typename NeighborIterType::NeighborhoodType imageN =
    imageItr.GetNeighborhood();
  typename NeighborIterType::NeighborhoodType priorN =
    priorItr.GetNeighborhood();

  assert( imageN.Size() == priorN.Size() );

  std::vector<float> normalizedImagePatch( imageN.Size() );
  std::vector<float> normalizedPriorPatch( imageN.Size() );

  // compute the mean intensity of the two patches
  float imageMean = std::accumulate( imageN.Begin(), imageN.End(), 0.0 )
    / imageN.Size();
  float priorMean = std::accumulate( priorN.Begin(), priorN.End(), 0.0 )
    / imageN.Size();

  std::transform( imageN.Begin(), imageN.End(),
                 normalizedImagePatch.begin(),
                 std::bind2nd( std::minus<float>(), imageMean ) );

  std::transform( priorN.Begin(), priorN.End(),
                 normalizedPriorPatch.begin(),
                 std::bind2nd( std::minus<float>(), priorMean ) );

  float imageStdDev = std::sqrt( std::inner_product(
    normalizedImagePatch.begin(), normalizedImagePatch.end(),
    normalizedImagePatch.begin(), 0.0 ) / ( imageN.Size() - 1 ) );

  float priorStdDev = std::sqrt( std::inner_product(
    normalizedPriorPatch.begin(), normalizedPriorPatch.end(),
    normalizedPriorPatch.begin(), 0.0 ) / ( priorN.Size() - 1 ) );

  float imageNorm = std::sqrt( std::inner_product(
    normalizedImagePatch.begin(), normalizedImagePatch.end(),
    normalizedImagePatch.begin(), 0.0 ) );

  float priorNorm = std::sqrt( std::inner_product(
    normalizedPriorPatch.begin(), normalizedPriorPatch.end(),
    normalizedPriorPatch.begin(), 0.0 ) );

  if( imageNorm > 0.0 )
    {
    std::transform( normalizedImagePatch.begin(),
      normalizedImagePatch.end(), normalizedImagePatch.begin(),
      std::bind2nd( std::divides<float>(), imageNorm ) );
    }

  if( priorNorm > 0.0 )
    {
    std::transform( normalizedPriorPatch.begin(),
      normalizedPriorPatch.end(), normalizedPriorPatch.begin(),
      std::bind2nd( std::divides<float>(), priorNorm ) );
    }

  float crossCorrelation = std::inner_product( normalizedImagePatch.begin(),
    normalizedImagePatch.end(), normalizedPriorPatch.begin(), 0.0 );

  std::vector<float> differencePatch( imageN.Size() );

  std::transform( imageN.Begin(), imageN.End(), priorN.Begin(),
    differencePatch.begin(), std::minus<float>() );

  float norm = std::inner_product( differencePatch.begin(),
    differencePatch.end(), differencePatch.begin(), 0.0 );

  float total = std::accumulate( differencePatch.begin(),
    differencePatch.end(), 0.0 );

  float maxdiff = *std::max_element( differencePatch.begin(),
    differencePatch.end() );

  float mindiff = *std::min_element( differencePatch.begin(),
    differencePatch.end() );

  size_t q1 = static_cast<size_t>( .25*differencePatch.size() );
  size_t q2 = static_cast<size_t>( .50*differencePatch.size() );
  size_t q3 = static_cast<size_t>( .75*differencePatch.size() );
  size_t q95 = static_cast<size_t>( .95*differencePatch.size() );

  std::sort( differencePatch.begin(), differencePatch.end() );

  float q1val = differencePatch[q1];
  float q2val = differencePatch[q2];
  float q3val = differencePatch[q3];
  float q95val = differencePatch[q95];

  assert( this->m_Features.size() == 11 );
  OutputType out( this->m_Features.size() );
  out[0] = imageMean-priorMean;
  out[1] = imageStdDev-priorStdDev;
  out[2] = crossCorrelation;
  out[3] = mindiff;
  out[4] = q1val;
  out[5] = q2val;
  out[6] = q3val;
  out[7] = q95val;
  out[8] = maxdiff;
  out[9] = total;
  out[10] = norm;

  return out;
}

}

#endif
