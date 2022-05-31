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

#ifndef __itktubeVotingResampleImageFunction_hxx
#define __itktubeVotingResampleImageFunction_hxx


#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhood.h>

#include <vnl/vnl_math.h>

#include <map>

namespace itk
{

namespace tube
{

/**
 * Define the number of neighbors
 */
template< class TInputImage, class TCoordRep >
const unsigned long
VotingResampleImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template< class TInputImage, class TCoordRep >
VotingResampleImageFunction< TInputImage, TCoordRep >
::VotingResampleImageFunction( void )
{
  m_Radius.Fill( 1 );
}


/**
 * PrintSelf
 */
template< class TInputImage, class TCoordRep >
void
VotingResampleImageFunction< TInputImage, TCoordRep >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );
}


/**
 * Evaluate at image index position
 */
template< class TInputImage, class TCoordRep >
typename VotingResampleImageFunction< TInputImage, TCoordRep >
::OutputType
VotingResampleImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index ) const
{
  typedef itk::ConstNeighborhoodIterator< TInputImage >
    NeighborhoodIteratorType;

  NeighborhoodIteratorType it( m_Radius, this->GetInputImage(),
    this->GetInputImage()->GetRequestedRegion() );

  IndexType newIndex;
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    newIndex[i] = ( int )(index[i]+0.5);  // Round to nearest int
    }

  it.SetLocation( newIndex );
  itk::Neighborhood<typename TInputImage::PixelType, ImageDimension> n =
    it.GetNeighborhood();
  std::map<typename TInputImage::PixelType, int> tally;
  for( unsigned int i = 0; i < n.Size(); i++ )
    {
    tally[n[i]] = 0;
    }
  int max = 0;
  OutputType ret = n[0];
  for( unsigned int i = 0; i < n.Size(); i++ )
    {
    if( ++tally[n[i]] > max )
      {
      max = tally[n[i]];
      ret = n[i];
      }
    }
  return ret;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeVotingResampleImageFunction_hxx )
