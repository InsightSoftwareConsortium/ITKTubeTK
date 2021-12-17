/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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
#ifndef __tubeSegmentConnectedComponents_hxx
#define __tubeSegmentConnectedComponents_hxx


namespace tube
{

template< class TImage, class TSeedMask >
SegmentConnectedComponents< TImage, TSeedMask >
::SegmentConnectedComponents( void )
{
  m_Filter = FilterType::New();

  m_SeedMask = nullptr;

  m_MinimumVolume = 0;
  m_NumberOfComponents = 0;
}

template< class TImage, class TSeedMask >
void
SegmentConnectedComponents< TImage, TSeedMask >
::Update( void )
{
  Superclass::Update();

  m_Filter->Update();

  m_NumberOfComponents = m_Filter->GetObjectCount();

  typename ImageType::Pointer curConnComp = m_Filter->GetOutput();

  itk::ImageRegionConstIterator< ImageType > inputIter( m_Filter->GetInput(),
    m_Filter->GetInput()->GetLargestPossibleRegion() );
  inputIter.GoToBegin();
  itk::ImageRegionIterator< ImageType > iter( curConnComp,
    curConnComp->GetLargestPossibleRegion() );
  iter.GoToBegin();

  while( !iter.IsAtEnd() )
    {
    if( inputIter.Get() == 0 )
      {
      iter.Set( 0 );
      }
    else
      {
      unsigned int c = iter.Get();
      iter.Set( c+1 );
      }
    ++inputIter;
    ++iter;
    }

  unsigned int numObjects = m_Filter->GetObjectCount()+1;
  std::vector< bool > cSize( numObjects, true );
  if( m_MinimumVolume > 0 )
    {
    // compute the size ( number of pixels ) of each connected component
    iter.GoToBegin();
    std::vector< unsigned int > cPixelCount( numObjects, 0 );
    while( !iter.IsAtEnd() )
      {
      unsigned int c = iter.Get();
      if( c > 0 && c < numObjects )
        {
        cPixelCount[c]++;
        }
      ++iter;
      }

    // compute voxelVolume
    double voxelVolume = 1;
    for( unsigned int i = 0; i < ImageDimension; i++ )
      {
      voxelVolume *= m_Filter->GetInput()->GetSpacing()[i];
      }
    double minimumVoxelCount = m_MinimumVolume / voxelVolume;
    for( unsigned int c = 0; c<numObjects; ++c )
      {
      if( cPixelCount[c] < minimumVoxelCount )
        {
        cSize[c] = false;
        --m_NumberOfComponents;
        }
      }

    // drop connected components of size ( physp ) below a user-specified cutoff
    iter.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      unsigned int c = iter.Get();
      if( c > 0 && c < numObjects && cSize[c] == false )
        {
        iter.Set( 0 );
        }
      ++iter;
      }
    }

  //
  if( m_SeedMask.IsNotNull() )
    {
    itk::ImageRegionIterator< SeedMaskType > seedIter( m_SeedMask,
      m_SeedMask->GetLargestPossibleRegion() );
    seedIter.GoToBegin();

    iter.GoToBegin();

    std::vector< bool > cSeeded( numObjects, false );
    while( !iter.IsAtEnd() )
      {
      if( seedIter.Get() != 0 )
        {
        unsigned int c = iter.Get();
        if( c > 0 && c < numObjects && !cSeeded[c] )
          {
          if( cSize[c] == true )
            {
            cSeeded[c] = true;
            }
          }
        }
      ++iter;
      ++seedIter;
      }
    m_NumberOfComponents = 0;
    for( unsigned int c = 0; c<numObjects; ++c )
      {
      if( cSeeded[c] )
        {
        ++m_NumberOfComponents;
        }
      }

    iter.GoToBegin();
    while( !iter.IsAtEnd() )
      {
      unsigned int c = iter.Get();
      if( c > 0 && c < numObjects )
        {
        if( !cSeeded[c] )
          {
          iter.Set( 0 );
          }
        }
      ++iter;
      }
    }
}

template< class TImage, class TSeedMask >
void
SegmentConnectedComponents< TImage, TSeedMask >
::PrintSelf( std::ostream & os, itk::Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Filter = " << m_Filter << std::endl;
  os << indent << "MinimumVolume = " << m_MinimumVolume << std::endl;
  os << indent << "SeedMask = " << m_SeedMask << std::endl;
}

}

#endif
