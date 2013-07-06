/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeLabelOverlapMeasuresImageFilter_hxx
#define __itktubeLabelOverlapMeasuresImageFilter_hxx

#include "itktubeLabelOverlapMeasuresImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkProgressReporter.h>

namespace itk
{

namespace tube
{

//NOTE: This class needs a mutex for gnu 2.95
#if defined(__GNUC__) && (__GNUC__ <= 2)
/** Used for mutex locking */
#define LOCK_HASHMAP this->m_Mutex.Lock()
#define UNLOCK_HASHMAP this->m_Mutex.Unlock()
#else
#define LOCK_HASHMAP
#define UNLOCK_HASHMAP
#endif

template< class TLabelImage >
LabelOverlapMeasuresImageFilter< TLabelImage >
::LabelOverlapMeasuresImageFilter( void )
{
  // this filter requires two input images
  this->SetNumberOfRequiredInputs( 2 );
}

template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::GenerateInputRequestedRegion( void )
{
  Superclass::GenerateInputRequestedRegion();
  if( this->GetSourceImage() )
    {
    typename LabelImageType::Pointer source = const_cast
      <LabelImageType *>( this->GetSourceImage() );
    source->SetRequestedRegionToLargestPossibleRegion();
    }
  if( this->GetTargetImage() )
    {
    typename LabelImageType::Pointer target = const_cast
      <LabelImageType *>( this->GetTargetImage() );
    target->SetRequestedRegionToLargestPossibleRegion();
    }
}

template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::EnlargeOutputRequestedRegion( DataObject *data )
{
  Superclass::EnlargeOutputRequestedRegion( data );
  data->SetRequestedRegionToLargestPossibleRegion();
}


template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::BeforeThreadedGenerateData( void )
{
  int numberOfThreads = this->GetNumberOfThreads();

  // Resize the thread temporaries
  this->m_LabelSetMeasuresPerThread.resize( numberOfThreads );

  // Initialize the temporaries
  for( int n = 0; n < numberOfThreads; n++ )
    {
    this->m_LabelSetMeasuresPerThread[n].clear();
    }

  // Initialize the final map
  this->m_LabelSetMeasures.clear();
}

template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::AfterThreadedGenerateData( void )
{
  // Run through the map for each thread and accumulate the set measures.
  for( ThreadIdType n = 0; n < this->GetNumberOfThreads(); n++ )
    {
    // iterate over the map for this thread
    for( MapConstIterator threadIt =
        this->m_LabelSetMeasuresPerThread[n].begin();
      threadIt != this->m_LabelSetMeasuresPerThread[n].end();
      ++threadIt )
      {
      // does this label exist in the cumulative stucture yet?
      MapIterator mapIt = this->m_LabelSetMeasures.find(
        ( *threadIt ).first );
      if( mapIt == this->m_LabelSetMeasures.end() )
        {
        // create a new entry
        typedef typename MapType::value_type MapValueType;
        mapIt = this->m_LabelSetMeasures.insert( MapValueType(
          (*threadIt).first, LabelSetMeasures() ) ).first;
        }

      // accumulate the information from this thread
      (*mapIt).second.m_Source += (*threadIt).second.m_Source;
      (*mapIt).second.m_Target += (*threadIt).second.m_Target;
      (*mapIt).second.m_Union += (*threadIt).second.m_Union;
      (*mapIt).second.m_Intersection +=
        (*threadIt).second.m_Intersection;
      (*mapIt).second.m_SourceComplement +=
        (*threadIt).second.m_SourceComplement;
      (*mapIt).second.m_TargetComplement +=
        (*threadIt).second.m_TargetComplement;
      } // end of thread map iterator loop
    } // end of thread loop
}

template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::ThreadedGenerateData( const RegionType& outputRegionForThread,
  ThreadIdType threadId )
{
  ImageRegionConstIterator<LabelImageType > ItS( this->GetSourceImage(),
    outputRegionForThread );
  ImageRegionConstIterator<LabelImageType > ItT( this->GetTargetImage(),
    outputRegionForThread );

  // support progress methods/callbacks
  ProgressReporter progress( this, threadId,
    2*outputRegionForThread.GetNumberOfPixels() );

  for( ItS.GoToBegin(), ItT.GoToBegin(); !ItS.IsAtEnd(); ++ItS, ++ItT )
    {
    LabelType sourceLabel = ItS.Get();
    LabelType targetLabel = ItT.Get();

    // is the label already in this thread?
    MapIterator mapItS =
      this->m_LabelSetMeasuresPerThread[threadId].find( sourceLabel );
    MapIterator mapItT =
      this->m_LabelSetMeasuresPerThread[threadId].find( targetLabel );

    if( mapItS == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItS = this->m_LabelSetMeasuresPerThread[threadId].insert(
        MapValueType( sourceLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    if( mapItT == this->m_LabelSetMeasuresPerThread[threadId].end() )
      {
      // create a new label set measures object
      typedef typename MapType::value_type MapValueType;
      LOCK_HASHMAP;
      mapItT = this->m_LabelSetMeasuresPerThread[threadId].insert(
        MapValueType( targetLabel, LabelSetMeasures() ) ).first;
      UNLOCK_HASHMAP;
      }

    (*mapItS).second.m_Source++;
    (*mapItT).second.m_Target++;

    if( sourceLabel == targetLabel )
      {
      (*mapItS).second.m_Intersection++;
      (*mapItS).second.m_Union++;
      }
    else
      {
      (*mapItS).second.m_Union++;
      (*mapItT).second.m_Union++;

      (*mapItS).second.m_SourceComplement++;
      (*mapItT).second.m_TargetComplement++;
      }

    progress.CompletedPixel();
    }
}

/**
 *  measures
 */
template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetTotalOverlap( void )
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType >::Zero )
      {
      continue;
      }
    numerator += static_cast< RealType >( (*mapIt).second.m_Intersection );
    denominator += static_cast< RealType >( (*mapIt).second.m_Target );
    }
  return ( numerator / denominator );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetTargetOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( << "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast< RealType >( (*mapIt).second.m_Intersection ) /
    static_cast< RealType >( (*mapIt).second.m_Target );
  return value;
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetUnionOverlap( void )
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType >::Zero )
      {
      continue;
      }
    numerator += static_cast< RealType >( (*mapIt).second.m_Intersection );
    denominator += static_cast< RealType >( (*mapIt).second.m_Union );
    }
  return ( numerator / denominator );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetUnionOverlap( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( << "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast< RealType >( (*mapIt).second.m_Intersection ) /
    static_cast< RealType >( (*mapIt).second.m_Union );
  return value;
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetMeanOverlap( void )
{
  RealType uo = this->GetUnionOverlap();
  return ( 2.0 * uo / ( 1.0 + uo ) );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetMeanOverlap( LabelType label )
{
  RealType uo = this->GetUnionOverlap( label );
  return ( 2.0 * uo / ( 1.0 + uo ) );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetSimilarity( void )
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType >::Zero )
      {
      continue;
      }
    numerator += ( ( static_cast< RealType >( (*mapIt).second.m_Source ) -
      static_cast< RealType >( (*mapIt).second.m_Target ) ) );
    denominator += ( ( static_cast< RealType >( (*mapIt).second.m_Source ) +
      static_cast< RealType >( (*mapIt).second.m_Target ) ) );
    }
  return ( 2.0 * numerator / denominator );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetSimilarity( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( << "Label " << label << " not found." );
    return 0.0;
    }
  RealType value = 2.0 *
    ( static_cast< RealType >( (*mapIt).second.m_Source ) -
      static_cast< RealType >( (*mapIt).second.m_Target ) ) /
    ( static_cast< RealType >( (*mapIt).second.m_Source ) +
      static_cast< RealType >( (*mapIt).second.m_Target ) );
  return value;
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetFalseNegativeError( void )
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType >::Zero )
      {
      continue;
      }
    numerator += static_cast< RealType >( (*mapIt).second.m_TargetComplement );
    denominator += static_cast< RealType >( (*mapIt).second.m_Target );
    }
  return ( numerator / denominator );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetFalseNegativeError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( << "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast< RealType >( (*mapIt).second.m_TargetComplement ) /
    static_cast< RealType >( (*mapIt).second.m_Target );
  return value;
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetFalsePositiveError( void )
{
  RealType numerator = 0.0;
  RealType denominator = 0.0;
  for( MapIterator mapIt = this->m_LabelSetMeasures.begin();
    mapIt != this->m_LabelSetMeasures.end(); ++mapIt )
    {
    // Do not include the background in the final value.
    if( (*mapIt).first == NumericTraits<LabelType >::Zero )
      {
      continue;
      }
    numerator += static_cast< RealType >( (*mapIt).second.m_SourceComplement );
    denominator += static_cast< RealType >( (*mapIt).second.m_Source );
    }
  return ( numerator / denominator );
}

template< class TLabelImage >
typename LabelOverlapMeasuresImageFilter< TLabelImage >::RealType
LabelOverlapMeasuresImageFilter< TLabelImage >
::GetFalsePositiveError( LabelType label )
{
  MapIterator mapIt = this->m_LabelSetMeasures.find( label );
  if( mapIt == this->m_LabelSetMeasures.end() )
    {
    itkWarningMacro( << "Label " << label << " not found." );
    return 0.0;
    }
  RealType value =
    static_cast< RealType >( (*mapIt).second.m_SourceComplement ) /
    static_cast< RealType >( (*mapIt).second.m_Source );
  return value;
}

template< class TLabelImage >
void
LabelOverlapMeasuresImageFilter< TLabelImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeLabelOverlapMeasuresImageFilter_hxx)
