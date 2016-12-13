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

#ifndef __itktubeMarkDuplicateFramesInvalidImageFilter_hxx
#define __itktubeMarkDuplicateFramesInvalidImageFilter_hxx

#include "itktubeMarkDuplicateFramesInvalidImageFilter.h"

#include <itkImageSliceConstIteratorWithIndex.h>
#include <itkMetaDataObject.h>

namespace itk
{

namespace tube
{

template< typename TAssociate >
void
MarkDuplicateFramesInvalidImageFilterThreader< TAssociate >
::BeforeThreadedExecution( void )
{
  const ThreadIdType numberOfThreads = this->GetNumberOfThreadsUsed();
  m_InvalidFramesPerThread.resize( numberOfThreads );
  for( ThreadIdType ii = 0; ii < numberOfThreads; ++ii )
    {
    m_InvalidFramesPerThread[ii].clear();
    }
}


template< typename TAssociate >
void
MarkDuplicateFramesInvalidImageFilterThreader< TAssociate >
::ThreadedExecution( const DomainType & regionForThread,
  const ThreadIdType threadId )
{
  typedef typename AssociateType::InputImageType       InputImageType;
  typedef typename AssociateType::InputImageRegionType InputImageRegionType;

  const InputImageType * inputImage = this->m_Associate->GetInput();
  const InputImageRegionType largestRegion =
    inputImage->GetLargestPossibleRegion();
  const typename InputImageRegionType::IndexType largestIndex =
    largestRegion.GetIndex();
  const typename InputImageRegionType::SizeType largestSize =
    largestRegion.GetSize();

  InputImageRegionType nextSliceRegion = regionForThread;
  typename InputImageType::IndexType nextSliceIndex =
    nextSliceRegion.GetIndex();
  ++nextSliceIndex[2];
  typename InputImageType::SizeType nextSliceSize =
    nextSliceRegion.GetSize();
  InputImageRegionType currentSliceRegion = regionForThread;
  const typename InputImageType::IndexType currentSliceIndex =
    currentSliceRegion.GetIndex();
  typename InputImageType::SizeType currentSliceSize =
    currentSliceRegion.GetSize();
  if( currentSliceIndex[2] + currentSliceSize[2] ==
        largestIndex[2] + largestSize[2] )
    {
    --nextSliceSize[2];
    --currentSliceSize[2];
    }
  nextSliceRegion.SetIndex( nextSliceIndex );
  nextSliceRegion.SetSize( nextSliceSize );
  currentSliceRegion.SetSize( currentSliceSize );

  typedef ImageSliceConstIteratorWithIndex< InputImageType > IteratorType;
  IteratorType nextSliceIt( inputImage, nextSliceRegion );
  IteratorType currentSliceIt( inputImage, currentSliceRegion );

  nextSliceIt.GoToBegin();
  nextSliceIt.SetFirstDirection( 0 );
  nextSliceIt.SetSecondDirection( 1 );
  currentSliceIt.GoToBegin();
  currentSliceIt.SetFirstDirection( 0 );
  currentSliceIt.SetSecondDirection( 1 );
  const float tolerance = static_cast< float >( this->m_Associate->m_Tolerance );
  const float numberOfSlicePixels =
    static_cast< float >( largestSize[0] * largestSize[1] );
  SizeValueType currentFrame = currentSliceIndex[2];
  while( !currentSliceIt.IsAtEnd() )
    {
    SizeValueType sameCount = 0;
    while( !currentSliceIt.IsAtEndOfSlice() )
      {
      while( !currentSliceIt.IsAtEndOfLine() )
        {
        const float absDifference =
          std::abs( static_cast< float >( currentSliceIt.Get() ) -
            static_cast< float >( nextSliceIt.Get() ) );
        if( absDifference <= tolerance )
          {
          ++sameCount;
          }
        ++currentSliceIt;
        ++nextSliceIt;
        }
      currentSliceIt.NextLine();
      nextSliceIt.NextLine();
      }
    //std::cout << sameCount / numberOfSlicePixels << std::endl;
    if( sameCount / numberOfSlicePixels > this->m_Associate->m_FractionalThreshold )
      {
      m_InvalidFramesPerThread[threadId].push_back( currentFrame );
      }
    currentSliceIt.NextSlice();
    nextSliceIt.NextSlice();
    ++currentFrame;
    }
}


template< typename TAssociate >
void
MarkDuplicateFramesInvalidImageFilterThreader< TAssociate >
::AfterThreadedExecution( void )
{
  MetaDataDictionary newDictionary =
    *( this->m_Associate->GetInputMetaDataDictionary() );
  const ThreadIdType numberOfThreads = this->GetNumberOfThreadsUsed();
  std::ostringstream keyPrefix;
  for( ThreadIdType ii = 0; ii < numberOfThreads; ++ii )
    {
    for( typename InvalidFramesType::const_iterator invalidIt =
      m_InvalidFramesPerThread[ii].begin();
      invalidIt != m_InvalidFramesPerThread[ii].end();
      ++invalidIt )
      {
      keyPrefix.str( "" );
      keyPrefix << "Seq_Frame";
      keyPrefix.fill( '0' );
      keyPrefix.width( 4 );
      keyPrefix << *invalidIt;
      EncapsulateMetaData< std::string >( newDictionary,
        keyPrefix.str() + "_ProbeToTrackerTransformStatus",
        std::string( "INVALID" ) );
      }
    }
  DataObject * outputDataObject =
    this->m_Associate->GetOutput( "MetaDataDictionaryObject" );
  outputDataObject->SetMetaDataDictionary( newDictionary );
}


template< typename TInputImage >
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::MarkDuplicateFramesInvalidImageFilter( void ):
  m_Tolerance( 2 ),
  m_FractionalThreshold( 0.7 ),
  m_InputMetaDataDictionary( NULL )
{
  m_Threader = ThreaderType::New();

  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );

  this->SetPrimaryOutputName( "MetaDataDictionaryObject" );
  this->ProcessObject::SetOutput( "MetaDataDictionaryObject",
    this->MakeOutput( 0 ) );
}


template< typename TInputImage >
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::~MarkDuplicateFramesInvalidImageFilter( void )
{
}


template< typename TInputImage >
void
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::SetInput( const InputImageType * image )
{
  this->ProcessObject::SetPrimaryInput( const_cast< InputImageType * >( image ) );
}


template< typename TInputImage >
const typename MarkDuplicateFramesInvalidImageFilter< TInputImage >::InputImageType *
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::GetInput( void ) const
{
  return itkDynamicCastInDebugMode< const TInputImage * >( this->GetPrimaryInput() );
}


template< typename TInputImage >
void
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::SetInputMetaDataDictionary( const MetaDataDictionary * dictionary )
{
  if( dictionary != m_InputMetaDataDictionary )
    {
    m_InputMetaDataDictionary = dictionary;
    this->Modified();
    }
}


template< typename TInputImage >
const MetaDataDictionary *
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::GetInputMetaDataDictionary( void ) const
{
  return m_InputMetaDataDictionary;
}


template< typename TInputImage >
const MetaDataDictionary &
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::GetOutputMetaDataDictionary( void ) const
{
  const DataObject * outputDataObject =
    this->GetOutput( "MetaDataDictionaryObject" );
  return outputDataObject->GetMetaDataDictionary();
}


template< typename TInputImage >
DataObject::Pointer
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::MakeOutput( DataObjectPointerArraySizeType )
{
  // Just use any type of Object that has a MetaDataDictionary.
  return ImageBase< 3 >::New().GetPointer();
}


template< typename TInputImage >
void
MarkDuplicateFramesInvalidImageFilter< TInputImage >
::GenerateData( void )
{
  const InputImageType * inputImage = this->GetInput();
  this->m_Threader->Execute( this, inputImage->GetBufferedRegion() );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMarkDuplicateFramesInvalidImageFilter_hxx )
