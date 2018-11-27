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

#ifndef __itktubeVectorImageToListGenerator_hxx
#define __itktubeVectorImageToListGenerator_hxx

#include "itktubeVectorImageToListGenerator.h"

#include <itkImageRegionConstIterator.h>

namespace itk
{

namespace tube
{

namespace Statistics
{

template< class TImage, class TMaskImage >
VectorImageToListGenerator< TImage, TMaskImage >
::VectorImageToListGenerator( void )
{
  m_UseSingleMaskValue = false;
  m_MaskValue = 1;
  this->SetNumberOfRequiredInputs( 1 );
  this->SetNumberOfRequiredOutputs( 1 );
  typename ListSampleOutputType::Pointer listSampleDecorator =
    static_cast< ListSampleOutputType * >( this->MakeOutput( 0 ).GetPointer() );
  this->ProcessObject::SetNthOutput( 0, listSampleDecorator.GetPointer() );
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "MaskValue: "
     << static_cast<typename NumericTraits<MaskPixelType>::PrintType>(
       m_MaskValue )
     << std::endl;
  if( m_UseSingleMaskValue )
    {
    os << indent << "UseSingleMaskValue: True" << std::endl;
    }
  else
    {
    os << indent << "UseSingleMaskValue: False" << std::endl;
    }
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::SetInput( const ImageType* image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
                                   const_cast< ImageType* >( image ) );
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::SetMaskImage( const MaskImageType* image )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 1,
                                   const_cast< MaskImageType* >( image ) );
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::SetMaskValue( const MaskPixelType maskValue )
{
  m_MaskValue = maskValue;
  m_UseSingleMaskValue = true;
  this->Modified();
}

template< class TImage, class TMaskImage >
const TImage*
VectorImageToListGenerator< TImage, TMaskImage >
::GetInput( void ) const
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return 0;
    }

  return static_cast<const ImageType * >
    ( this->ProcessObject::GetInput( 0 ) );
}

template< class TImage, class TMaskImage >
const TMaskImage*
VectorImageToListGenerator< TImage, TMaskImage >
::GetMaskImage( void ) const
{
  if( this->GetNumberOfInputs() < 2 )
    {
    return 0;
    }

  return static_cast<const MaskImageType * >
    ( this->ProcessObject::GetInput( 1 ) );
}

template< class TImage, class TMaskImage >
typename VectorImageToListGenerator< TImage, TMaskImage >::DataObjectPointer
VectorImageToListGenerator< TImage, TMaskImage >
::MakeOutput( DataObjectPointerArraySizeType itkNotUsed( idx ) )
{
  typename ListSampleOutputType::Pointer decoratedOutput =
    ListSampleOutputType::New();
  decoratedOutput->Set( ListSampleType::New() );
  return static_cast< DataObject * >( decoratedOutput.GetPointer() );
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::GenerateData( void )
{
  ListSampleOutputType * decoratedOutput =
    static_cast< ListSampleOutputType * >(
      this->ProcessObject::GetOutput( 0 ) );

  ListSampleType *output =
    const_cast< ListSampleType * >( decoratedOutput->Get() );

  const ImageType *input = this->GetInput();
  MaskImageType *maskImage = NULL;

  output->Clear();

  if( this->GetNumberOfInputs() > 1 )
    {
    maskImage = const_cast< MaskImageType * >( this->GetMaskImage() );
    }

  typedef ImageRegionConstIterator< ImageType >     IteratorType;
  IteratorType it( input, input->GetBufferedRegion() );
  it.GoToBegin();

  if( maskImage ) // mask specified
    {
    if( m_UseSingleMaskValue )
      {
      typedef ImageRegionConstIterator< MaskImageType > MaskIteratorType;
      MaskIteratorType mit( maskImage, maskImage->GetBufferedRegion() );
      mit.GoToBegin();
      while( !it.IsAtEnd() )
        {
        if( mit.Get() == this->m_MaskValue )
          {
          MeasurementVectorType m;
          m = it.Get();
          output->PushBack( m );
          }
        ++mit;
        ++it;
        }
      }
    else
      {
      typedef ImageRegionConstIterator< MaskImageType > MaskIteratorType;
      MaskIteratorType mit( maskImage, maskImage->GetBufferedRegion() );
      mit.GoToBegin();
      while( !it.IsAtEnd() )
        {
        if( mit.Get() != 0 )
          {
          MeasurementVectorType m;
          m = it.Get();
          output->PushBack( m );
          }
        ++mit;
        ++it;
        }
      }
    }
  else // no mask specified
    {
    while( !it.IsAtEnd() )
      {
      MeasurementVectorType m;
      m = it.Get();
      output->PushBack( m );
      ++it;
      }
    }
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::GenerateOutputInformation( void )
{
  Superclass::GenerateOutputInformation();

  ListSampleOutputType * decoratedOutput =
    static_cast< ListSampleOutputType * >(
      this->ProcessObject::GetOutput( 0 ) );
  ListSampleType *output =
   const_cast< ListSampleType *>( decoratedOutput->Get() );
  output->SetMeasurementVectorSize(
    itkGetStaticConstMacro( MeasurementVectorSize ) );
}

template< class TImage, class TMaskImage >
void
VectorImageToListGenerator< TImage, TMaskImage >
::GenerateInputRequestedRegion() throw( InvalidRequestedRegionError )
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // Make sure that the mask's requested region, if specified is at least
  // as large as the input image's buffered region. If not funny things can
  // happen such as the mask iterator going out of bounds etc..
  //
  // TODO: Why don't most other ITK filters that take multiple inputs check
  // for this ?
  //
  if( this->GetNumberOfInputs() > 1 )
    {
    MaskImageType *maskImage = const_cast< MaskImageType * >(
      this->GetMaskImage() );
    ImageType     *image = const_cast< ImageType * >(
      this->GetInput() );
    if( !image->GetBufferedRegion().IsInside(
      maskImage->GetBufferedRegion() ) )
      {
      maskImage->SetRequestedRegion( image->GetBufferedRegion() );
      }
    }
}

template< class TImage, class TMaskImage >
const typename VectorImageToListGenerator< TImage, TMaskImage >
::ListSampleType *
VectorImageToListGenerator< TImage, TMaskImage >
::GetListSample( void ) const
{
  const ListSampleOutputType * decoratedOutput =
    static_cast< const ListSampleOutputType * >(
      this->ProcessObject::GetOutput( 0 ) );
  return decoratedOutput->Get();
}

} // End namespace Statistics

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeVectorImageToListGenerator_hxx )
