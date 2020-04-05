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

#ifndef __itktubeSpatialObjectToSpatialObjectFilter_hxx
#define __itktubeSpatialObjectToSpatialObjectFilter_hxx

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

#include <itkTextOutput.h>

namespace itk
{

namespace tube
{

template< class TInputSpatialObject, class TOutputSpatialObject >
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::SpatialObjectToSpatialObjectFilter( void )
{
  this->SetNumberOfRequiredInputs( 1 );

  itk::OutputWindow::SetInstance( itk::TextOutput::New() );
}


template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::SetInput( const InputSpatialObjectType * input )
{
  // Process object is not const-correct so the const_cast is required here
  this->SpatialObjectSource< TOutputSpatialObject >::SetNthInput( 0,
    const_cast< TInputSpatialObject * >( input ) );
}


template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::SetInput( unsigned int index, const InputSpatialObjectType * input )
{
  // Process object is not const-correct so the const_cast is required here
  this->SpatialObjectSource< TOutputSpatialObject >::SetNthInput( index,
    const_cast< TInputSpatialObject * >( input ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::SetInput( const itk::ProcessObject::DataObjectIdentifierType & index,
  itk::DataObject * input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetInput( index, input );
}


template< class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter < TInputSpatialObject,
TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::GetInput( void ) const
{
  return itkDynamicCastInDebugMode< const TInputSpatialObject * >(
    this->GetPrimaryInput() );
}


template< class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter < TInputSpatialObject,
TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
TOutputSpatialObject >
::GetInput( unsigned int index ) const
{
  const TInputSpatialObject * input = dynamic_cast< const
    TInputSpatialObject * >( this->SpatialObjectSource<
    TOutputSpatialObject >::GetInput( index ) );

  if( input == NULL && this->SpatialObjectSource< TOutputSpatialObject >::
    GetInput( input ) != NULL )
    {
    itkWarningMacro( << "Unable to convert input number " << index
      << " to type " << typeid( InputSpatialObjectType ).name () );
    }
  return input;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSpatialObjectToSpatialObjectFilter_hxx )
