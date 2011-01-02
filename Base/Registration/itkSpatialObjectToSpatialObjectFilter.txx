/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkSpatialObjectToSpatialObjectFilter_txx
#define __itkSpatialObjectToSpatialObjectFilter_txx
#include "itkSpatialObjectToSpatialObjectFilter.h"

namespace itk
{

/**
 *
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::SpatialObjectToSpatialObjectFilter()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs( 1 );
}

/**
 *
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::~SpatialObjectToSpatialObjectFilter()
{
}

/**
 *
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::SetInput( const InputSpatialObjectType *input )
{
  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( 0,
    const_cast< InputSpatialObjectType * >( input ) );
}

/**
 * Connect one of the operands for pixel-wise addition
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::SetInput( unsigned int index, const TInputSpatialObject * object )
{
  if( index+1 > this->GetNumberOfInputs() )
    {
    this->SetNumberOfRequiredInputs( index + 1 );
    }

  // Process object is not const-correct so the const_cast is required here
  this->ProcessObject::SetNthInput( index,
    const_cast< TInputSpatialObject * >( object ) );
}

/**
 *
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::GetInput( void )
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return 0;
    }

  return static_cast< const TInputSpatialObject * >(
    this->ProcessObject::GetInput( 0 ) );
}

/**
 *
 */
template < class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::GetInput( unsigned int idx )
{
  return static_cast< const TInputSpatialObject * >(
    this->ProcessObject::GetInput( idx ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject,
  TOutputSpatialObject >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
}

} // end namespace itk

#endif
