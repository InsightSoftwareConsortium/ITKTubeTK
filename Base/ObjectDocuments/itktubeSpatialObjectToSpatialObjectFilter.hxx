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

/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itktubeSpatialObjectToSpatialObjectFilter.hxx,v $
  Language:  C++
  Date:      $Date: 2004/08/11 20:59:14 $
  Version:   $Revision: 1.2 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef __itktubeSpatialObjectToSpatialObjectFilter_hxx
#define __itktubeSpatialObjectToSpatialObjectFilter_hxx

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

namespace itk
{

namespace tube
{

template< class TInputSpatialObject, class TOutputSpatialObject >
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::SpatialObjectToSpatialObjectFilter( void )
{
  this->SetNumberOfRequiredInputs( 1 );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::~SpatialObjectToSpatialObjectFilter( void )
{
}

template< class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::GetInput( void )
{
  if( this->GetNumberOfInputs() < 1 )
    {
    return NULL;
    }

  return static_cast< const InputSpatialObjectType * >( this->Superclass::GetInput( 0 ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
const typename SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >::InputSpatialObjectType *
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::GetInput( unsigned int index )
{
  return static_cast< const InputSpatialObjectType * >( this->Superclass::GetInput( index ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::SetInput( const InputSpatialObjectType * input )
{
  this->Superclass::SetNthInput( 0, const_cast< InputSpatialObjectType * >( input ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::SetInput( unsigned int index, const InputSpatialObjectType * object )
{
  if( index + 1 > this->GetNumberOfInputs() )
    {
    this->SetNumberOfRequiredInputs( index + 1 );
    }

  this->Superclass::SetNthInput( index, const_cast< InputSpatialObjectType * >( object ) );
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::GenerateOutputInformation( void )
{
}

template< class TInputSpatialObject, class TOutputSpatialObject >
void
SpatialObjectToSpatialObjectFilter< TInputSpatialObject, TOutputSpatialObject >
::GenerateData( void )
{
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeSpatialObjectToSpatialObjectFilter_hxx)
