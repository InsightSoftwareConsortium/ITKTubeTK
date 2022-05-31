/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __itktubeSpatialObjectFilter_hxx
#define __itktubeSpatialObjectFilter_hxx


namespace itk
{

namespace tube
{

template< unsigned int ObjectDimension >
SpatialObjectFilter< ObjectDimension >
::SpatialObjectFilter( void )
{
}

template< unsigned int ObjectDimension >
void
SpatialObjectFilter< ObjectDimension >
::SetInput( const SpatialObject<ObjectDimension> * input )
{
  // Process object is not const-correct so the const_cast is required here
  this->SpatialObjectSource< SpatialObject<ObjectDimension> >::SetNthInput( 0,
   const_cast<SpatialObject<ObjectDimension> *>(input) );

  typename SpatialObject<ObjectDimension>::Pointer output =
   static_cast<SpatialObject<ObjectDimension> *>(this->MakeOutput( 0 ).GetPointer());
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}

template< unsigned int ObjectDimension >
ProcessObject::DataObjectPointer
SpatialObjectFilter< ObjectDimension >
::MakeOutput( ProcessObject::DataObjectPointerArraySizeType itkNotUsed(
  idx ) )
{
  return this->GetInput()->Clone().GetPointer();
}



} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeSpatialObjectFilter_hxx )
