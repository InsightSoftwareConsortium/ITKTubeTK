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

#ifndef __itktubeExtractTubePointsSpatialObjectFilter_hxx
#define __itktubeExtractTubePointsSpatialObjectFilter_hxx

#include "itktubeExtractTubePointsSpatialObjectFilter.h"

namespace itk
{

namespace tube
{

template< class TTubeSpatialObject >
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::ExtractTubePointsSpatialObjectFilter( void )
{
  // Create the output. We use static_cast<> here because we know the default
  // output must be of type PointsContainerDecoratorType.
  typename PointsContainerDecoratorType::Pointer output =
    static_cast< PointsContainerDecoratorType * >( this->MakeOutput( 0 ).GetPointer() );
  this->m_PointsContainer = PointsContainerType::New();
  output->Set( this->m_PointsContainer );

  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->ProcessObject::SetNthOutput( 0, output.GetPointer() );
}


template< class TTubeSpatialObject >
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::~ExtractTubePointsSpatialObjectFilter( void )
{
}


template< class TTubeSpatialObject >
ProcessObject::DataObjectPointer
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::MakeOutput( ProcessObject::DataObjectPointerArraySizeType itkNotUsed( idx ) )
{
  return PointsContainerDecoratorType::New().GetPointer();
}


template< class TTubeSpatialObject >
void
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::SetInput( const GroupSpatialObjectType * group )
{
  // ProcessObject is not const-correct so a const_cast is required here.
  this->ProcessObject::SetNthInput( 0,
                                    const_cast< GroupSpatialObjectType * >( group ) );
}


template< class TTubeSpatialObject >
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >::GroupSpatialObjectType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetInput( void ) const
{
  return itkDynamicCastInDebugMode< const GroupSpatialObjectType * >( this->GetPrimaryInput() );
}


template< class TTubeSpatialObject >
typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >::PointsContainerDecoratorType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainerOutput( void )
{
  // we assume that the first output is of the templated type
  return itkDynamicCastInDebugMode< PointsContainerDecoratorType * >
    ( this->GetPrimaryOutput() );
}


template< class TTubeSpatialObject >
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >::PointsContainerDecoratorType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainerOutput( void ) const
{
  // we assume that the first output is of the templated type
  return itkDynamicCastInDebugMode< const PointsContainerDecoratorType * >
    ( this->GetPrimaryOutput() );
}


template< class TTubeSpatialObject >
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >::PointsContainerType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainer( void ) const
{
  const PointsContainerDecoratorType * output = this->GetPointsContainerOutput();
  return output->Get();
}


template< class TTubeSpatialObject >
void
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GenerateData( void )
{
  typename PointsContainerType::STLContainerType & pointsContainer =
    m_PointsContainer->CastToSTLContainer();
  pointsContainer.clear();

  const GroupSpatialObjectType * inputGroup = this->GetInput();

  char childName[] = "Tube";
  typedef typename TubeSpatialObjectType::ChildrenListType ChildrenListType;
  ChildrenListType * childrenList =
    inputGroup->GetChildren( inputGroup->GetMaximumDepth(), childName );

  for( typename ChildrenListType::const_iterator childrenIt = childrenList->begin();
       childrenIt != childrenList->end();
       ++childrenIt )
    {
    TubeSpatialObjectType * tube =
      dynamic_cast< TubeSpatialObjectType * >( ( *childrenIt ).GetPointer() );
    if( tube != NULL )
      {
      tube->RemoveDuplicatePoints();
      tube->ComputeTangentAndNormals();
      const typename TubeSpatialObjectType::PointListType pointsForThisTube =
        tube->GetPoints();
      pointsContainer.insert( pointsContainer.end(),
        pointsForThisTube.begin(),
        pointsForThisTube.end() );
      }
    }

  delete childrenList;
}


} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeExtractTubePointsSpatialObjectFilter_hxx)
