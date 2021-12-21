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

#ifndef __itktubeExtractTubePointsSpatialObjectFilter_hxx
#define __itktubeExtractTubePointsSpatialObjectFilter_hxx


#include "tubeTubeMathFilters.h"

namespace itk
{

namespace tube
{

template< class TTubeSpatialObject >
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::ExtractTubePointsSpatialObjectFilter( void )
{
  this->ProcessObject::SetNumberOfRequiredOutputs( 1 );
  this->MakeOutput( 0 );
}


template< class TTubeSpatialObject >
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::~ExtractTubePointsSpatialObjectFilter( void )
{
}


template< class TTubeSpatialObject >
ProcessObject::DataObjectPointer
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::MakeOutput( ProcessObject::DataObjectPointerArraySizeType
  itkNotUsed( idx ) )
{
  this->m_PointsContainerDecorator = PointsContainerDecoratorType::New();
  this->m_PointsContainer = PointsContainerType::New();

  this->m_PointsContainerDecorator->Set( this->m_PointsContainer );

  this->ProcessObject::SetNthOutput( 0,
    m_PointsContainerDecorator.GetPointer() );

  return this->m_PointsContainerDecorator.GetPointer();
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
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
  ::GroupSpatialObjectType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetInput( void ) const
{
  return itkDynamicCastInDebugMode< const GroupSpatialObjectType * >
    ( this->GetPrimaryInput() );
}


template< class TTubeSpatialObject >
typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
  ::PointsContainerDecoratorType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainerOutput( void )
{
  return this->m_PointsContainerDecorator.GetPointer();
}


template< class TTubeSpatialObject >
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
  ::PointsContainerDecoratorType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainerOutput( void ) const
{
  return this->m_PointsContainerDecorator.GetPointer();
}


template< class TTubeSpatialObject >
const typename ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
  ::PointsContainerType *
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GetPointsContainer( void ) const
{
  return this->m_PointsContainer.GetPointer();
}


template< class TTubeSpatialObject >
void
ExtractTubePointsSpatialObjectFilter< TTubeSpatialObject >
::GenerateData( void )
{
  const GroupSpatialObjectType * inputGroup = this->GetInput();

  char childName[] = "Tube";
  typedef typename TubeSpatialObjectType::ChildrenListType ChildrenListType;
  ChildrenListType * childrenList =
    inputGroup->GetChildren( inputGroup->GetMaximumDepth(), childName );

  for( typename ChildrenListType::const_iterator childrenIt =
    childrenList->begin(); childrenIt != childrenList->end();
    ++childrenIt )
    {
    typename TubeSpatialObjectType::Pointer tube;
    tube = dynamic_cast< TubeSpatialObjectType * >(
      childrenIt->GetPointer() );
    if( tube.IsNotNull() )
      {
      tube->RemoveDuplicatePointsInObjectSpace();
      tube->ComputeTangentsAndNormals();
      const typename TubeSpatialObjectType::TubePointListType
        pointsForThisTube = tube->GetPoints();
      unsigned int count = m_PointsContainer->Size();
      m_PointsContainer->Reserve( count +
        pointsForThisTube.size() );
      for( unsigned int i=0; i<pointsForThisTube.size(); ++i )
        {
        m_PointsContainer->SetElement( count++, pointsForThisTube[i] );
        }
      }
    }

  delete childrenList;
}


} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeExtractTubePointsSpatialObjectFilter_hxx )
