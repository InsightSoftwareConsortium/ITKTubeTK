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

#ifndef __itktubePointBasedSpatialObjectTransformFilter_hxx
#define __itktubePointBasedSpatialObjectTransformFilter_hxx

#include "itktubePointBasedSpatialObjectTransformFilter.h"

#include <itkSpatialObjectFactory.h>

namespace itk
{

namespace tube
{

template< class TTransformType, unsigned int TDimension >
PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
::PointBasedSpatialObjectTransformFilter( void )
{
  m_OutputObjectToParentTransform = 0;
  m_Transform = 0;

  SpatialObjectFactoryBase::RegisterDefaultSpatialObjects();
  SpatialObjectFactory< SpatialObject< TDimension > >::
    RegisterSpatialObject();
  SpatialObjectFactory< TubeType >::RegisterSpatialObject();
  SpatialObjectFactory< SurfaceType >::RegisterSpatialObject();
  SpatialObjectFactory< LineType >::RegisterSpatialObject();
  SpatialObjectFactory< DTITubeType >::RegisterSpatialObject();
  SpatialObjectFactory< ContourType >::RegisterSpatialObject();
}

/**
 * Apply the transformation to the tube
 */
template< class TTransformType, unsigned int TDimension >
void
PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
::GenerateData( void )
{

  typename SpatialObjectType::Pointer output = this->GetInput()->Clone();
  this->Transform( this->GetInput(), output );

  SpatialObjectType * soOutput =
    static_cast<SpatialObjectType *>( this->ProcessObject::GetOutput(0) );
  //soOutput->Set( output.GetPointer() );

  typedef typename SpatialObject< TDimension >::ChildrenListType
    ChildrenListType;
  ChildrenListType * children = this->GetInput()->GetChildren();
  typename ChildrenListType::const_iterator it = children->begin();
  while( it != children->end() )
    {
    this->UpdateLevel( *it, output );
    ++it;
    }
  delete children;
}

/**
 * Apply the transformation to the tube
 */
template< class TTransformType, unsigned int TDimension >
void
PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
::UpdateLevel( SpatialObject< TDimension > * inputSO,
  SpatialObject< TDimension > * parentSO )
{
  typename SpatialObject< TDimension >::Pointer outputSO =
    inputSO->Clone();
  if( outputSO.IsNull() )
    {
    itkExceptionMacro( << "Could not create an instance of "
      << outputSO->GetTypeName() << ". The usual cause of this error is not"
      << "registering the SpatialObject with SpatialFactory." );
    }

  this->Transform( inputSO, outputSO );

  parentSO->AddChild( outputSO );

  typedef typename SpatialObject< TDimension >::ChildrenListType
    ChildrenListType;
  ChildrenListType * children = inputSO->GetChildren();
  typename ChildrenListType::const_iterator it = children->begin();
  while( it != children->end() )
    {
    this->UpdateLevel( *it, outputSO );
    ++it;
    }
  delete children;
}

template< class TTransformType, unsigned int TDimension >
void
PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
::Transform( const SpatialObject< TDimension > * inputSO,
  SpatialObject< TDimension > * outputSO )
{
  // We make the copy and sub-sample if it is a tube.
  const PointBasedType * inputSOAsPointBased =
    dynamic_cast< const PointBasedType * >( inputSO );
  if( inputSOAsPointBased != NULL )
    {
    PointBasedType * outputSOAsPointBased = dynamic_cast< PointBasedType * >(
      outputSO );

    Point<double, TDimension> worldPoint;
    Point<double, TDimension> transformedWorldPoint;

    //inputSOAsPointBased->Update();

    outputSOAsPointBased->CopyInformation( inputSOAsPointBased );
    outputSOAsPointBased->Clear();

    if( m_OutputObjectToParentTransform.IsNotNull() )
      {
      typename SpatialObjectTransformType::Pointer tfm =
        SpatialObjectTransformType::New();
      tfm->SetIdentity();
      tfm->SetMatrix( m_OutputObjectToParentTransform->GetMatrix() );
      tfm->SetOffset( m_OutputObjectToParentTransform->GetOffset() );
      outputSO->SetObjectToParentTransform( tfm );
      }

    //inputSOAsPointBased->Update();

    const TubeType * inputSOAsTube =
      dynamic_cast< const TubeType * >( inputSO );
    if( inputSOAsTube != NULL )
      {
      TubeType * outputSOAsTube = dynamic_cast< TubeType * >( outputSO );
      typedef typename TubeType::TubePointListType      TubePointListType;
      TubePointListType tubePointList = inputSOAsTube->GetPoints();
      typename TubePointListType::const_iterator tubePointIterator =
        tubePointList.begin();

      while( tubePointIterator != tubePointList.end() )
        {
        worldPoint = ( *tubePointIterator ).GetPositionInWorldSpace();

        transformedWorldPoint = m_Transform->TransformPoint( worldPoint );

        TubeSpatialObjectPoint<TDimension> pnt;

        pnt.SetId( tubePointIterator->GetId() );
        pnt.SetColor( tubePointIterator->GetColor() );

        pnt.SetPositionInWorldSpace( transformedWorldPoint );

        // get both normals
        typename TubeType::CovariantVectorType n1 = tubePointIterator
          ->GetNormal1InWorldSpace();
        typename TubeType::CovariantVectorType n2 = tubePointIterator
          ->GetNormal2InWorldSpace();

        // only try transformation of normals if both are non-zero
        if( !n1.GetVnlVector().is_zero() && !n2.GetVnlVector().is_zero() )
          {
          n1 = m_Transform->TransformCovariantVector( n1, worldPoint );
          n2 = m_Transform->TransformCovariantVector( n2, worldPoint );
          n1.Normalize();
          n2.Normalize();
          pnt.SetNormal1InWorldSpace( n1 );
          pnt.SetNormal2InWorldSpace( n2 );
          }

        typename TubeType::VectorType tang = tubePointIterator->
          GetTangentInWorldSpace();
        if( !tang.GetVnlVector().is_zero() )
          {
          tang = m_Transform->TransformVector( tang, worldPoint );
          tang.Normalize();
          pnt.SetTangentInWorldSpace( tang );
          }

        typename TubeType::VectorType radi;
        for( unsigned int i=0; i<TDimension; ++i )
          {
          radi[i] = tubePointIterator->GetRadiusInWorldSpace();
          }
        radi = m_Transform->TransformVector( radi, worldPoint );
        pnt.SetRadiusInWorldSpace( radi[0] );

        pnt.SetMedialness( ( *tubePointIterator ).GetMedialness() );
        pnt.SetRidgeness( ( *tubePointIterator ).GetRidgeness() );
        pnt.SetBranchness( ( *tubePointIterator ).GetBranchness() );

        outputSOAsTube->AddPoint( pnt );

        ++tubePointIterator;
        }
      }
    else
      {
      const SurfaceType * inputSOAsSurface = dynamic_cast< const SurfaceType * >( inputSO );
      if( inputSOAsSurface != NULL )
        {
        SurfaceType * outputSOAsSurface = dynamic_cast< SurfaceType * >( 
          outputSO );
        typedef typename SurfaceType::SurfacePointListType SurfacePointListType;
        SurfacePointListType surfacePointList = inputSOAsSurface->GetPoints();
        typename SurfacePointListType::const_iterator surfacePointIterator =
          surfacePointList.begin();

        while( surfacePointIterator != surfacePointList.end() )
          {
          worldPoint = ( *surfacePointIterator ).GetPositionInWorldSpace();
  
          transformedWorldPoint = m_Transform->TransformPoint( worldPoint );
  
          SurfaceSpatialObjectPoint<TDimension> pnt;
  
          pnt.SetId( surfacePointIterator->GetId() );
          pnt.SetColor( surfacePointIterator->GetColor() );
  
          pnt.SetPositionInWorldSpace( transformedWorldPoint );
  
          // get normals
          typename SurfaceType::CovariantVectorType n1 = surfacePointIterator
            ->GetNormalInWorldSpace();
  
          // only try transformation of normals if both are non-zero
          if( !n1.GetVnlVector().is_zero() )
            {
            n1 = m_Transform->TransformCovariantVector( n1, worldPoint );
            n1.Normalize();
            pnt.SetNormalInWorldSpace( n1 );
            }
  
          outputSOAsSurface->AddPoint( pnt );
  
          ++surfacePointIterator;
          }
        }
      else
        {
        std::cerr <<
          "WARNING: Transforms currently only applied to Tubes and Surfaces. " 
          << std::endl;
        }
      }
    }
}

template< class TTransformType, unsigned int TDimension >
void
PointBasedSpatialObjectTransformFilter< TTransformType, TDimension >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Transformation: " << m_Transform << std::endl;
  os << indent << "OutputObjectToParent Transform: " <<
    m_OutputObjectToParentTransform << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePointBasedSpatialObjectTransformFilter_hxx )
