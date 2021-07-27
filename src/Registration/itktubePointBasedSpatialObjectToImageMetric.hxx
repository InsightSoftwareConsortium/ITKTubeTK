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

#ifndef __itktubePointBasedSpatialObjectToImageMetric_hxx
#define __itktubePointBasedSpatialObjectToImageMetric_hxx

#include "itktubePointBasedSpatialObjectToImageMetric.h"
#include "itktubeNJetImageFunction.h"

namespace itk
{

namespace tube
{

template< unsigned int ObjectDimension, class TFixedImage >
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::PointBasedSpatialObjectToImageMetric( void )
{
  m_SamplingRatio = 0.1;

  m_Kappa = 1.0;
  m_Extent = 3.0;

  m_TubePriorityRadius = 0;
  m_TubeSamplingRadiusMin = 0;
  m_TubeSamplingRadiusMax = 0;

  m_CenterOfRotation.Fill( 0.0 );

  m_SubsampledSingularPoints.clear();
  m_SubsampledSingularPointsWeights.clear();
  m_SubsampledTubePoints.clear();
  m_SubsampledTubePointsWeights.clear();
  m_SubsampledSurfacePoints.clear();
  m_SubsampledSurfacePointsWeights.clear();
}


template< unsigned int ObjectDimension, class TFixedImage >
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::~PointBasedSpatialObjectToImageMetric( void )
{
}

template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::Initialize( void )
{
  if( !this->m_MovingSpatialObject || !this->m_FixedImage )
    {
    itkExceptionMacro( << "No tube/image net plugged in." );
    }
}


template< unsigned int ObjectDimension, class TFixedImage >
unsigned int
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetMaximumNumberOfPoints( void )
{
  m_MaximumNumberOfPoints = 0;
  m_MaximumNumberOfTubePoints = 0;
  m_MaximumNumberOfSurfacePoints = 0;

  typename PointBasedSpatialObjectType::ChildrenListType * pbsoList =
    this->GetPointBasedChildren( m_MovingSpatialObject );
  typename PointBasedSpatialObjectType::ChildrenListType::const_iterator
    pbsoIter;
  for( pbsoIter = pbsoList->begin(); pbsoIter != pbsoList->end();
    ++pbsoIter )
    {
    TubeType* currentTube =
      dynamic_cast<TubeType*>((*pbsoIter).GetPointer());
    SurfaceType* currentSurface =
      dynamic_cast<SurfaceType*>((*pbsoIter).GetPointer());
    PointBasedSpatialObjectType* currentSingular =
      dynamic_cast<SingularObjectType*>((*pbsoIter).GetPointer());
    if( currentTube != nullptr )
      {
      const typename TubeType::TubePointListType & currentPoints =
        currentTube->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          ++m_MaximumNumberOfTubePoints;
          }
        ++pIter;
        }
      }
    else if( currentSurface != NULL )
      {
      const typename SurfaceType::SurfacePointListType & currentPoints =
        currentSurface->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          ++m_MaximumNumberOfSurfacePoints;
          }
        ++pIter;
        }
      }
    else
      {
      const typename SingularType::SpatialObjectPointListType &
        currentPoints = currentSingular->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          ++m_MaximumNumberOfPoints;
          }
        ++pIter;
        }
      }
    }
  delete pbsoList;

  return m_MaximumNumberOfPoints +
    m_MaximumNumberOfTubePoints +
    m_MaximumNumberOfSurfacePoints;
}


template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::ComputeCenterOfRotation( void )
{
  unsigned int pointCount = 0;

  typename PointBasedSpatialObjectType::ChildrenConstListType * pbsoList =
    this->GetPointBasedChildren( m_MovingSpatialObject );
  auto pbsoIter = pbsoList->begin();
  while( pbsoIter != pbsoList->end() )
  {
    const TubeType* currentTube =
      dynamic_cast<const TubeType*>((*pbsoIter).GetPointer());
    const SurfaceType* currentSurface =
      dynamic_cast<const SurfaceType*>((*pbsoIter).GetPointer());
    const PointBasedSpatialObjectType* currentSingular =
      dynamic_cast<const PointBasedSpatialObjectType*>((*pbsoIter).GetPointer());
    if( currentTube != nullptr )
      {
      const typename TubeType::TubePointListType & currentPoints =
        currentTube->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
            {
            this->m_CenterOfRotation[ii] +=
              pIter->GetPositionInWorldSpace()[ii];
            }
          ++pointCount;
          }
        ++pIter;
        }
      }
    else if( currentSurface != nullptr )
      {
      const typename SurfaceType::SurfacePointListType & currentPoints =
        currentSurface->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
            {
            this->m_CenterOfRotation[ii] +=
              (pIter->GetPositionInWorldSpace())[ii];
            }
          ++pointCount;
          }
        ++pIter;
        }
      }
    else
      {
      const typename SingularType::SpatialObjectPointListType &
        currentPoints = currentSingular->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
            {
            this->m_CenterOfRotation[ii] +=
              pIter->GetPositionInWorldSpace()[ii];
            }
          ++pointCount;
          }
        ++pIter;
        }
      }
    ++pbsoIter;
    }
  delete pbsoList;

  for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
    {
    this->m_CenterOfRotation[ii] /= pointCount;
    }
}


template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::ComputeSubsampledPoints( void ) 
{
  unsigned int maxPointCount = this->GetMaximumNumberOfPoints();

  unsigned int targetPointCount = maxPointCount * m_SamplingRatio;
  
  m_SubsampledSingularPoints.clear();
  m_SubsampledSingularPointsWeights.clear();
  m_SubsampledTubePoints.clear();
  m_SubsampledTubePointsWeights.clear();
  m_SubsampledSurfacePoints.clear();
  m_SubsampledSurfacePointsWeights.clear();

  unsigned int pointCount = 0;
  typename PointBasedSpatialObjectType::ChildrenConstListType * pbsoList =
    this->GetPointBasedSpatialObjects();
  auto pbsoIter = pbsoList->begin();
  while( pbsoIter != pbsoList->end() )
    {
    const TubeType* currentTube =
      dynamic_cast<const TubeType*>((*pbsoIter).GetPointer());
    const SurfaceType* currentSurface =
      dynamic_cast<const SurfaceType*>((*pbsoIter).GetPointer());
    const PointBasedSpatialObjectType* currentSingular =
      dynamic_cast<const SingularType*>((*pbsoIter).GetPointer());
    unsigned int preCount;
    unsigned int postCount;
    PointIteratorType pointIter = currentPoints.begin();
    if( currentTube != nullptr )
      {
      const TubePointListType & currentPoints = currentTube->GetPoints();
      auto pointIter = currentPoints.begin();
      while( pointIter != currentPoints.end() )
        {
        m_SubsampledTubePoints.push_back( *pointIter );
        m_SubsampledTubePointsWeights.push_back(1);
        preCount = pointCount * m_SamplingRatio;
        postCount = preCount;
        while( preCount == postCount && pointIter != currentPoints.end() )
          {
          ++pointIter;
          ++pointCount;
          postCount = pointCount * m_SamplingRatio;
          }
        }
      }
    else if( currentSurface != nullptr )
      {
      const SurfacePointListType & currentPoints = currentSurface->GetPoints();
      auto pointIter = currentPoints.begin();
      while( pointIter != currentPoints.end() )
        {
        m_SubsampledSurfacePoints.push_back( *pointIter );
        m_SubsampledSurfacePointsWeights.push_back(1);
        preCount = pointCount * m_SamplingRatio;
        postCount = preCount;
        while( preCount == postCount && pointIter != currentPoints.end() )
          {
          ++pointIter;
          ++pointCount;
          postCount = pointCount * m_SamplingRatio;
          }
        }
      }
    else
      {
      const SingularPointListType & currentPoints = currentSingular->GetPoints();
      auto pointIter = currentPoints.begin();
      while( pointIter != currentPoints.end() )
        {
        m_SubsampledSingularPoints.push_back( *pointIter );
        m_SubsampledSingularPointsWeights.push_back(1);
        preCount = pointCount * m_SamplingRatio;
        postCount = preCount;
        while( preCount == postCount && pointIter != currentPoints.end() )
          {
          ++pointIter;
          ++pointCount;
          postCount = pointCount * m_SamplingRatio;
          }
        }
      }
    ++pbsoIter;
    }
  delete pbsoList;
}

template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::ComputeSubsampledPointsWeights( void ) 
{
  m_SubsampledTubePointsWeights.clear();
  auto tubePointIter = m_SubsampledTubePoints.begin();
  while( tubePointIter != m_SubsampledTubePoints.end() )
  {
    double radius = tubePointIter->GetRadiusInWorldSpace();
    if( m_TubePriorityRadius != 0 &&
      m_TubePriorityRadius > m_TubeSamplingRadiusMin)
      {
      if( radius > m_TubePriorityRadius )
        {
        m_SubsampledTubePointsWeights.push_back(1);
        }
      else
        {
        double weight = ( radius - m_TubeSamplingRadiusMin ) /
          (m_TubePriorityRadius - m_TubeSamplingRadiusMin);
        m_SubsampledTubePointsWeights.push_back( weight );
        }
      }
    else
      {
      if( m_TubeSamplingRadiusMin < m_TubeSamplingRadiusMax )
        {
        double weight = ( radius - m_TubeSamplingRadiusMin ) /
          (m_TubeSamplingRadiusMax - m_TubeSamplingRadiusMin);
        m_SubsampledTubePointsWeights.push_back( weight );
        }
      else
        {
        m_SubsampledTubePointsWeights.push_back( 1 );
        }
      }
    ++tubePointIter;
  }
  m_SubsampledSurfacePointsWeights.clear();
  auto surfacePointIter = m_SubsampledSurfacePoints.begin();
  while( surfacePointIter != m_SubsampledSurfacePoints.end() )
  {
    m_SubsampledSurfacePointsWeights.push_back( 1 );
    ++surfacePointIter;
  }
  m_SubsampledSingularPointsWeights.clear();
  auto pointIter = m_SubsampledSingularPoints.begin();
  while( pointIter != m_SubsampledSingularPoints.end() )
  {
    m_SubsampledSingularPointsWeights.push_back( 1 );
    ++pointIter;
  }
}

template< unsigned int ObjectDimension, class TFixedImage >
typename
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >::MovingSpatialObjectType::ChildrenConstListType*
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetPointBasedChildren( typename MovingSpatialObjectType::ConstPointer & parentSO, 
  typename MovingSpatialObjectType::ChildrenConstListType * childrenSO ) const
{
  if( parentSO.IsNull() )
  {
    return childrenSO;
  }

  if( childrenSO == nullptr )
  {
    childrenSO = new MovingSpatialObjectType::ChildrenConstListType;
  }

  auto tmpChildren = parentSO->GetConstChildren();
  auto it = tmpChildren->begin();
  while (it != tmpChildren->end())
  {
    if (dynamic_cast< const PointBasedSpatialObjectType * >((*it).GetPointer()) != nullptr)
    {
      childrenSO->push_back(*it);
    }
    it++;
  }

  it = tmpChildren->begin();
  while (it != tmpChildren->end())
  {
    this->GetPointBasedChildren( *it, childrenSO );
    it++;
  }

  delete tmpChildren;

  return childrenSO;
}

template< unsigned int ObjectDimension, class TFixedImage >
typename PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >::MeasureType
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetValue( const ParametersType & parameters ) const
{
  typedef vnl_matrix_fixed<double, ImageDimension, ImageDimension> VnlMatrixType;
  typedef vnl_vector_fixed<double, ImageDimension> VnlVectorType;

  itkDebugMacro( << "**** Get Value ****" );
  itkDebugMacro( << "Parameters = " << parameters );

  double value = 0;

  double weightSum = 0;

  // Create a copy of the transform to keep true const correctness (
  // thread-safe )
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  NJetImageFunction< FixedImageType >::Pointer imFunc =
    NJetImageFunction< FixedImageType >::New();
  imFunc->SetInputImage( m_FixedImage );

  VnlVectorType n1T, n2T;

  typename TubePointWeightListType::const_iterator pointWeightIter
    = m_SubsampledTubePointsWeights.begin();
  typename TubePointListType::const_iterator pointIter =
    m_SubsampledTubePoints.begin();
  while( pointIter != m_SubsampledTubePoints.end() )
  {
    MovingPointType movingPoint = pointIter->GetPositionInWorldSpace();
    FixedPointType fixedPoint = transformCopy->TransformPoint( movingPoint );
    if( this->IsValidFixedPoint( fixedPoint ) )
    {
      double radius = pointIter->GetRadiusInWorldSpace();
      radius = radius * m_Kappa;
      MovingPointType radiusPoint;
      radiusPoint.Fill(radius);
      FixedPointType fixedRadiusPoint =
        transformCopy->TransformPoint( radiusPoint );
      double fixedScale = radiusPoint[0];
      for( unsigned int i=1; i<ImageDimension; ++i )
      {
        fixedScale += radiusPoint[i];
      }
      fixedScale /= ImageDimension;

      typename TubeType::CovariantVectorType n1 = pointIter->GetNormal1InWorldSpace();
      n1 = transformCopy->TransformCovariantVector( n1 );
      Vector<double,ImageDimension> n1v;
      n1v.SetVnlVector(n1.GetVnlVector());
      n1T = n1.GetVnlVector();
      if( ObjectDimension > 2 )
      {
        typename TubeType::CovariantVectorType n2 = pointIter->GetNormal2InWorldSpace();
        n2 = transformCopy->TransformCovariantVector( n2 );
        Vector<double,ImageDimension> n2v;
        n2v.SetVnlVector(n2.GetVnlVector());
        n2T = n2.GetVnlVector();
        value += imFunc->Evaluate( fixedPoint, n1v, n2v, fixedScale );
      }
      else
      {
        value += imFunc->Evaluate( fixedPoint, n1v, fixedScale );
      }
      value = *pointWeightIter * imFunc->GetMostRecentIntensity();
      weightSum += *pointWeightIter;
    }
    ++pointIter;
    ++pointWeightIter;
  }

  if( weightSum > 0 )
  {
    value /= weightSum;
  }

  return value;
}

template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetDerivative( const ParametersType & parameters, DerivativeType & derivative ) const
{
  MeasureType value;
  this->GetValueAndDerivative(parameters, value, derivative );
}

template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetValueAndDerivative( const ParametersType & parameters, MeasureType & value, DerivativeType & derivative ) const
{
  typedef vnl_matrix_fixed<double, ImageDimension, ImageDimension> VnlMatrixType;
  typedef vnl_vector_fixed<double, ImageDimension> VnlVectorType;

  itkDebugMacro( << "**** Get Derivative ****" );
  itkDebugMacro( << "Parameters = " << parameters );

  value = 0;

  double weightSum = 0;

  // Create a copy of the transform to keep true const correctness (
  // thread-safe )
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  NJetImageFunction< FixedImageType >::Pointer imFunc =
    NJetImageFunction< FixedImageType >::New();
  imFunc->SetInputImage( m_FixedImage );

  VnlMatrixType biasV;
  VnlVectorType n1T, n2T;
  VnlMatrixType tM;

  std::list< FixedPointType > fixedPoints;
  std::list< FixedVectorType > fixedPointsDerivs;

  typename TubePointWeightListType::const_iterator pointWeightIter
    = m_SubsampledTubePointsWeights.begin();
  typename TubePointListType::const_iterator pointIter =
    m_SubsampledTubePoints.begin();
  while( pointIter != m_SubsampledTubePoints.end() )
  {
    MovingPointType movingPoint = pointIter->GetPositionInWorldSpace();
    FixedPointType fixedPoint = transformCopy->TransformPoint( movingPoint );
    if( this->IsValidFixedPoint( fixedPoint ) )
    {
      fixedPoints.push_back( fixedPoint );

      FixedVectorType fixedDeriv;

      double radius = pointIter->GetRadiusInWorldSpace();
      radius = radius * m_Kappa;
      MovingPointType radiusPoint;
      radiusPoint.Fill(radius);
      FixedPointType fixedRadiusPoint =
        transformCopy->TransformPoint( radiusPoint );
      double fixedScale = radiusPoint[0];
      for( unsigned int i=1; i<ImageDimension; ++i )
      {
        fixedScale += radiusPoint[i];
      }
      fixedScale /= ImageDimension;

      typename TubeType::CovariantVectorType n1 = pointIter->GetNormal1InWorldSpace();
      n1 = transformCopy->TransformCovariantVector( n1 );
      Vector<double,ImageDimension> n1v;
      n1v.SetVnlVector(n1.GetVnlVector());
      n1T = n1.GetVnlVector();
      tM = outer_product( n1T, n1T );
      if( ObjectDimension > 2 )
      {
        typename TubeType::CovariantVectorType n2 = pointIter->GetNormal2InWorldSpace();
        n2 = transformCopy->TransformCovariantVector( n2 );
        Vector<double,ImageDimension> n2v;
        n2v.SetVnlVector(n2.GetVnlVector());
        n2T = n2.GetVnlVector();
        tM = tM + outer_product( n2T, n2T );
        Vector<double, ImageDimension> tmpDeriv;
        imFunc->Derivative( fixedPoint, n1v, n2v, fixedScale, tmpDeriv );
        fixedDeriv.fill(0);
        for( unsigned int ii = 0; ii < ImageDimension; ++ii )
        {
          fixedDeriv[ii] += tmpDeriv[0] * n1T[ii];
          fixedDeriv[ii] += tmpDeriv[1] * n2T[ii];
        }
      }
      else
      {
        Vector<double, ImageDimension> tmpDeriv;
        imFunc->Derivative( fixedPoint, n1v, fixedScale, tmpDeriv );
        fixedDeriv.fill(0);
        for( unsigned int ii = 0; ii < ImageDimension; ++ii )
        {
          fixedDeriv[ii] += tmpDeriv[0] * n1T[ii];
        }
      }
      fixedPointsDerivs.push_back( fixedDeriv );

      tM = *pointWeightIter * tM;
      biasV += tM;

      weightSum += *pointWeightIter;
    }
    ++pointIter;
    ++pointWeightIter;
  }

  VnlMatrixType biasVI = vnl_matrix_inverse< double >( biasV ).inverse();

  derivative.fill(0);

  auto fixedPointsIter = fixedPoints.begin();
  auto fixedPointsDerivsIter = fixedPointsDerivs.begin();
  while( fixedPointsIter != fixedPoints.end() )
  {
    FixedPointType fixedPoint = *fixedPointsIter;

    VnlVectorType dXT = (*fixedPointsDerivsIter);

    dXT = dXT * biasVI;

    typename TransformType::JacobianType jacobian;
    m_Transform->ComputeJacobianWithRespectToParameters( fixedPoint,
      jacobian );

    for( unsigned int p = 0; p<m_Transform->GetNumberOfParameters(); ++p)
    {
      for( unsigned int d=0; d<ImageDimension; ++d)
      {
        derivative[p] += jacobian[p][d] * dXT[d];
      }
    }

    ++fixedPointsIter;
    ++fixedPointsDerivsIter;
  }
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const TubePointType & pnt ) const
{
  if( m_TubeSamplingRadiusMin != 0 || m_TubeSamplingRadiusMax != 0 )
    {
    double radius = pnt.GetRadiusInWorldSpace();
    if( radius < m_TubeSamplingRadiusMin || radius > m_TubeSamplingRadiusMax )
      {
      return false;
      }
    }
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInsideInWorldSpace(
      pnt.GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
  return true;
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const SurfacePointType & pnt ) const
{
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInsideInWorldSpace(
      pnt.GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
  return true;
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const SingularPointType & pnt ) const
{
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInsideInWorldSpace(
      pnt.GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
  return true;
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidFixedPoint( const FixedPointType & pnt ) const
{
  if( m_FixedImageMaskObject.IsNotNull() && m_UseFixedImageMaskObject )
    {
    if( !m_FixedImageMaskObject->IsInsideInWorldSpace( pnt ) )
      {
      return false;
      }
    }
  return true;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_hxx )
