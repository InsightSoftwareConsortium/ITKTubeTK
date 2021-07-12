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

  m_SubsampledPoints.clear();
  m_SubsampledPointsWeights.clear();
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
    PointBasedSpatialObjectType* currentPBSO =
      dynamic_cast<PointBasedSpatialObjectType*>((*pbsoIter).GetPointer());
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
      const typename PointBasedSpatialObjectType::SpatialObjectPointListType &
        currentPoints = currentPBSO->GetPoints();
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
    PointBasedSpatialObjectType* currentPBSO =
      dynamic_cast<PointBasedSpatialObjectType*>((*pbsoIter).GetPointer());
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
              (tubePointIter->GetPositionInWorldSpace())[ii];
            }
          ++pointCount;
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
          ++maximumNumberOfPoints;
          for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
            {
            this->m_CenterOfRotation[ii] +=
              (tubePointIter->GetPositionInWorldSpace())[ii];
            }
          ++pointCount;
          }
        ++pIter;
        }
      }
    else
      {
      const typename PointBasedSpatialObjectType::SpatialObjectPointListType &
        currentPoints = currentPBSO->GetPoints();
      auto pIter = currentPoints.begin();
      while( pIter != currentPoints.end() )
        {
        if( this->IsValidMovingPoint( *pIter ) )
          {
          for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
            {
            this->m_CenterOfRotation[ii] +=
              (tubePointIter->GetPositionInWorldSpace())[ii];
            }
          ++pointCount;
          }
        ++pIter;
        }
      }
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
  
  m_SubsampledPoints.clear();
  m_SubsampledPointsWeights.clear();
  m_SubsampledTubePoints.clear();
  m_SubsampledTubePointsWeights.clear();
  m_SubsampledSurfacePoints.clear();
  m_SubsampledSurfacePointsWeights.clear();

  unsigned int pointCount = 0;
  typename PointBasedSpatialObjectType::ChildrenListType * pbsoList =
    this->GetPointBasedSpatialObjects();
  typename PointBasedSpatialObjectType::ChildrenListType::const_iterator
    pbsoIter;
  pbsoIter = pbsoList->begin();
  while( pbsoIter != pbsoList->end() )
    {
    TubeType* currentTube =
      dynamic_cast<TubeType*>((*pbsoIter).GetPointer());
    SurfaceType* currentSurface =
      dynamic_cast<SurfaceType*>((*pbsoIter).GetPointer());
    PointBasedSpatialObjectType* currentPBSO =
      dynamic_cast<PointBasedSpatialObjectType*>((*pbsoIter).GetPointer());
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
      const SpatialObjectPointListType & currentPoints = currentPBSO->GetPoints();
      auto pointIter = currentPoints.begin();
      while( pointIter != currentPoints.end() )
        {
        m_SubsampledPoints.push_back( *pointIter );
        m_SubsampledPointsWeights.push_back(1);
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
  m_SubsampledPointsWeights.clear();
  auto pointIter = m_SubsampledPoints.begin();
  while( pointIter != m_SubsampledPoints.end() )
  {
    m_SubsampledPointsWeights.push_back( 1 );
    ++pointIter;
  }
}

template< unsigned int ObjectDimension, class TFixedImage >
typename
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >::SpatialObjectType::ChildrenListType*
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetPointBasedChildren( typename SpatialObjectType::Pointer & parentSO, 
  typename SpatialObjectType::ChildrenListType * childrenSO ) const
{
  if( parentSO.IsNull() )
  {
    return childrenSO;
  }

  if( childrenSO == nullptr )
  {
    childrenSO = new SpatialObjectType::ChildrenListType;
  }

  auto tmpChildren = parentSO->GetChildren();

  auto it = tmpChildren->begin();
  while (it != tmpChildren->end())
  {
    if (dynamic_cast< PointBasedSpatialObjectType * >(*it) != nullptr)
    {
      childrenSO->push_back((*it));
    }
    it++;
  }

  auto it = tmpChildren->begin();
  while (it != tmpChildren->end())
  {
    this->AddPointBasedChildrenToList( *it, childrenSO );
    it++;
  }

  delete tmpChildren;

  return childrenSO;
}

template< unsigned int ObjectDimension, class TFixedImage >
typename PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >::
  MeasureType
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro( << "**** Get Value ****" );
  itkDebugMacro( << "Parameters = " << parameters );

  double matchMeasure;
  double weightSum;

  // Create a copy of the transform to keep true const correctness (
  // thread-safe )
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  typename NJetImageFunction< FixedImageType >::Pointer imFunc =
    NJetImageFunction< FixedImageType >::New();
  imFunc->SetInput( m_FixedImage );

  auto tubePointIter = m_SubsampledTubePoints.begin();
  auto tubePointWeightIter = m_SubsampledTubePointsWeights.begin();
  while( tubePointIter != m_SubsampledTubePoints.end() )
    {
    if( this->IsValidMovingPoint( inputPoint ) )
      {
      const InputPointType inputPoint = tubePointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double radius = pointIter->GetRadiusInWorldSpace();
        radius = radius * m_Kappa;
        InputPointType radiusPoint.Fill(radius);
        OutputPointType currentRadiusPoint =
          transformCopy->TransformPoint( radiusPoint );
        double currentScale = radiusPoint[0];
        for( unsigned int i=1; i<ImageDimension; ++i )
          {
          currentScale += radiusPoint[i];
          }
        currentScale /= ImageDimension;

        CovariantVectorType n1;
        TransformedCovariantVectorType n1T;
        CovariantVectorType n2;
        TransformedCovariantVectorType n2T;
        for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
          {
          n1[ii] = pointIter->GetNormal1InWorldSpace()[ii];
          if( ObjectDimension > 2 )
            {
            n2[ii] = pointIter->GetNormal2InWorldSpace()[ii];
            }
          }

        for( unsigned int ii = 0; ii < ImageDimension; ++ii )
          {
          n1T[ii] = transformCopy->TransformCovariantVector( n1 )[ii];
          if( ObjectDimension > 2 )
            {
            n2T[ii] = transformCopy->TransformCovariantVector( n2 )[ii];
            }
          }

        if( ObjectDimension > 2 )
          {
          matchMeasure += *tubePointWeightIter * imFunc->Evaluate( currentPoint,
            n1, n2, currentScale );
          }
        else
          {
          matchMeasure += *tubePointWeightIter * imFunc->Evaluate( currentPoint,
            n1, currentScale );
          }
        weightSum += *tubePointWeightIter;
        }
      }
    ++tubePointIter;
    ++tubePointWeightIter;
    }

  auto surfacePointIter = m_SubsampledSurfacePoints.begin();
  auto surfacePointWeightIter = m_SubsampledSurfacePointsWeights.begin();
  while( surfacePointIter != m_SubsampledSurfacePoints.end() )
    {
    if( this->IsValidMovingPoint( inputPoint ) )
      {
      const InputPointType inputPoint =
        surfacePointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double radius = pointIter->GetRadiusInWorldSpace();
        radius = radius * m_Kappa;
        InputPointType radiusPoint.Fill(radius);
        OutputPointType currentRadiusPoint =
          transformCopy->TransformPoint( radiusPoint );
        double currentScale = radiusPoint[0];
        for( unsigned int i=1; i<ImageDimension; ++i )
          {
          currentScale += radiusPoint[i];
          }
        currentScale /= ImageDimension;

        CovariantVectorType n1;
        TransformedCovariantVectorType n1T;
        for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
          {
          n1[ii] = pointIter->GetNormal1InWorldSpace()[ii];
          }

        for( unsigned int ii = 0; ii < ImageDimension; ++ii )
          {
          n1T[ii] = transformCopy->TransformCovariantVector( n1 )[ii];
          }

        matchMeasure += *surfacePointWeightIter * imFunc->Evaluate(
          currentPoint, n1, currentScale );
        weightSum += *surfacePointWeightIter;
        }
      }
    ++surfacePointIter;
    ++surfacePointWeightIter;
    }

  auto pointIter = m_SubsampledPoints.begin();
  auto pointWeightIter = m_SubsampledPointsWeights.begin();
  while( pointIter != m_SubsampledSurfacePoints.end() )
    {
    if( this->IsValidMovingPoint( inputPoint ) )
      {
      const InputPointType inputPoint =
        pointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double radius = pointIter->GetRadiusInWorldSpace();
        radius = radius * m_Kappa;
        InputPointType radiusPoint.Fill(radius);
        OutputPointType currentRadiusPoint =
          transformCopy->TransformPoint( radiusPoint );
        double currentScale = radiusPoint[0];
        for( unsigned int i=1; i<ImageDimension; ++i )
          {
          currentScale += radiusPoint[i];
          }
        currentScale /= ImageDimension;

        matchMeasure += *pointWeightIter * imFunc->Evaluate(
          currentPoint, currentScale );
        weightSum += *pointWeightIter;
        }
      }
    ++pointIter;
    ++pointWeightIter;
    }

  if( weightSum == NumericTraits< double >::Zero )
    {
    itkWarningMacro(
      << "GetValue: All the transformed tube points are outside the image." );
    matchMeasure = 0;
    }
  else
    {
    matchMeasure /= weightSum;
    }

  itkDebugMacro( << "matchMeasure = " << matchMeasure );

  delete tubeList;

  return matchMeasure;
}

template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  derivative.SetSize( this->GetNumberOfParameters() );
  derivative.fill( 0.0 );

  // Create a copy of the transform to keep true const correctness
  // ( thread-safe )
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform =
    this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  itkDebugMacro( << "**** Get Derivative ****" );
  itkDebugMacro( << "parameters = "<< parameters );

  vnl_matrix<double> biasV( ObjectDimension, ObjectDimension,
    NumericTraits< double >::Zero );
  vnl_matrix<double> biasVI( ObjectDimension, ObjectDimension,
    NumericTraits< double >::Zero );

  double weightSum = 0;
  double dMagSum = 0;

  vnl_matrix<double> tM( ObjectDimension, ObjectDimension );
  vnl_vector<double> n1T( ObjectDimension );
  vnl_vector<double> n2T( ObjectDimension );

  int totalNumberOfPoints = m_SubsampledPoints.GetSize()
    + m_SubsampledTubePoints.GetSize() + m_SubsampledSurfacePoints.GetSize();

  std::vector< OutputPointType > FixedPointListType;
  FixedPointListType fixedPoints;
  fixedPoints.reserve( totalNumberOfPoints );

  typedef std::vector< outputVectorType > FixedPointDerivListType;
  FixedPointDerivListType  fixedPointsDerivs;
  fixedPointsDerivs.reserve( totalNumberOfPoints );

  typename NJetImageFunction< FixedImageType >::Pointer imFunc =
    NJetImageFunction< FixedImageType >::New();
  imFunc->SetInput( m_FixedImage );

  typename PointWeightListType::const_iterator pointWeightIter;
  typename TubePointListType::const_iterator pointIter =
    m_TubePointList.begin();
  while( pointIter != m_TubePointList.end() )
    {
    InputPointType inputPoint = pointIter->GetPositionInWorldSpace();
    OutputPointType currentPoint = transformCopy->TransformPoint( inputPoint );
    if( this->IsValidFixedPoint( currentPoint ) )
      {
      VectorType currentDeriv;

      fixedPoints.push_back( currentPoint );

      double radius = pointIter->GetRadiusInWorldSpace();
      radius = radius * m_Kappa;
      InputPointType radiusPoint.Fill(radius);
      OutputPointType currentRadiusPoint =
        transformCopy->TransformPoint( radiusPoint );
      double currentScale = radiusPoint[0];
      for( unsigned int i=1; i<ImageDimension; ++i )
        {
        currentScale += radiusPoint[i];
        }
      currentScale /= ImageDimension;

      CovariantVectorType n1 = pointIter->GetNormal1InWorldSpace();
      n1 = transformCopy->TransformCovariantVector( n1 );
      n1T = n1.GetVnlVector();
      tM = outer_product( n1T, n1T );
      if( ObjectDimension > 2 )
        {
        CovariantVectorType n2 = pointIter->GetNormal2InWorldSpace();
        n2 = transformCopy->TransformCovariantVector( n2 );
        n2T = n2.GetVnlVector();
        tM = tM + outer_product( n2T, n2T );
        VectorType tmpDeriv;
        dMagSum += *pointWeightIter * imFunc->Derivative( currentPoint,
          n1T, n2T, currentScale, tmpDeriv );
        currentDeriv.fill(0);
        for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
          {
          currentDeriv[ii] += tmpDeriv[0] * n1T[ii];
          currentDeriv[ii] += tmpDeriv[1] * n2T[ii];
          }
        }
      else
        {
        VectorType tmpDeriv;
        dMagSum += *pointWeightIter * imFunc->Derivative( currentPoint,
          n1T, currentScale, tmpDeriv );
        currentDeriv.fill(0);
        for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
          {
          currentDeriv[ii] += tmpDeriv[0] * n1T[ii];
          }
        }
      fixedPointsDerivs.push_back( currentDeriv );

      tM = *pointWeightIter * tM;
      biasV += tM;

      weightSum += *pointWeightIter;
      }
    ++pointIter;
    ++pointWeightIter;
    }

  biasVI = vnl_matrix_inverse< double >( biasV ).inverse();

  auto fixedPointsIter = fixedPoints.begin();
  auto fixedPointsDerivsIter = fixedPointsDerivs.begin();
  while( fixedPointsIter != fixedPoints.end() )
    {
    OutputPointType currentPoint = *fixedPointIter;

    VnlVectorType dXT( ObjectDimension );
    for( unsigned int ii = 0; ii < ObjectDimension; ++ii )
      {
      dXT[ii] = ( *fixedPointsDerivsIter )[ii];
      }

    dXT = dXT * biasVI;

    typename TransformType::JacobianType jacobian;
    m_Transform->ComputeJacobianWithRespectToParameters( currentPoint,
      jacobian );

    DerivativeType tmpDeriv = jacobian * dXT;

    derivative += tmpDeriv;

    ++fixedPointsIter;
    ++fixedPointsDerivsIter;
    }
}


template< unsigned int ObjectDimension, class TFixedImage >
void
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType & value,
                         DerivativeType & derivative ) const
{
  value = this->GetValue( parameters );
  this->GetDerivative( parameters, derivative );
}


template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const TubePointType & pnt ) const
{
  if( m_TubeSamplingRadiusMin != 0 || m_TubeSamplingRadiusMax != 0 )
    {
    double radius = pnt->GetRadiusInWorldSpace();
    if( radius < m_TubeSamplingRadiusMin || radius > m_TubeSamplingRadiusMax )
      {
      return false;
      }
    }
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInside(
      pnt->GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const SurfacePointType & pnt ) const
{
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInside(
      pnt->GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidMovingPoint( const SpatialObjectPointType & pnt ) const
{
  if( m_MovingSpatialObjectMaskObject.IsNotNull() &&
    m_UseMovingSpatialObjectMaskObject )
    {
    if( !m_MovingSpatialObjectMaskObject->IsInside(
      pnt->GetPositionInWorldSpace() ) )
      {
      return false;
      }
    }
}

template< unsigned int ObjectDimension, class TFixedImage >
bool
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::IsValidFixedPoint( const OutputPointType & pnt ) const
{
  if( m_FixedImageMaskObject.IsNotNull() && m_UseFixedImageMaskObject )
    {
    if( !m_FixedImageMaskObject->IsInside( pnt ) )
      {
      return false;
      }
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_hxx )
