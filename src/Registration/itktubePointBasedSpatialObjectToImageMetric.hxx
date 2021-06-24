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

  this->ComputeCenterOfRotation();
}


template< unsigned int ObjectDimension, class TFixedImage >
unsigned int
PointBasedSpatialObjectToImageMetric< ObjectDimension, TFixedImage >
::GetMaximumNumberOfPoints( void )
{
  unsigned int maximumNumberOfPoints = 0;

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
          ++maximumNumberOfPoints;
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
          ++maximumNumberOfPoints;
          }
        ++pIter;
        }
      }
    }
  delete pbsoList;

  return maximumNumberOfPoints;
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
::SubsamplePoints( void ) 
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
      double weight = ( radius - m_TubeSamplingRadiusMin ) /
        (m_TubePriorityRadius - m_TubeSamplingRadiusMin);
      m_SubsampledTubePointsWeights.push_back( weight );
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
::GetPointBasedChildren( SpatialObjectType::Pointer & parentSO, 
  SpatialObjectType::ChildrenListType * childrenSO ) const
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
    if( this->IsValidMovingPoint( inputPoint )
      {
      const InputPointType inputPoint = tubePointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double scale = tubePointIter->GetRadiusInWorldSpace() * m_Kappa;

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
            n1, n2, scale );
          }
        else
          {
          matchMeasure += *tubePointWeightIter * imFunc->Evaluate( currentPoint,
            n1, scale );
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
    if( this->IsValidMovingPoint( inputPoint )
      {
      const InputPointType inputPoint =
        surfacePointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double scale = surfacePointIter->GetRadiusInWorldSpace() * m_Kappa;

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
          currentPoint, n1, scale );
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
    if( this->IsValidMovingPoint( inputPoint )
      {
      const InputPointType inputPoint =
        pointIter->GetPositionInWorldSpace();
      OutputPointType currentPoint;
      transformCopy->TransformPoint( inputPoint, currentPoint );
      if( IsValidFixedPoint( currentPoint ) )
        {
        double scale = pointIter->GetRadiusInWorldSpace() * m_Kappa;

        matchMeasure += *pointWeightIter * imFunc->Evaluate(
          currentPoint, scale );
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

  VnlMatrixType biasV( TubeDimension, TubeDimension,
    NumericTraits< double >::Zero );
  VnlMatrixType biasVI( TubeDimension, TubeDimension,
    NumericTraits< double >::Zero );

  SizeValueType weightCount = 0;

  CompensatedSummationType dPosition[TubeDimension];

  VnlMatrixType tM( TubeDimension, TubeDimension );
  VnlVectorType v1T( TubeDimension );
  VnlVectorType v2T( TubeDimension );

  typedef std::vector< OutputPointType > TubePointsContainerType;
  TubePointsContainerType transformedTubePoints;
  transformedTubePoints.reserve( m_FeatureWeights.GetSize() );

  typedef std::vector< VectorType >      dTubePointsContainerType;
  dTubePointsContainerType dtransformedTubePoints;
  dtransformedTubePoints.reserve( m_FeatureWeights.GetSize() );

  derivative.fill( 0.0 );

  typename TubeTreeType::ChildrenListType * tubeList = GetTubes();
  typename TubeTreeType::ChildrenListType::const_iterator tubeIter;
  for( tubeIter = tubeList->begin();
       tubeIter != tubeList->end();
       ++tubeIter )
    {
    TubeType* currentTube = dynamic_cast<TubeType*>(
      ( *tubeIter ).GetPointer() );

    if( currentTube != NULL )
      {
      typename TubeType::TubePointListType::const_iterator pointIter;
      for( pointIter = currentTube->GetPoints().begin();
           pointIter != currentTube->GetPoints().end();
           ++pointIter )
        {
        InputPointType inputPoint = pointIter->GetPositionInObjectSpace();
        OutputPointType currentPoint;
        if( this->IsInside( inputPoint, currentPoint, transformCopy ) )
          {
          transformedTubePoints.push_back( currentPoint );

          CovariantVectorType v1;
          CovariantVectorType v2;
          for( unsigned int ii = 0; ii < TubeDimension; ++ii )
            {
            v1[ii] = pointIter->GetNormal1InObjectSpace()[ii];
            v2[ii] = pointIter->GetNormal2InObjectSpace()[ii];
            }

          for( unsigned int ii = 0; ii < TubeDimension; ++ii )
            {
            v1T[ii] = transformCopy->TransformCovariantVector( v1 )[ii];
            v2T[ii] = transformCopy->TransformCovariantVector( v2 )[ii];
            }
          tM = outer_product( v1T, v1T );
          tM = tM + outer_product( v2T, v2T );
          tM = m_FeatureWeights[weightCount] * tM;
          biasV += tM;

          v1 = transformCopy->TransformCovariantVector( v1 );
          v2 = transformCopy->TransformCovariantVector( v2 );

          double scalingRadius = pointIter->GetRadiusInObjectSpace();
          scalingRadius = std::max( scalingRadius, m_MinimumScalingRadius );

          const double scale = scalingRadius * m_Kappa;

          const double dXProj1
            = this->ComputeThirdDerivatives( v1, scale, currentPoint );
          const double dXProj2
            = this->ComputeThirdDerivatives( v2, scale, currentPoint );

          VectorType dtransformedTubePoint;
          for( unsigned int ii = 0; ii < TubeDimension; ++ii )
            {
            dtransformedTubePoint[ii] =
              ( dXProj1 * v1[ii] + dXProj2 * v2[ii] );
            dPosition[ii] +=
              m_FeatureWeights[weightCount] * ( dtransformedTubePoint[ii] );
            }
          dtransformedTubePoints.push_back( dtransformedTubePoint );
          }
        ++weightCount;
        }
      }
    }

  biasVI = vnl_matrix_inverse< double >( biasV ).inverse();

  VnlVectorType tV( TubeDimension );
  for( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    tV[ii] = dPosition[ii].GetSum();
    }

  tV *= biasVI;

  for( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    dPosition[ii] = tV[ii];
    }

  CompensatedSummationType dAngle[TubeDimension];
  VectorType offsets;
  for( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    offsets[ii] = dPosition[ii].GetSum();
    }
  typename TubePointsContainerType::const_iterator transformedTubePointsIt =
    transformedTubePoints.begin();
  typename dTubePointsContainerType::const_iterator
    dtransformedTubePointsIt = dtransformedTubePoints.begin();
  weightCount = 0;
  for( tubeIter = tubeList->begin();
       tubeIter != tubeList->end();
       ++tubeIter )
    {
    TubeType* currentTube = dynamic_cast<TubeType*>(
      ( *tubeIter ).GetPointer() );

    if( currentTube != NULL )
      {
      typename TubeType::TubePointListType::const_iterator pointIter;
      for( pointIter = currentTube->GetPoints().begin();
           pointIter != currentTube->GetPoints().end();
           ++pointIter )
        {
        InputPointType inputPoint = pointIter->GetPositionInObjectSpace();
        OutputPointType currentPoint;
        //! \todo this checking is not necessary, but we have to make sure
        //have
        //the corresponding weight -- should cache the result and use it
        //here.
        if( this->IsInside( inputPoint, currentPoint, transformCopy ) )
          {
          VnlVectorType dXT( TubeDimension );
          for( unsigned int ii = 0; ii < TubeDimension; ++ii )
            {
            dXT[ii] = ( *dtransformedTubePointsIt )[ii];
            }

          dXT = dXT * biasVI;

          double angleDelta[TubeDimension];
          this->GetDeltaAngles( *transformedTubePointsIt, dXT, offsets,
            angleDelta );
          for( unsigned int ii = 0; ii < TubeDimension; ++ii )
            {
            dAngle[ii] += m_FeatureWeights[weightCount] * angleDelta[ii];
            }
          ++dtransformedTubePointsIt;
          ++transformedTubePointsIt;
          }
        ++weightCount;
        }
      }
    }

  derivative[0] = dAngle[0].GetSum();
  derivative[1] = dAngle[1].GetSum();
  derivative[2] = dAngle[2].GetSum();
  derivative[3] = dPosition[0].GetSum();
  derivative[4] = dPosition[1].GetSum();
  derivative[5] = dPosition[2].GetSum();

  delete tubeList;
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
::IsValidMovingPoint( TubePointType * pnt )
{
  
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubePointBasedSpatialObjectToImageMetric_hxx )
