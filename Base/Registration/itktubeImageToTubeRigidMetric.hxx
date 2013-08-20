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

#ifndef __itktubeImageToTubeRigidMetric_hxx
#define __itktubeImageToTubeRigidMetric_hxx

#include "itktubeImageToTubeRigidMetric.h"

#include <itkLinearInterpolateImageFunction.h>

namespace itk
{

namespace tube
{

template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::ImageToTubeRigidMetric( void )
{
  m_Kappa = 1.0;
  m_Extent = 3.0;
  m_ImageMin = 0.0;
  m_ImageMax = 0.0;

  m_CenterOfRotation.Fill( 0.0 );

  m_DerivativeImageFunction = DerivativeImageFunctionType::New();

  typedef LinearInterpolateImageFunction< FixedImageType > DefaultInterpolatorType;
  this->m_Interpolator = DefaultInterpolatorType::New();
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
                        TResolutionWeightFunction >
::~ImageToTubeRigidMetric( void )
{
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject,
  TTubeSpatialObject, TResolutionWeightFunction >::ResolutionWeightFunctionType &
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetResolutionWeightFunction( void )
{
  return this->m_ResolutionWeightFunction;
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
const typename ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject,
  TTubeSpatialObject, TResolutionWeightFunction >::ResolutionWeightFunctionType &
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetResolutionWeightFunction( void ) const
{
  return this->m_ResolutionWeightFunction;
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::SetResolutionWeightFunction( const ResolutionWeightFunctionType & function )
{
  if( this->m_ResolutionWeightFunction != function )
    {
    this->m_ResolutionWeightFunction = function;
    this->Modified();
    }
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::Initialize( void ) throw ( ExceptionObject )
{
  m_ResolutionWeights.clear();

  if( !this->m_MovingSpatialObject || !this->m_FixedImage )
    {
    itkExceptionMacro( << "No tube/image net plugged in." );
    }

  this->ComputeTubePointResolutionWeights();

  this->m_Interpolator->SetInputImage( this->m_FixedImage );
  this->m_DerivativeImageFunction->SetInputImage( this->m_FixedImage );
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeTubePointResolutionWeights( void )
{
  typename TubeNetType::ChildrenListType* tubeList = this->GetTubes();

  CompensatedSummationType resolutionWeightSum;

  typedef typename TubeNetType::ChildrenListType::iterator TubesIteratorType;
  for( TubesIteratorType tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeType* currentTube =
      dynamic_cast< TubeType * >( ( *tubeIterator ).GetPointer() );
    if( currentTube != NULL )
      {
      currentTube->RemoveDuplicatePoints();
      currentTube->ComputeTangentAndNormals();

      const typename TubeType::PointListType & currentTubePoints
        = currentTube->GetPoints();
      typedef typename TubeType::PointListType::const_iterator
        TubePointIteratorType;
      for( TubePointIteratorType tubePointIterator = currentTubePoints.begin();
            tubePointIterator != currentTubePoints.end();
            ++tubePointIterator )
        {
        const ScalarType weight =
          this->m_ResolutionWeightFunction( *tubePointIterator );

        this->m_ResolutionWeights.push_back( weight );

        for( unsigned int ii = 0; ii < TubeDimension; ++ii )
          {
          this->m_CenterOfRotation[ii] +=
            weight * ( tubePointIterator->GetPosition() )[ii];
          }
        resolutionWeightSum += weight;
        }
      }
    }
  delete tubeList;

  for( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    this->m_CenterOfRotation[ii] /= resolutionWeightSum.GetSum();
    }

  //! \todo partial template specialization for transforms with a center?
  TransformType * transform = static_cast< TransformType * >
    ( this->m_Transform.GetPointer() );
  transform->SetCenter( this->m_CenterOfRotation );
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject,
  TTubeSpatialObject, TResolutionWeightFunction >::TubeNetType::ChildrenListType*
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetTubes( void ) const
{
  if( !this->m_MovingSpatialObject )
    {
    return NULL;
    }

  char childName[] = "Tube";
  return this->m_MovingSpatialObject->GetChildren(
    this->m_MovingSpatialObject->GetMaximumDepth(), childName );
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject,
  TTubeSpatialObject, TResolutionWeightFunction >::MeasureType
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro( << "**** Get Value ****" );
  itkDebugMacro( << "Parameters = " << parameters );

  CompensatedSummationType matchMeasure;
  CompensatedSummationType weightSum;

  ResolutionWeightsContainerType::const_iterator weightIterator;
  weightIterator = this->m_ResolutionWeights.begin();

  // Create a copy of the transform to keep true const correctness (thread-safe)
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  typename TubeNetType::ChildrenListType::iterator tubeIterator;
  typename TubeNetType::ChildrenListType * tubeList = GetTubes();
  for( tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeType * currentTube =
      dynamic_cast< TubeType * >( ( *tubeIterator ).GetPointer() );
    if( currentTube != NULL )
      {
      typename TubeType::PointListType::iterator pointIterator;
      for( pointIterator = currentTube->GetPoints().begin();
           pointIterator != currentTube->GetPoints().end();
           ++pointIterator )
        {
        const InputPointType inputPoint = pointIterator->GetPosition();
        OutputPointType currentPoint;
        if( this->IsInside( inputPoint, currentPoint, transformCopy ) )
          {
          weightSum += *weightIterator;
          ScalarType scalingRadius = pointIterator->GetRadius();
          // !TODO 0.5 should be a parameter of the class
          scalingRadius = std::max( scalingRadius, 0.5 );

          const ScalarType scale = scalingRadius * m_Kappa;

          matchMeasure += *weightIterator * vnl_math_abs(
            this->ComputeLaplacianMagnitude( pointIterator->GetNormal1(),
              scale,
              currentPoint ) );
          }
        ++weightIterator;
        }
      } // end is a tube
    }

  if( weightSum.GetSum() == NumericTraits< ScalarType >::Zero )
    {
    itkWarningMacro(
      << "GetValue: All the transformed tube points are outside the image." );
    matchMeasure = NumericTraits< ScalarType >::min();
    }
  else
    {
    const ScalarType normalizedMeasure =
      static_cast< ScalarType >( matchMeasure.GetSum() / weightSum.GetSum() );
    matchMeasure = normalizedMeasure;
    }

  itkDebugMacro( << "matchMeasure = " << matchMeasure );

  delete tubeList;
  return matchMeasure.GetSum();
}


/** Compute the Laplacian magnitude */
// TODO FACTORIZE CODE --> See ComputeThirdDerivative
template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >::ScalarType
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeLaplacianMagnitude(
  const typename TubePointType::CovariantVectorType & tubeNormal,
  const ScalarType scale,
  const OutputPointType & currentPoint ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // currentPoint with a second derivative of a Gaussian
  const ScalarType scaleSquared = scale * scale;
  const ScalarType scaleExtentProduct = scale * m_Extent;
  CompensatedSummationType kernelSum;
  SizeValueType numberOfKernelPoints = 0;

  for( ScalarType distance = -scaleExtentProduct;
       distance <= scaleExtentProduct;
       //! \todo better calculation of the increment instead of just +1
       ++distance )
    {
    typename FixedImageType::PointType point;
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      point[ii] = currentPoint[ii] + distance * tubeNormal.GetElement(ii);
      }

    if( this->m_Interpolator->IsInsideBuffer( point ) )
      {
      const ScalarType distanceSquared = distance * distance;
      //! \todo cache this value, so it can be used in the next loop
      kernelSum += ( -1.0 + ( distanceSquared / scaleSquared ) )
        * vcl_exp( -0.5 * distanceSquared / scaleSquared );
      ++numberOfKernelPoints;
      }
    }

  //! \todo check this normalization (where is the 1/(scale * sqrt( 2 pi ))
  //term?
  const ScalarType error = kernelSum.GetSum() / numberOfKernelPoints;
  CompensatedSummationType result;
  for( ScalarType distance = -scaleExtentProduct;
       distance <= scaleExtentProduct;
       //! \todo better calculation of the increment instead of just +1
       ++distance )
    {
    const ScalarType distanceSquared = distance * distance;
    const ScalarType kernelValue =
      ( -1.0 + ( distanceSquared / scaleSquared ) )
      * vcl_exp( -0.5 * distanceSquared / scaleSquared ) - error;

    typename FixedImageType::PointType point;
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      point[ii] = currentPoint[ii] + distance * tubeNormal.GetElement(ii);
      }

    if( this->m_Interpolator->IsInsideBuffer( point ) )
      {
      const ScalarType value =
        static_cast< ScalarType >(
          this->m_Interpolator->Evaluate( point ) );
      result += value * kernelValue;
      }
    }

  return result.GetSum();
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetDeltaAngles( const OutputPointType & tubePoint,
  const VnlVectorType & dx,
  const VectorType & offsets,
  ScalarType dAngle[TubeDimension] ) const
{
  VectorType radius = ( tubePoint - offsets ) - ( this->m_CenterOfRotation );
  radius.Normalize();

  dAngle[0] = dx[1] * -radius[2] + dx[2] *  radius[1];
  dAngle[1] = dx[0] *  radius[2] + dx[2] * -radius[0];
  dAngle[2] = dx[0] * -radius[1] + dx[1] *  radius[0];
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  derivative.SetSize( this->GetNumberOfParameters() );

  // Create a copy of the transform to keep true const correctness (thread-safe)
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform =
    this->m_Transform->CreateAnother();
  TransformType * transformCopy =
    static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  itkDebugMacro( << "**** Get Derivative ****" );
  itkDebugMacro( << "parameters = "<< parameters )

  VnlMatrixType biasV( TubeDimension, TubeDimension,
    NumericTraits< ScalarType >::Zero );
  VnlMatrixType biasVI( TubeDimension, TubeDimension,
    NumericTraits< ScalarType >::Zero );

  ResolutionWeightsContainerType::const_iterator weightIterator;
  weightIterator = m_ResolutionWeights.begin();

  CompensatedSummationType dPosition[TubeDimension];

  VnlMatrixType tM( TubeDimension, TubeDimension );
  VnlVectorType v1T( TubeDimension );
  VnlVectorType v2T( TubeDimension );

  typedef std::vector< OutputPointType > TubePointsContainerType;
  TubePointsContainerType transformedTubePoints;
  transformedTubePoints.reserve( this->m_ResolutionWeights.size() );

  typedef std::vector< VectorType >      dTubePointsContainerType;
  dTubePointsContainerType dtransformedTubePoints;
  dtransformedTubePoints.reserve( this->m_ResolutionWeights.size() );

  derivative.fill(0.0);

  typename TubeNetType::ChildrenListType * tubeList = GetTubes();
  typename TubeNetType::ChildrenListType::const_iterator tubeIterator;
  for( tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeType* currentTube = static_cast<TubeType*>(
      ( *tubeIterator ).GetPointer() );

    typename TubeType::PointListType::const_iterator pointIterator;
    for( pointIterator = currentTube->GetPoints().begin();
         pointIterator != currentTube->GetPoints().end();
         ++pointIterator )
      {
      InputPointType inputPoint = pointIterator->GetPosition();
      OutputPointType currentPoint;
      if( this->IsInside( inputPoint, currentPoint, transformCopy ) )
        {
        transformedTubePoints.push_back( currentPoint );

        //! \todo: these should be CovariantVectors?
        VectorType v1;
        VectorType v2;
        for( unsigned int ii = 0; ii < TubeDimension; ++ii )
          {
          v1[ii] = pointIterator->GetNormal1()[ii];
          v2[ii] = pointIterator->GetNormal2()[ii];
          }

        for( unsigned int ii = 0; ii < TubeDimension; ++ii )
          {
          v1T[ii] = transformCopy->TransformVector( v1 )[ii];
          v2T[ii] = transformCopy->TransformVector( v2 )[ii];
          }
        tM = outer_product( v1T, v1T );
        tM = tM + outer_product( v2T, v2T );
        tM = *weightIterator * tM;
        biasV += tM;

        v1 = transformCopy->TransformVector( v1 );
        v2 = transformCopy->TransformVector( v2 );

        ScalarType scalingRadius = pointIterator->GetRadius();
        scalingRadius = std::max( scalingRadius, 0.5 );

        const ScalarType scale = scalingRadius * m_Kappa;

        const ScalarType dXProj1 = this->ComputeThirdDerivatives( v1, scale, currentPoint );
        const ScalarType dXProj2 = this->ComputeThirdDerivatives( v2, scale, currentPoint );

        VectorType dtransformedTubePoint;
        for( unsigned int ii = 0; ii < TubeDimension; ++ii )
          {
          dtransformedTubePoint[ii] = ( dXProj1 * v1[ii] + dXProj2 * v2[ii] );
          dPosition[ii] += *weightIterator * ( dtransformedTubePoint[ii] );
          }
        dtransformedTubePoints.push_back( dtransformedTubePoint );
        }
      ++weightIterator;
      }
    }
  delete tubeList;

  biasVI = vnl_matrix_inverse< ScalarType >( biasV ).inverse();

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
  weightIterator = m_ResolutionWeights.begin();
  typename dTubePointsContainerType::const_iterator  dtransformedTubePointsIt =
    dtransformedTubePoints.begin();
  while( dtransformedTubePointsIt != dtransformedTubePoints.end() )
    {
    VnlVectorType dXT( TubeDimension );
    for( unsigned int ii = 0; ii < TubeDimension; ++ii )
      {
      dXT[ii] = ( *dtransformedTubePointsIt )[ii];
      }

    dXT = dXT * biasVI;

    ScalarType angleDelta[TubeDimension];
    this->GetDeltaAngles( *transformedTubePointsIt, dXT, offsets, angleDelta );
    for( unsigned int ii = 0; ii < TubeDimension; ++ii )
      {
      dAngle[ii] += *weightIterator * angleDelta[ii];
      }
    ++weightIterator;
    ++dtransformedTubePointsIt;
    ++transformedTubePointsIt;
    }

  derivative[0] = dAngle[0].GetSum();
  derivative[1] = dAngle[1].GetSum();
  derivative[2] = dAngle[2].GetSum();
  derivative[3] = dPosition[0].GetSum();
  derivative[4] = dPosition[1].GetSum();
  derivative[5] = dPosition[2].GetSum();
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType & value,
                         DerivativeType & derivative ) const
{
  value = this->GetValue( parameters );
  this->GetDerivative( parameters, derivative );
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject,
  TTubeSpatialObject, TResolutionWeightFunction >::ScalarType
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeThirdDerivatives(
  const VectorType & tubeNormal,
  const ScalarType scale,
  const OutputPointType & currentPoint ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // currentPoint with a second derivative of a Gaussian
  CompensatedSummationType result;
  CompensatedSummationType kernelSum;

  const ScalarType scaleSquared = scale * scale;
  const ScalarType scaleExtentProduct = scale * m_Extent;

  for( ScalarType distance = -scaleExtentProduct;
       distance <= scaleExtentProduct;
       //! \todo better way to calculate this increment
       distance += 0.1 )
    {
    const ScalarType distanceSquared = distance * distance;
    const ScalarType kernelValue =
      2.0 * distance * vcl_exp( -0.5 * distanceSquared / scaleSquared );

    kernelSum += vnl_math_abs( kernelValue );

    typename FixedImageType::PointType point;
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      point[ii] = currentPoint[ii] + distance * tubeNormal.GetElement(ii);
      }

    if( this->m_Interpolator->IsInsideBuffer( point ) )
      {
      const ScalarType value =
        static_cast< ScalarType >(
          this->m_Interpolator->Evaluate( point ) );
      result += value * kernelValue;
      }
    }

  return result.GetSum() / kernelSum.GetSum();
}


template< class TFixedImage, class TMovingSpatialObject,
          class TTubeSpatialObject, class TResolutionWeightFunction >
bool
ImageToTubeRigidMetric< TFixedImage, TMovingSpatialObject, TTubeSpatialObject,
  TResolutionWeightFunction >
::IsInside( const InputPointType & inputPoint,
  OutputPointType & outputPoint,
  const TransformType * transform ) const
{
  outputPoint = transform->TransformPoint( inputPoint );
  return ( this->m_Interpolator->IsInsideBuffer( outputPoint ) );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeImageToTubeRigidMetric_hxx)
