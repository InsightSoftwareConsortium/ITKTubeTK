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

#ifndef __itktubeMeanSquareRegistrationFunction_hxx
#define __itktubeMeanSquareRegistrationFunction_hxx


#include <vnl/vnl_math.h>

namespace itk
{

namespace tube
{

/**
 * Default constructor
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
MeanSquareRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::MeanSquareRegistrationFunction( void )
{

  RadiusType r;
  unsigned int j;
  for( j = 0; j < ImageDimension; j++ )
    {
    r[j] = 0;
    }
  this->SetRadius( r );

  this->SetEnergy( 0.0 );
  m_TimeStep = 1.0;
  m_DenominatorThreshold = 1e-9;
  m_IntensityDifferenceThreshold = 0.001;
  this->SetMovingImage( NULL );
  this->SetFixedImage( NULL );
  m_FixedImageGradientCalculator = GradientCalculatorType::New();


  typename DefaultInterpolatorType::Pointer interp =
    DefaultInterpolatorType::New();

  m_MovingImageInterpolator = static_cast<InterpolatorType*>(
    interp.GetPointer() );

  m_BackgroundIntensity = 0.0;

}


/*
 * Standard "PrintSelf" method.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
MeanSquareRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
/*
  os << indent << "MovingImageIterpolator: ";
  os << m_MovingImageInterpolator.GetPointer() << std::endl;
  os << indent << "FixedImageGradientCalculator: ";
  os << m_FixedImageGradientCalculator.GetPointer() << std::endl;
  os << indent << "DenominatorThreshold: ";
  os << m_DenominatorThreshold << std::endl;
  os << indent << "IntensityDifferenceThreshold: ";
  os << m_IntensityDifferenceThreshold << std::endl;
*/
}


/*
 * Set the function state values before each iteration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
MeanSquareRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::InitializeIteration( void )
{
  if( !this->GetMovingImage() || !this->GetFixedImage()
    || !m_MovingImageInterpolator )
    {
    itkExceptionMacro( << "MovingImage, FixedImage and/or Interpolator not "
                       << "set." );
    }

  // cache fixed image information
  m_FixedImageSpacing    = this->GetFixedImage()->GetSpacing();

  // setup gradient calculator
  m_FixedImageGradientCalculator->SetInputImage( this->GetFixedImage() );

  // setup moving image interpolator
  m_MovingImageInterpolator->SetInputImage( this->GetMovingImage() );

  this->SetEnergy( 0.0 );
}


/**
 * Compute update at a non boundary neighborhood
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename MeanSquareRegistrationFunction<TFixedImage, TMovingImage,
  TDeformationField>::PixelType
MeanSquareRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeUpdate( const NeighborhoodType &it, void * itkNotUsed(
    globalData ), const FloatOffsetType& itkNotUsed( offset ) )
{
  // Get fixed image related information
  // Note: no need to check the index is within
  // fixed image buffer. This is done by the external filter.
  const IndexType index = it.GetIndex();
  const CovariantVectorType fixedGradient =
    m_FixedImageGradientCalculator->EvaluateAtIndex( index );
  double fixedGradientSquaredMagnitude = 0.0;
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    fixedGradientSquaredMagnitude +=
      vnl_math::sqr( fixedGradient[j] ) * m_FixedImageSpacing[j];
    }

  // Compute update
  DeformationFieldPixelType itvec = this->GetDisplacementField()->GetPixel(
    index );
  const double speedValue = this->ComputeIntensityDifference( index, itvec );

  const bool normalizemetric=this->GetNormalizeGradient();
  double denominator = 1.0;
  if( normalizemetric )
    {
    denominator = speedValue*speedValue *fixedGradientSquaredMagnitude;
    denominator = std::sqrt( denominator );
    }
  if( denominator == 0 )
    {
    denominator=1.0;
    }
  PixelType update;
  if( std::fabs( speedValue ) < m_IntensityDifferenceThreshold ||
    denominator < m_DenominatorThreshold )
    {
    update.Fill( 0.0 );
    return update;
    }

  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
    update[j] = speedValue * fixedGradient[j] /
      denominator*this->m_GradientStep;
    if( normalizemetric )
      {
      update[j] *= vnl_math::sqr( m_FixedImageSpacing[j] );
      }
    }
  return update;
}

/**
 * Compute intensity difference between fixed and moving images at a
 * non boundary neighborhood
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
double
MeanSquareRegistrationFunction<TFixedImage, TMovingImage, TDeformationField>
::ComputeIntensityDifference( const IndexType & index,
                             const DeformationFieldPixelType & itvec )
{
  const double fixedValue = ( double ) this->GetFixedImage()->GetPixel(
    index );

  // Get moving image related information
  PointType mappedPoint;
  this->GetFixedImage()->TransformIndexToPhysicalPoint( index, mappedPoint );
  for( unsigned int j = 0; j < ImageDimension; j++ )
    {
      mappedPoint[j] += itvec[j];
    }

  double movingValue = fixedValue;
  if( m_MovingImageInterpolator->IsInsideBuffer( mappedPoint ) )
    {
    movingValue = m_MovingImageInterpolator->Evaluate( mappedPoint );
    }

  // Compute update
  double difference = fixedValue - movingValue;

  return difference;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMeanSquareRegistrationFunction_hxx )
