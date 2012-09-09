/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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
#ifndef __itkImageToTubeRigidMetric_txx
#define __itkImageToTubeRigidMetric_txx

#include "itkImageToTubeRigidMetric.h"

namespace itk
{

template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::ImageToTubeRigidMetric()
{
  m_Iteration = 1;
  m_Kappa = 1;
  m_Extent = 3;     // TODO Check depedencies --> enum { ImageDimension = 3 };

  m_InitialScale = 1;
  m_Factors.fill(1.0);
  m_CenterOfRotation.Fill( 0.0 );

  m_DerivativeImageFunction = DerivativeImageFunctionType::New();

  this->m_FixedImage = 0;           // has to be provided by the user.
  this->m_MovingSpatialObject = 0;  // has to be provided by the user.
  this->m_Transform = 0;            // has to be provided by the user.
  this->m_Interpolator = 0;         // has to be provided by the user.
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::~ImageToTubeRigidMetric()
{
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >::ResolutionWeightFunctionType &
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetResolutionWeightFunction()
{
  return this->m_ResolutionWeightFunction;
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
const typename ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >::ResolutionWeightFunctionType &
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetResolutionWeightFunction() const
{
  return this->m_ResolutionWeightFunction;
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::SetResolutionWeightFunction( const ResolutionWeightFunctionType & function )
{
  if( this->m_ResolutionWeightFunction != function )
    {
    this->m_ResolutionWeightFunction = function;
    this->Modified();
    }
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::Initialize( void ) throw ( ExceptionObject )
{
  m_ResolutionWeights.clear();

  if( !this->m_MovingSpatialObject || !this->m_FixedImage )
    {
    itkExceptionMacro( "No tube/image net plugged in !" );
    }

  this->ComputeImageRange();
  this->ComputeTubePointResolutionWeights();

  this->m_Interpolator->SetInputImage( this->m_FixedImage );
  this->m_DerivativeImageFunction->SetInputImage( this->m_FixedImage );
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeImageRange( void )
{
  m_RangeCalculator = RangeCalculatorType::New();
  m_RangeCalculator->SetImage( this->m_FixedImage );
  m_RangeCalculator->Compute();
  m_ImageMin = m_RangeCalculator->GetMinimum();
  m_ImageMax = m_RangeCalculator->GetMaximum();

  itkDebugMacro( "ImageMin = " << m_ImageMin );
  itkDebugMacro( "ImageMax = " << m_ImageMax );
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeTubePointResolutionWeights()
{
  typename TubeNetType::ChildrenListType* tubeList = this->GetTubes();

  //! \todo Use CompensatedSummation from ITKv4.
  InternalComputationValueType resolutionWeightSum =
    NumericTraits< InternalComputationValueType >::Zero;

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
      for ( TubePointIteratorType tubePointIterator = currentTubePoints.begin();
            tubePointIterator != currentTubePoints.end();
            ++tubePointIterator )
        {
        const InternalComputationValueType weight =
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

  for ( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    this->m_CenterOfRotation[ii] /= resolutionWeightSum;
    }

  //! \todo partial template specialization for transforms with a center?
  TransformType * transform = static_cast< TransformType * >
    ( this->m_Transform.GetPointer() );
  transform->SetCenter( this->m_CenterOfRotation );
}


/** Get tubes contained within the Spatial Object */
// WARNING:
// Method might use GetMaximumDepth from ITK.
// Patch pushed in ITKv4, waiting for validation.
template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >::TubeNetType::ChildrenListType*
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetTubes() const
{
  if (!this->m_MovingSpatialObject)
    {
    return NULL;
    }

  char childName[] = "Tube";
  return this->m_MovingSpatialObject->GetChildren( 999999, childName );
  //return this->m_MovingSpatialObject
  // ->GetChildren( this->m_MovingSpatialObject->GetMaximumDepth(), childName );
}


/** Get the match Measure */
// TODO Do not pass the parameter as arguments use instead
// the transform parameters previously set.
template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
typename ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >::MeasureType
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro( "**** Get Value ****" );
  itkDebugMacro( "Parameters = " << parameters );

  MeasureType matchMeasure = 0.0;
  // \todo replace with ITKv4 CompensatedSummation
  InternalComputationValueType sumWeight =
    NumericTraits< InternalComputationValueType >::Zero;
  InternalComputationValueType scale = this->m_InitialScale;

  ResolutionWeightsContainerType::const_iterator weightIterator;
  weightIterator = this->m_ResolutionWeights.begin();

  // Create a copy of the transform to keep true const correctness (thread-safe)
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy = static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  typename TubeNetType::ChildrenListType::iterator tubeIterator;
  typename TubeNetType::ChildrenListType * tubeList = GetTubes();
  for( tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       tubeIterator++ )
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
            sumWeight += *weightIterator;
            InternalComputationValueType scalingRadius = pointIterator->GetRadius();
            // !TODO 0.5 should be a parameter of the class
            scalingRadius = std::max( scalingRadius, 0.5 );

            scale = scalingRadius * m_Kappa;

            Vector<InternalComputationValueType, TubeDimension> v2;
            for( unsigned int ii = 0; ii < TubeDimension; ++ii )
              {
              v2[ii] = pointIterator->GetNormal1()[ii];
              }

              matchMeasure += *weightIterator * fabs(
                ComputeLaplacianMagnitude( &v2, scale, currentPoint ) );
              }
          else
            {
            matchMeasure -= m_ImageMax;
            }

          weightIterator++;
          }
      } // end is a tube
    }

  if( sumWeight == 0 )
    {
    std::cerr << "GetValue: All the mapped image is outside ! " << std::endl;
    matchMeasure = -1;
    }
  else
    {
    float imageRatio = static_cast<float>( m_ImageMin / m_ImageMax );
    float normalizedMeasure = static_cast<float>( matchMeasure / sumWeight );
    matchMeasure = normalizedMeasure
                   - imageRatio
                   - static_cast<float>( m_ImageMin );
    }

  itkDebugMacro( "matchMeasure = " << matchMeasure );

  delete tubeList;
  return matchMeasure;
}


/** Compute the Laplacian magnitude */
// TODO FACTORIZE CODE --> See computeThirdDerivative
template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
double
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeLaplacianMagnitude( Vector< InternalComputationValueType, 3> *v,
  const InternalComputationValueType & scale,
  const OutputPointType & currentPoint ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // currentPoint with a second derivative of a gaussian
  const InternalComputationValueType scaleSquared = scale * scale;
  const InternalComputationValueType scaleExtentProduct = scale * m_Extent;
  InternalComputationValueType result = 0.0;
  InternalComputationValueType wI = 0.0;
  unsigned int n = 0;
  typename FixedImageType::IndexType index;

  for( double dist = -scaleExtentProduct; dist <= scaleExtentProduct; ++dist )
    {
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      index[ii] =
        static_cast< IndexValueType >( currentPoint[ii] + dist * v->GetElement(ii) );
      }

    double distSquared = dist * dist;
    const typename FixedImageType::SizeType size =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if( index[0] >= 0 && ( index[0] < static_cast< IndexValueType >( size[0] ) )
      && index[1] >= 0 && ( index[1] < static_cast< IndexValueType >( size[1] ) )
      && index[2] >= 0 && ( index[2] < static_cast< IndexValueType >( size[2] ) ) )
      {
      wI += ( -1 + ( distSquared / scaleSquared ) )
        * exp( -0.5 * distSquared / scaleSquared );
      n++;
      }
    }

  const double error = wI / n;
  for( double dist = -scaleExtentProduct; dist <= scaleExtentProduct; ++dist )
    {
    const double distSquared = dist * dist;
    wI = ( -1 + ( distSquared / scaleSquared ) )
      * exp( -0.5 * distSquared / scaleSquared ) - error;

    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      index[ii] =
        static_cast< IndexValueType >( currentPoint[ii] + dist * v->GetElement(ii) );
      }

    const typename FixedImageType::SizeType size =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if( index[0] >= 0 && ( index[0] < static_cast< IndexValueType >( size[0] ) )
      && index[1] >= 0 && ( index[1] < static_cast< IndexValueType >( size[1] ) )
      && index[2] >= 0 && ( index[2] < static_cast< IndexValueType >( size[2] ) ) )
      {
      const double value = this->m_FixedImage->GetPixel( index );
      result += value * wI;
      }
    }

  return result;
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetDeltaAngles( const Point<double, 3> & x,
  const vnl_vector_fixed<double, 3> & dx,
  const vnl_vector_fixed<double, 3> & offsets,
  double angle[3] ) const
{
  vnl_vector_fixed<double, 3> tempV;
  vnl_vector_fixed<double, 3> pos;
  pos[0] = x[0];
  pos[1] = x[1];
  pos[2] = x[2];

  tempV = ( pos - offsets ) - ( this->m_CenterOfRotation.GetVnlVector() );
  tempV.normalize();

  angle[0] = dx[1] * ( -tempV[2] ) + dx[2] * tempV[1];
  angle[1] = dx[0] * tempV[2] + dx[2] * ( -tempV[0] );
  angle[2] = dx[0] * ( -tempV[1] ) + dx[1] * tempV[0];
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  derivative = DerivativeType( SpaceDimension );

  // Create a copy of the transform to keep true const correctness (thread-safe)
  // Set the parameters on the copy and uses the copy.
  LightObject::Pointer anotherTransform = this->m_Transform->CreateAnother();
  TransformType * transformCopy = static_cast< TransformType * >( anotherTransform.GetPointer() );
  transformCopy->SetFixedParameters( this->m_Transform->GetFixedParameters() );
  transformCopy->SetParameters( parameters );

  //! \todo remove me
  this->m_Transform->SetParameters( parameters );

  InternalComputationValueType scale = this->m_InitialScale;
  InternalComputationValueType opR;

  itkDebugMacro( "**** Get Derivative ****" );
  itkDebugMacro( "parameters = "<< parameters )

  vnl_matrix<double> biasV( 3, 3, 0 );
  vnl_matrix<double> biasVI( 3, 3, 0 );


  ResolutionWeightsContainerType::const_iterator weightIterator;
  weightIterator = m_ResolutionWeights.begin();

  double dPosition[3] = { 0, 0, 0 };
  double dAngle[3] = { 0, 0, 0 };
  double dXProj1, dXProj2;

  vnl_vector<double> tV( 3 );
  vnl_matrix<double> tM( 3, 3 );
  vnl_vector<double> v1T( 3 );
  vnl_vector<double> v2T( 3 );
  double angleDelta[3];

  //! \todo replace with ITKv4 CompensatedSummation
  InternalComputationValueType sumWeight =
    NumericTraits< InternalComputationValueType >::Zero;

  typedef itk::Vector<double, 3>    ITKVectorType;
  typedef std::list<ITKVectorType>  ListType;
  ListType dXTlist;
  FixedArray<Point<double, 3>, 5000> XTlist;

  unsigned int listindex = 0;
  typename TubeNetType::ChildrenListType* tubeList = GetTubes();
  typename TubeNetType::ChildrenListType::iterator tubeIterator = tubeList->begin();
  for( ; tubeIterator != tubeList->end(); ++tubeIterator )
    {
    TubeType* currentTube = static_cast<TubeType*>(
      ( *tubeIterator ).GetPointer() );

    typename std::vector<TubePointType>::iterator pointIterator;
    for( pointIterator = currentTube->GetPoints().begin();
         pointIterator != currentTube->GetPoints().end();
         ++pointIterator )
      {
      InputPointType inputPoint = pointIterator->GetPosition();
      InputPointType point;
      OutputPointType currentPoint;
      typename TransformType::MatrixType matrix =
        this->GetTransform()->GetMatrix();

      point =  matrix * inputPoint + this->GetTransform()->GetOffset();

      CenterOfRotationType rotationOffset = matrix * this->m_CenterOfRotation;

      for( unsigned int ii = 0; ii < TubeDimension; ++ii )
        {
        point[ii] += this->m_CenterOfRotation[ii] - rotationOffset[ii];
        }

      typename FixedImageType::IndexType index;
      for( unsigned int ii = 0; ii < TubeDimension; ++ii )
        {
        index[ii] = static_cast< IndexValueType >( point[ii] );
        }

      if( this->IsInside( inputPoint, currentPoint, transformCopy ) )
        {
        XTlist[listindex++] = currentPoint;
        sumWeight += *weightIterator;
        opR = pointIterator->GetRadius();
        opR = std::max( opR, 0.5 );

        scale = opR * m_Kappa;

        //! \todo: these should be CovariantVectors?
        Vector<double, 3> v1;
        Vector<double, 3> v2;
        for( unsigned int i = 0; i < 3; ++i )
          {
          v1[i] = pointIterator->GetNormal1()[i];
          v2[i] = pointIterator->GetNormal2()[i];
          }

        for( unsigned int i = 0; i < 3; ++i )
          {
          v1T[i] = this->m_Transform->TransformVector( v1 )[i];
          v2T[i] = this->m_Transform->TransformVector( v2 )[i];
          }

        v1 = this->m_Transform->TransformVector( v1 );
        v2 = this->m_Transform->TransformVector( v2 );

        dXProj1 = ComputeThirdDerivatives( &v1, scale, currentPoint );
        dXProj2 = ComputeThirdDerivatives( &v2, scale, currentPoint );

        Vector<double, 3> dXT;

        for( unsigned int i = 0; i < 3; ++i )
          {
          dXT[i] = ( dXProj1 * v1[i] + dXProj2 * v2[i] );
          }

        tM = outer_product( v1T, v1T );
        tM = tM + outer_product( v2T, v2T );

        tM = *weightIterator * tM;

        biasV += tM;

        dPosition[0] += *weightIterator * ( dXT[0] );
        dPosition[1] += *weightIterator * ( dXT[1] );
        dPosition[2] += *weightIterator * ( dXT[2] );

        dXTlist.push_back( dXT );
        }
      weightIterator++;
      }
    }

  biasVI = vnl_matrix_inverse<double>( biasV ).inverse();

  tV( 0 ) = dPosition[0];
  tV( 1 ) = dPosition[1];
  tV( 2 ) = dPosition[2];

  tV *= biasVI;

  dPosition[0] = tV( 0 );
  dPosition[1] = tV( 1 );
  dPosition[2] = tV( 2 );

  if( sumWeight == NumericTraits< InternalComputationValueType >::Zero )
    {
    biasV = 0;
    dAngle[0] = dAngle[1] = dAngle[2] = 0;
    derivative.fill(0);

    itkWarningMacro( "GetDerivative : sumWeight == 0 !" );
    return;
    }
  else
    {
    biasV = 1.0 / sumWeight * biasV;
    }

  biasVI = vnl_matrix_inverse<double>( biasV ).inverse();

  weightIterator = m_ResolutionWeights.begin();
  ListType::iterator  dXTIterator = dXTlist.begin();

  listindex  = 0;

  // ImageDimension correct here?
  vnl_vector_fixed< double, ImageDimension > offsets;
  offsets[0] = dPosition[0];
  offsets[1] = dPosition[1];
  offsets[2] = dPosition[2];

  while( dXTIterator != dXTlist.end() )
    {
    vnl_vector<double> dXT( 3 );
    dXT[0] = ( *dXTIterator )[0];
    dXT[1] = ( *dXTIterator )[1];
    dXT[2] = ( *dXTIterator )[2];

    dXT = dXT * biasVI;
    const Point<double, 3> & xT = XTlist[listindex++];

    this->GetDeltaAngles( xT, dXT, offsets, angleDelta );
    dAngle[0] += *weightIterator * angleDelta[0];
    dAngle[1] += *weightIterator * angleDelta[1];
    dAngle[2] += *weightIterator * angleDelta[2];
    weightIterator++;
    dXTIterator++;
    }

  dAngle[0] /= sumWeight * dXTlist.size();
  dAngle[1] /= sumWeight * dXTlist.size();
  dAngle[2] /= sumWeight * dXTlist.size();

  itkDebugMacro( "dA = " << dAngle[0] );
  itkDebugMacro( "dB = " << dAngle[1] );
  itkDebugMacro( "dG = " << dAngle[2] );
  itkDebugMacro( "dX = " << dPosition[0] );
  itkDebugMacro( "dY = " << dPosition[1] );
  itkDebugMacro( "dZ = " << dPosition[2] );

  if( m_Iteration > 0 )
    {
    derivative[0] = dAngle[0];
    derivative[1] = dAngle[1];
    derivative[2] = dAngle[2];
    derivative[3] = dPosition[0];
    derivative[4] = dPosition[1];
    derivative[5] = dPosition[2];
    }

  delete tubeList;
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
void
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType & value,
                         DerivativeType & derivative ) const
{
  value = this->GetValue( parameters );
  this->GetDerivative( parameters, derivative );
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
double
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::ComputeThirdDerivatives( Vector< InternalComputationValueType, 3> *v,
  const InternalComputationValueType & scale,
  const OutputPointType & currentPoint ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // currentPoint with a second derivative of a gaussian
  InternalComputationValueType result = 0.0;
  InternalComputationValueType wI = 0.0;

  InternalComputationValueType wTotalX = 0.0;
  const InternalComputationValueType scaleSquared = scale * scale;
  const InternalComputationValueType scaleExtentProduct = scale * m_Extent;

  for( double dist = -scaleExtentProduct;
       dist <= scaleExtentProduct;
       dist += 0.1 )
    {
    wI = 2 * dist * exp( -0.5 * std::pow( dist, 2 ) / scaleSquared );

    wTotalX += fabs( wI );

    typename FixedImageType::IndexType index;
    for( unsigned int ii = 0; ii < ImageDimension; ++ii )
      {
      index[ii] =
        static_cast< IndexValueType >( currentPoint[ii] + dist * v->GetElement(ii) );
      }

    const typename FixedImageType::SizeType size =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if( index[0] >= 0 && ( index[0] < static_cast<unsigned int>( size[0] ) )
      && index[1] >= 0 && ( index[1] < static_cast<unsigned int>( size[1] ) )
      && index[2] >= 0 && ( index[2] < static_cast<unsigned int>( size[2] ) ) )
      {
      result += this->m_FixedImage->GetPixel( index ) * wI;
      }
    }

  return result / wTotalX;
}


template < class TFixedImage,
  class TMovingSpatialObject,
  class TTubeSpatialObject,
  class TResolutionWeightFunction >
bool
ImageToTubeRigidMetric< TFixedImage,
  TMovingSpatialObject,
  TTubeSpatialObject,
  TResolutionWeightFunction >
::IsInside( const InputPointType & inputPoint,
  OutputPointType & outputPoint,
  const TransformType * transform ) const
{
  outputPoint = transform->TransformPoint( inputPoint );
  return ( this->m_Interpolator->IsInsideBuffer( outputPoint ) );
}

} // end namespace itk

#endif
