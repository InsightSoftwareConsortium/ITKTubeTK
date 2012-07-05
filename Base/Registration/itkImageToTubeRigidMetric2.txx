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
#ifndef __itkImageToTubeRigidMetric2_txx
#define __itkImageToTubeRigidMetric2_txx

#include "itkImageToTubeRigidMetric2.h"

namespace itk
{

/** Constructor */
template < class TFixedImage, class TMovingSpatialObject>
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::ImageToTubeRigidMetric2()
{
  m_Iteration = 1;
  m_Kappa = 1;
  m_Sampling = 30;
  m_CachedValue = 0;

  m_Factors.fill(1.0);
  m_Offset.fill( std::numeric_limits<double>::min() );
  m_RotationCenter.fill(0);

  m_Extent = 3;
  m_Verbose = true;
  m_DerivativeImageFunction = DerivativeImageFunctionType::New();

  this->m_FixedImage = 0;           // has to be provided by the user.
  this->m_MovingSpatialObject = 0;  // has to be provided by the user.
  this->m_Transform = 0;            // has to be provided by the user.
  this->m_Interpolator = 0;         // has to be provided by the user.
}

/** SetImageRange */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::ComputeImageRange( void )
{
  m_RangeCalculator = RangeCalculatorType::New();
  m_RangeCalculator->SetImage( this->m_FixedImage );
  m_RangeCalculator->Compute();
  m_ImageMin = m_RangeCalculator->GetMinimum();
  m_ImageMax = m_RangeCalculator->GetMaximum();

  if( m_Verbose )
    {
    std::cout << "ImageMin = " << m_ImageMin << std::endl;
    std::cout << "ImageMax = " << m_ImageMax << std::endl;
    }
}

/** Initialize the metric  */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::Initialize( void ) throw ( ExceptionObject )
{
  m_SumWeight = 0;
  m_NumberOfPoints = 0;

  if( !this->m_MovingSpatialObject || !this->m_FixedImage)
    {
    std::cout << "SubSampleTube : No tube/image net plugged in ! " << std::endl;
    return;
    }

  // Get the number of points of the SpatialObject to reserve the right space.
  unsigned int totalNumberOfPoints = 0;
  TubeNetType::ChildrenListType* tubeList = GetTubes();
  TubeNetType::ChildrenListType::iterator tubeIterator = tubeList->begin();
  for ( ; tubeIterator != tubeList->end(); ++tubeIterator )
    {
    TubeType* currTube =
      static_cast<TubeType*>( ( *tubeIterator ).GetPointer() );

    totalNumberOfPoints += currTube->GetPoints().size();
    }
  m_TransformedPoints.reserve( totalNumberOfPoints / m_Sampling );

  this->ComputeImageRange();
  this->SubSampleTube();
  this->ComputeCenterRotation();

  this->m_Interpolator->SetInputImage( this->m_FixedImage );
  m_DerivativeImageFunction->SetInputImage( this->m_FixedImage );
}

/** Subsample the MovingSpatialObject tubenet  */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::SubSampleTube()
{
  double weight = 0.0;
  unsigned int tubeSize = 0;
  unsigned int step = m_Sampling / 2 - 1;
  TransformedPointType newPoint;

  TubeNetType::Pointer newTubeNet = TubeNetType::New();
  TubeNetType::ChildrenListType* tubeList = GetTubes();
  TubeNetType::ChildrenListType::iterator tubeIterator = tubeList->begin();
  for ( ; tubeIterator != tubeList->end(); ++tubeIterator )
    {
    TubeType::Pointer newTube = TubeType::New();
    TubeType* currTube =
      static_cast<TubeType*>( ( *tubeIterator ).GetPointer() );

    currTube->RemoveDuplicatePoints();
    currTube->ComputeTangentAndNormals();

    tubeSize = currTube->GetPoints().size();
    if ( tubeSize > m_Sampling )
      {
      typename std::vector<TubePointType>::iterator  tubePointIterator;
      int loopIndex = 0;
      unsigned int skippedPoints = 0;
      for ( tubePointIterator = currTube->GetPoints().begin();
            tubePointIterator != currTube->GetPoints().end();
            tubePointIterator++, ++loopIndex )
        {
        if( m_Sampling != 1 )
          {
          if ( std::distance(tubePointIterator, currTube->GetPoints().end() -1 )
               > step )
            {
            tubePointIterator +=
              (tubePointIterator == currTube->GetPoints().begin()) ? 0 : step;
            skippedPoints = ( loopIndex * ( step + 1 ) * 2 ) + 1;
            }
          else
            {
            skippedPoints = step + 2;
            tubePointIterator = currTube->GetPoints().end();
            }
          }

        // TODO Why +10 ?!
        if( tubePointIterator != currTube->GetPoints().end()
            && ( ( skippedPoints + 10 ) < tubeSize ) )
          {
          newTube->GetPoints().push_back( *( tubePointIterator ) );
          double val = -2 * ( tubePointIterator->GetRadius() );
          weight = 2.0 / ( 1.0 + exp( val ) );

          for( unsigned int i = 0; i < ImageDimension; ++i )
            {
            m_RotationCenter[i] +=
              weight * ( tubePointIterator->GetPosition() )[i];
            }
          m_SumWeight += weight;

          newPoint.TransformedPoint = *tubePointIterator;
          newPoint.Weight = weight;
          newPoint.Scale = std::max( static_cast<double>(
            tubePointIterator->GetRadius() ), 0.5 ) * m_Kappa;
          m_TransformedPoints.push_back( newPoint );

          m_NumberOfPoints++;
          }

        if( m_Sampling > 1 )
          {
          if ( std::distance(tubePointIterator, currTube->GetPoints().end() -1 )
               > step )
            {
            tubePointIterator += step;
            }
          else
            {
            tubePointIterator = currTube->GetPoints().end() - 1;
            }
          }
        }
      }
    newTubeNet->AddSpatialObject( newTube );
    }

  if ( m_Verbose )
    {
    std::cout << "Number of Points for the metric = "
              << m_NumberOfPoints << std::endl;
    }

  this->SetMovingSpatialObject( newTubeNet );
  delete tubeList;
}

/** Compute the Center of Rotation */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::ComputeCenterRotation()
{
  m_RotationCenter /= m_SumWeight;

  if( m_Verbose )
    {
    std::cout << "Center of Rotation = "
              << m_RotationCenter[0] << "  " \
              << m_RotationCenter[1] << "  " \
              << m_RotationCenter[2] << std::endl;
    std::cout << "Extent = " << m_Extent << std::endl;
    }
}

/** Get tubes contained within the Spatial Object */
// WARNING:
// Method might use GetMaximumDepth from ITKv4.
template < class TFixedImage, class TMovingSpatialObject>
typename GroupSpatialObject<3>::ChildrenListType*
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::GetTubes() const
{
  if (!this->m_MovingSpatialObject)
    {
    return 0;
    }

  char childName[] = "Tube";
  return this->m_MovingSpatialObject->GetChildren( 999999, childName );
  //return this->m_MovingSpatialObject
  // ->GetChildren( this->m_MovingSpatialObject->GetMaximumDepth(), childName );
}

/** Get the match Measure */
template < class TFixedImage, class TMovingSpatialObject>
typename ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>::MeasureType
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::GetValue( const ParametersType& parameters ) const
{
  if( m_Verbose )
    {
    std::cout << "**** Get Value ****" << std::endl;
    std::cout << "Parameters = " << parameters << std::endl;
    }

  this->UpdateTransformedPoints( parameters );

  MeasureType matchMeasure = 0;
  typename std::vector<TransformedPointType>::iterator transformedPointIterator;
  for ( transformedPointIterator = m_TransformedPoints.begin();
        transformedPointIterator != m_TransformedPoints.end();
        ++transformedPointIterator )
    {
    matchMeasure += transformedPointIterator->Weight
      * fabs( ComputeLaplacianMagnitude( *transformedPointIterator ) );
    }

  matchMeasure -= m_ImageMax * (m_NumberOfPoints - m_TransformedPoints.size());

  if( m_SumWeight == 0 )
    {
    std::cout << "GetValue: All the mapped image is outside ! " << std::endl;
    matchMeasure = -1;
    }
  else
    {
    float imageRatio = static_cast<float>( m_ImageMin / m_ImageMax );
    float normalizedMeasure = static_cast<float>( matchMeasure / m_SumWeight );
    matchMeasure = normalizedMeasure
                   - imageRatio
                   - static_cast<float>( m_ImageMin );
    }

  if( m_Verbose )
    {
    std::cout << "matchMeasure = " << matchMeasure << std::endl;
    }

  m_CachedValue = matchMeasure;
  return matchMeasure;
}

/** Compute the Laplacian magnitude */
template < class TFixedImage, class TMovingSpatialObject>
double
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::ComputeLaplacianMagnitude( TransformedPointType& point ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // m_CurrentPoint with a second derivative of a gaussian
  const double scaleSquared = std::pow( point.Scale, 2 );
  const double scaleExtentProduct = point.Scale * m_Extent;
  double result = 0;
  double wI = 0;

  itk::Index<3> index;
  for( double dist = -scaleExtentProduct; dist <= scaleExtentProduct; ++dist )
    {
    double distSquared = std::pow( dist, 2 );
    wI = exp( -0.5 * distSquared / scaleSquared );

    index[0] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[0]
          + dist * point.TransformedPoint.GetNormal1().GetElement(0) );
    index[1] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[1]
          + dist * point.TransformedPoint.GetNormal1().GetElement(1) );
    index[2] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[2]
          + dist * point.TransformedPoint.GetNormal1().GetElement(2) );

    itk::Size<3> size =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if( index[0] >= 0 && ( index[0] < static_cast<unsigned int>( size[0] ) )
      && index[1] >= 0 && ( index[1] < static_cast<unsigned int>( size[1] ) )
      && index[2] >= 0 && ( index[2] < static_cast<unsigned int>( size[2] ) ) )
      {
      double value = this->m_FixedImage->GetPixel( index );
      result += value * wI;
      }
    }

  return result;
}

/** GetDeltaAngles */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::GetDeltaAngles( const Point<double, 3> & x,
  const vnl_vector_fixed<double, 3> & dx,
  double angle[3] ) const
{
  vnl_vector_fixed<double, 3> tempV;
  vnl_vector_fixed<double, 3> pos;
  pos[0] = x[0];
  pos[1] = x[1];
  pos[2] = x[2];

  tempV = ( pos - m_Offset ) - m_RotationCenter;
  tempV.normalize();

  angle[0] = dx[1] * ( -tempV[2] ) + dx[2] * tempV[1];
  angle[1] = dx[0] * tempV[2] + dx[2] * ( -tempV[0] );
  angle[2] = dx[0] * ( -tempV[1] ) + dx[1] * tempV[0];
}

/** Set the offset */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::SetOffset( double oX, double oY, double oZ ) const
{
  if( m_Offset[0] == oX && m_Offset[1] == oY && m_Offset[2] == oZ )
    {
    return;
    }

  m_Offset[0] = oX;
  m_Offset[1] = oY;
  m_Offset[2] = oZ;
}

/** Get the Derivative Measure */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  if( m_Verbose )
    {
    std::cout << "**** Get Derivative ****" << std::endl;
    std::cout <<"parameters = "<< parameters << std::endl;
    }

  derivative = DerivativeType( SpaceDimension );
  this->UpdateTransformedPoints( parameters );

  vnl_matrix<double> tM( 3, 3 );        // Weighted normal-plane component
  vnl_matrix<double> biasV( 3, 3, 0 );  // Bias Matrix
  double dAngle[3] = { 0, 0, 0 };
  double dXProj1, dXProj2;
  vnl_vector<double> dPosition( 3 );

  Vector<double, 3> v1;         // Co-Vector that defines the direction normal1
  Vector<double, 3> v2;         // Co-Vector that defines the direction normal2
  vnl_vector<double> v1T( 3 );  // Co-Vector that defines the direction normal1
  vnl_vector<double> v2T( 3 );  // Co-Vector that defines the direction normal2

  PointsType  dXTlist;
  typename std::vector<TransformedPointType>::iterator transformedPointIterator;
  for ( transformedPointIterator = m_TransformedPoints.begin();
        transformedPointIterator != m_TransformedPoints.end();
        ++transformedPointIterator )
    {
    for( unsigned int i = 0; i < 3; ++i )
      {
      v1[i] = transformedPointIterator->TransformedPoint.GetNormal1()[i];
      v2[i] = transformedPointIterator->TransformedPoint.GetNormal2()[i];
      }

    v1 = this->m_Transform->TransformVector( v1 );
    v2 = this->m_Transform->TransformVector( v2 );

    dXProj1 = ComputeThirdDerivatives( &v1, *transformedPointIterator );
    dXProj2 = ComputeThirdDerivatives( &v2, *transformedPointIterator );

    PointType dXT;
    for( unsigned int i = 0; i < 3; ++i )
      {
      dXT[i] = ( dXProj1 * v1[i] + dXProj2 * v2[i] );
      v1T[i] = v1[i];
      v2T[i] = v2[i];
      }

    tM = outer_product( v1T, v1T ) + outer_product( v2T, v2T );
    tM *= transformedPointIterator->Weight;

    biasV += tM;

    dPosition[0] += transformedPointIterator->Weight * ( dXT[0] );
    dPosition[1] += transformedPointIterator->Weight * ( dXT[1] );
    dPosition[2] += transformedPointIterator->Weight * ( dXT[2] );
    dXTlist.push_back( dXT );
    }

  dPosition *= vnl_matrix_inverse<double>( biasV ).inverse();

  if( m_SumWeight == 0 )
    {
    biasV = 0;
    dAngle[0] = dAngle[1] = dAngle[2] = 0;
    derivative.Fill(0);
    m_CachedDerivative.Fill(0);

    if( m_Verbose )
      {
      std::cout << "GetDerivative : sumWeight == 0 !" << std::endl;
      }
    return;
    }
  else
    {
    biasV = 1.0 / m_SumWeight * biasV;
    }

  this->SetOffset( dPosition[0], dPosition[1], dPosition[2] );
  PointsType::iterator  dXTIterator = dXTlist.begin();

  for ( transformedPointIterator = m_TransformedPoints.begin();
        transformedPointIterator != m_TransformedPoints.end();
        ++transformedPointIterator, ++dXTIterator )
    {
    vnl_vector<double> dXT( 3 );
    dXT[0] = ( *dXTIterator )[0];
    dXT[1] = ( *dXTIterator )[1];
    dXT[2] = ( *dXTIterator )[2];

    dXT = dXT * vnl_matrix_inverse<double>( biasV ).inverse();
    const Point<double, 3> & xT =
      transformedPointIterator->TransformedPoint.GetPosition();

    double angleDelta[3];
    GetDeltaAngles( xT, dXT, angleDelta );
    dAngle[0] += transformedPointIterator->Weight * angleDelta[0];
    dAngle[1] += transformedPointIterator->Weight * angleDelta[1];
    dAngle[2] += transformedPointIterator->Weight * angleDelta[2];
    }

  dAngle[0] /= m_SumWeight * dXTlist.size();
  dAngle[1] /= m_SumWeight * dXTlist.size();
  dAngle[2] /= m_SumWeight * dXTlist.size();

  if( m_Verbose )
    {
    std::cout << "dA = " << dAngle[0] << std::endl;
    std::cout << "dB = " << dAngle[1] << std::endl;
    std::cout << "dG = " << dAngle[2] << std::endl;
    std::cout << "dX = " << dPosition[0] << std::endl;
    std::cout << "dY = " << dPosition[1] << std::endl;
    std::cout << "dZ = " << dPosition[2] << std::endl;
    }

  if( m_Iteration > 0 )
    {
    derivative[0] = dAngle[0];
    derivative[1] = dAngle[1];
    derivative[2] = dAngle[2];
    derivative[3] = dPosition[0];
    derivative[4] = dPosition[1];
    derivative[5] = dPosition[2];
    m_CachedDerivative = derivative;
    }
}

/** Get both the match Measure and theDerivative Measure */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType&  Value,
                         DerivativeType& Derivative ) const
{
  if ( !this->UpdateTransformedPoints( parameters ) )
    {
    Value = m_CachedValue;
    Derivative = m_CachedDerivative;
    return;
    }

  Value = GetValue( parameters );
  GetDerivative( parameters, Derivative );
}

template < class TFixedImage, class TMovingSpatialObject>
double
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::ComputeThirdDerivatives( Vector<double, 3> *v,
                           TransformedPointType& point ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // m_CurrentPoint with a second derivative of a gaussian
  double result = 0;
  double wI = 0;
  itk::Index<3> index;

  double wTotalX = 0;
  const double scaleSquared = std::pow( point.Scale, 2 );
  const double scaleExtentProduct = point.Scale * m_Extent;

  for( double dist = -scaleExtentProduct;
       dist <= scaleExtentProduct;
       dist += 0.1 )
    {
    double distSquared = std::pow( dist, 2 );
    wI = - dist / scaleSquared * exp( -0.5 * distSquared / scaleSquared );
    wTotalX += fabs( wI );

    index[0] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[0] + dist * v->GetElement(0) );
    index[1] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[1] + dist * v->GetElement(1) );
    index[2] = static_cast<unsigned int>(
          point.TransformedPoint.GetPosition()[2] + dist * v->GetElement(2) );

    itk::Size<3> size =
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

/** Update the transformed points location given the transformation parameters
 * \warning This method cannot be safely used in more than one thread at
 * a time.  */
template < class TFixedImage, class TMovingSpatialObject>
bool
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::UpdateTransformedPoints( const ParametersType& parameters ) const
{
  if ( m_Offset[0] == parameters[3] &&
       m_Offset[1] == parameters[4] &&
       m_Offset[2] == parameters[5] &&
       parameters == GetTransform()->GetParameters() )
    {
    return false;
    }

  SetOffset( parameters[3], parameters[4], parameters[5] );
  this->m_Transform->SetParameters( parameters );

  this->UpdateTransformedPoints();
  return true;
}

/** Update the transformed points location given the transformation parameters
 * \warning This method cannot be safely used in more than one thread at
 * a time.  */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric2<TFixedImage, TMovingSpatialObject>
::UpdateTransformedPoints() const
{
  m_SumWeight = 0;
  unsigned int indexInsideTransformedPoint = 0;
  TubeNetType::ChildrenListType* tubeList = GetTubes();
  TubeNetType::ChildrenListType::iterator tubeIterator = tubeList->begin();
  for( ; tubeIterator != tubeList->end(); ++tubeIterator )
    {
    TubeType* currTube = static_cast<TubeType*>(
      ( *tubeIterator ).GetPointer() );

    std::vector<TubePointType>::iterator pointIterator;
    for( pointIterator = currTube->GetPoints().begin();
         pointIterator != currTube->GetPoints().end();
         ++pointIterator )
      {
      OutputPointType tempPoint;
      Matrix<double, 3, 3> matrix = GetTransform()->GetRotationMatrix();
      vnl_vector<double> rotationOffset = matrix * m_RotationCenter;

      tempPoint = matrix * pointIterator->GetPosition()
                  + GetTransform()->GetOffset();

      tempPoint[0] += m_RotationCenter[0] - rotationOffset[0];
      tempPoint[1] += m_RotationCenter[1] - rotationOffset[1];
      tempPoint[2] += m_RotationCenter[2] - rotationOffset[2];

      if( this->m_Interpolator->IsInsideBuffer( tempPoint ) )
        {
        // Lazy resizment
        if (indexInsideTransformedPoint > m_TransformedPoints.size() + 1)
          {
          m_TransformedPoints.resize( m_NumberOfPoints );
          }

        m_TransformedPoints[indexInsideTransformedPoint++].
            TransformedPoint.SetPosition( tempPoint );
        m_SumWeight += m_TransformedPoints[indexInsideTransformedPoint].Weight;
        }
      }
    }

  m_TransformedPoints.resize( indexInsideTransformedPoint );
}

} // end namespace itk

#endif
