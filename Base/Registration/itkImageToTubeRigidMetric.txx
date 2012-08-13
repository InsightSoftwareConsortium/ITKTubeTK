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

template < class TFixedImage, class TMovingSpatialObject>
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ImageToTubeRigidMetric()
{
  m_Iteration = 1;
  m_Kappa = 1;
  m_Sampling = 30;

  m_InitialScale = 1;
  m_Factors.fill(1.0);
  m_CenterOfRotation.Fill( 0.0 );

  m_Extent = 3;     // TODO Check depedencies --> enum { ImageDimension = 3 };
  m_DerivativeImageFunction = DerivativeImageFunctionType::New();

  this->m_FixedImage = 0;           // has to be provided by the user.
  this->m_MovingSpatialObject = 0;  // has to be provided by the user.
  this->m_Transform = 0;            // has to be provided by the user.
  this->m_Interpolator = 0;         // has to be provided by the user.
}


template < class TFixedImage, class TMovingSpatialObject>
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::~ImageToTubeRigidMetric()
{
}


template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
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


template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::Initialize( void ) throw ( ExceptionObject )
{
  m_Weight.clear();
  m_NumberOfPoints = 0;
  m_SumWeight = 0.0;

  if( !this->m_MovingSpatialObject || !this->m_FixedImage )
    {
    itkExceptionMacro( "No tube/image net plugged in !" );
    }

  this->ComputeImageRange();
  this->SubSampleTube();
  this->ComputeTubePointScalesAndWeights();
  this->ComputeCenterRotation();

  this->m_Interpolator->SetInputImage( this->m_FixedImage );
  this->m_DerivativeImageFunction->SetInputImage( this->m_FixedImage );
}

/** Subsample the MovingSpatialObject tubenet  */
// WARNING: I think there is a bug from Stephen method with the
// index incrementation. I currently replicated the exact same behaviour,
// but the method must be validated.
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::SubSampleTube()
{
  SizeValueType tubeSize = 0;
  const OffsetValueType step = m_Sampling / 2 - 1;

  typename TubeNetType::Pointer newTubeNet = TubeNetType::New();
  typename TubeNetType::ChildrenListType* tubeList = this->GetTubes();
  typedef typename TubeNetType::ChildrenListType::iterator TubesIteratorType;
  for( TubesIteratorType tubeIterator = tubeList->begin();
    tubeIterator != tubeList->end();
    ++tubeIterator )
    {
    typename TubeType::Pointer newTube = TubeType::New();
    TubeType* currentTube =
      static_cast< TubeType * >( ( *tubeIterator ).GetPointer() );

    currentTube->RemoveDuplicatePoints();
    currentTube->ComputeTangentAndNormals();

    tubeSize = currentTube->GetPoints().size();
    if ( static_cast< OffsetValueType >( tubeSize ) > this->m_Sampling )
      {
      OffsetValueType loopIndex = 0;
      OffsetValueType skippedPoints = 0;
      const typename TubeType::PointListType & currentTubePoints
        = currentTube->GetPoints();
      typedef typename TubeType::PointListType::const_iterator
        TubePointIteratorType;
      for ( TubePointIteratorType tubePointIterator = currentTubePoints.begin();
            tubePointIterator != currentTubePoints.end();
            ++tubePointIterator, ++loopIndex )
        {
        if( m_Sampling > 1 )
          {
          if ( std::distance(tubePointIterator, currentTubePoints.end() - 1 )
               > step )
            {
            if( tubePointIterator != currentTubePoints.begin() )
              {
              tubePointIterator += step;
              }
            skippedPoints = ( loopIndex * ( step + 1 ) * 2 ) + 1;
            }
          else
            {
            skippedPoints = step + 2;
            tubePointIterator = currentTubePoints.end();
            }
          }

        // TODO Why +10 ?!
        if( tubePointIterator != currentTubePoints.end()
            && ( ( skippedPoints + 10 ) < static_cast< OffsetValueType >( tubeSize )) )
          {
          newTube->GetPoints().push_back( *( tubePointIterator ) );
          }

        if( m_Sampling > 1 )
          {
          if ( std::distance(tubePointIterator, currentTubePoints.end() -1 )
               > step )
            {
            tubePointIterator += step;
            }
          else
            {
            tubePointIterator = currentTubePoints.end() - 1;
            }
          }
        }
      }
    newTubeNet->AddSpatialObject( newTube );
    }
  itkDebugMacro( "Number of Points for the metric = " << m_NumberOfPoints );

  this->SetMovingSpatialObject( newTubeNet );
  delete tubeList;
}

template < class TFixedImage, class TMovingSpatialObject >
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ComputeCenterRotation()
{
  for ( unsigned int i = 0; i<ImageDimension; ++i )
    {
    this->m_CenterOfRotation[i] /= m_SumWeight;
    }

  itkDebugMacro( "Center of Rotation = "
              << this->m_CenterOfRotation[0] << "  " \
              << this->m_CenterOfRotation[1] << "  " \
              << this->m_CenterOfRotation[2] );
  itkDebugMacro( "Extent = " << m_Extent );
}

template < class TFixedImage, class TMovingSpatialObject >
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ComputeTubePointScalesAndWeights()
{
  typename TubeNetType::ChildrenListType* tubeList = this->GetTubes();
  typedef typename TubeNetType::ChildrenListType::iterator TubesIteratorType;

  for( TubesIteratorType tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       ++tubeIterator )
    {
    TubeType* currentTube =
      static_cast< TubeType * >( ( *tubeIterator ).GetPointer() );
    const typename TubeType::PointListType & currentTubePoints
      = currentTube->GetPoints();
    typedef typename TubeType::PointListType::const_iterator
      TubePointIteratorType;
    for ( TubePointIteratorType tubePointIterator = currentTubePoints.begin();
          tubePointIterator != currentTubePoints.end();
          ++tubePointIterator )
      {
      const InternalComputationValueType val =
        -2.0 * ( tubePointIterator->GetRadius() );
      const InternalComputationValueType weight =
        2.0 / ( 1.0 + exp( val ) );

      m_Weight.push_back( weight );

      for( unsigned int ii = 0; ii < ImageDimension; ++ii )
        {
        this->m_CenterOfRotation[ii] +=
          weight * ( tubePointIterator->GetPosition() )[ii];
        }
      m_SumWeight += weight;
      m_NumberOfPoints++;
      }
    }

  delete tubeList;
}

/** Get tubes contained within the Spatial Object */
// WARNING:
// Method might use GetMaximumDepth from ITK.
// Patch pushed in ITKv4, waiting for validation.
template < class TFixedImage, class TMovingSpatialObject>
typename ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>::TubeNetType::ChildrenListType*
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
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
// TODO Do not pass the parameter as arguments use instead
// the transform parameters previously set.
template < class TFixedImage, class TMovingSpatialObject>
typename ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>::MeasureType
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::GetValue( const ParametersType & parameters ) const
{
  itkDebugMacro( "**** Get Value ****" );
  itkDebugMacro( "Parameters = " << parameters );

  MeasureType matchMeasure = 0.0;
  InternalComputationValueType sumWeight = 0.0;
  InternalComputationValueType count = 0.0;
  InternalComputationValueType opR;
  InternalComputationValueType scale = this->m_InitialScale;

  std::list< InternalComputationValueType >::const_iterator weightIterator;
  weightIterator = m_Weight.begin();

  // TODO change the place where you set the parameters !
  //this->m_Transform->SetParameters( parameters );

  GroupSpatialObject<3>::ChildrenListType::iterator tubeIterator;
  typename TubeNetType::ChildrenListType* tubeList = GetTubes();
  for( tubeIterator = tubeList->begin();
       tubeIterator != tubeList->end();
       tubeIterator++ )
    {
    TubeType* currTube = static_cast<TubeType*>(
      ( *tubeIterator ).GetPointer() );

    typename std::vector<TubePointType>::iterator pointIterator;
    for( pointIterator = currTube->GetPoints().begin();
         pointIterator != currTube->GetPoints().end();
         ++pointIterator )
        {
        itk::Point<double, 3> inputPoint = pointIterator->GetPosition();
        /*static itk::Point<double, 3> point;
        Matrix<double, 3, 3> matrix =  GetTransform()->GetRotationMatrix();

        point =  matrix * inputPoint + GetTransform()->GetOffset();

        vnl_vector_fixed<double, 3> rotationOffset = matrix * ;

        point[0] += this->m_CenterOfRotation[0] - rotationOffset[0];
        point[1] += this->m_CenterOfRotation[1] - rotationOffset[1];
        point[2] += this->m_CenterOfRotation[2] - rotationOffset[2];

        // TODO
        // Need to use interpolator intead

        itk::Index<3> index;
        index[0] = static_cast<unsigned int>( point[0] );
        index[1] = static_cast<unsigned int>( point[1] );
        index[2] = static_cast<unsigned int>( point[2] );*/

        if( this->IsInside( inputPoint ) )
          {
          sumWeight += *weightIterator;
          count++;
          opR = pointIterator->GetRadius();
          opR = std::max( opR, 0.5 );

          scale = opR * m_Kappa;

          Vector<double, 3> v2;
          for( unsigned int i = 0; i < 3; ++i )
            {
            v2[i] = pointIterator->GetNormal1()[i];
            }

            matchMeasure += *weightIterator * fabs(
              ComputeLaplacianMagnitude( &v2, scale ) );
            }
        else
          {
          matchMeasure -= m_ImageMax;
          }

        weightIterator++;
        }
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
template < class TFixedImage, class TMovingSpatialObject>
double
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ComputeLaplacianMagnitude( Vector< InternalComputationValueType, 3> *v,
  const InternalComputationValueType & scale ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // m_CurrentPoint with a second derivative of a gaussian
  const InternalComputationValueType scaleSquared = scale * scale;
  const InternalComputationValueType scaleExtentProduct = scale * m_Extent;
  InternalComputationValueType result = 0.0;
  InternalComputationValueType wI = 0.0;
  unsigned int n = 0;

  itk::Index<3> index;
  for( double dist = -scaleExtentProduct; dist <= scaleExtentProduct; ++dist )
    {
    index[0] =
      static_cast<unsigned int>( m_CurrentPoint[0] + dist * v->GetElement(0) );
    index[1] =
      static_cast<unsigned int>( m_CurrentPoint[1] + dist * v->GetElement(1) );
    index[2] =
      static_cast<unsigned int>( m_CurrentPoint[2] + dist * v->GetElement(2) );

    double distSquared = std::pow( dist, 2 );
    itk::Size<3> size =
      this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if( index[0] >= 0 && ( index[0] < static_cast<unsigned int>( size[0] ) )
      && index[1] >= 0 && ( index[1] < static_cast<unsigned int>( size[1] ) )
      && index[2] >= 0 && ( index[2] < static_cast<unsigned int>( size[2] ) ) )
      {
      wI += ( -1 + ( distSquared / scaleSquared ) )
        * exp( -0.5 * distSquared / scaleSquared );
      n++;
      }
    }

  double error = wI / n;
  for( double dist = -scaleExtentProduct; dist <= scaleExtentProduct; ++dist )
    {
    double distSquared = std::pow( dist, 2 );
    wI = ( -1 + ( distSquared / scaleSquared ) )
      * exp( -0.5 * distSquared / scaleSquared ) - error;

    index[0] =
      static_cast<unsigned int>( m_CurrentPoint[0] + dist * v->GetElement(0) );
    index[1] =
      static_cast<unsigned int>( m_CurrentPoint[1] + dist * v->GetElement(1) );
    index[2] =
      static_cast<unsigned int>( m_CurrentPoint[2] + dist * v->GetElement(2) );

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

template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
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


/** Get the Derivative Measure */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::GetDerivative( const ParametersType & parameters,
                 DerivativeType & derivative ) const
{
  derivative = DerivativeType( SpaceDimension );

  this->m_Transform->SetParameters( parameters );

  InternalComputationValueType scale = this->m_InitialScale;
  InternalComputationValueType opR;

  itkDebugMacro( "**** Get Derivative ****" );
  itkDebugMacro( "parameters = "<< parameters )

  vnl_matrix<double> biasV( 3, 3, 0 );
  vnl_matrix<double> biasVI( 3, 3, 0 );


  std::list<double>::const_iterator weightIterator;
  weightIterator = m_Weight.begin();

  unsigned int count = 0;

  double dPosition[3] = { 0, 0, 0 };
  double dAngle[3] = { 0, 0, 0 };
  double dXProj1, dXProj2;

  vnl_vector<double> tV( 3 );
  vnl_matrix<double> tM( 3, 3 );
  vnl_vector<double> v1T( 3 );
  vnl_vector<double> v2T( 3 );
  double angleDelta[3];

  double sumWeight = 0;

  typedef itk::Vector<double, 3>    ITKVectorType;
  typedef std::list<ITKVectorType>  ListType;
  ListType dXTlist;
  itk::FixedArray<Point<double, 3>, 5000> XTlist;

  unsigned int listindex = 0;
  typename TubeNetType::ChildrenListType* tubeList = GetTubes();
  typename TubeNetType::ChildrenListType::iterator tubeIterator = tubeList->begin();
  for( ; tubeIterator != tubeList->end(); ++tubeIterator )
    {
    TubeType* currTube = static_cast<TubeType*>(
      ( *tubeIterator ).GetPointer() );

    typename std::vector<TubePointType>::iterator pointIterator;
    for( pointIterator = currTube->GetPoints().begin();
         pointIterator != currTube->GetPoints().end();
         ++pointIterator )
      {
      InputPointType inputPoint = pointIterator->GetPosition();
      itk::Point<double, 3> point;
      typename TransformType::MatrixType matrix =
        this->GetTransform()->GetMatrix();

      point =  matrix * inputPoint + GetTransform()->GetOffset();

      CenterOfRotationType rotationOffset = matrix * this->m_CenterOfRotation;

      point[0] += this->m_CenterOfRotation[0] - rotationOffset[0];
      point[1] += this->m_CenterOfRotation[1] - rotationOffset[1];
      point[2] += this->m_CenterOfRotation[2] - rotationOffset[2];

      itk::Index<3> index;
      index[0] = static_cast<unsigned int>( point[0] );
      index[1] = static_cast<unsigned int>( point[1] );
      index[2] = static_cast<unsigned int>( point[2] );

      if( this->IsInside( inputPoint ) )
        {
        XTlist[listindex++] = m_CurrentPoint;
        sumWeight += *weightIterator;
        count++;
        opR = pointIterator->GetRadius();
        opR = std::max( opR, 0.5 );

        scale = opR * m_Kappa;

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

        dXProj1 = ComputeThirdDerivatives( &v1, scale );
        dXProj2 = ComputeThirdDerivatives( &v2, scale );

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

  if( sumWeight == 0 )
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

  weightIterator = m_Weight.begin();
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

/** Get both the match Measure and theDerivative Measure */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::GetValueAndDerivative( const ParametersType & parameters,
                         MeasureType&  Value,
                         DerivativeType  & Derivative ) const
{
  Value = GetValue( parameters );
  GetDerivative( parameters, Derivative );
}

template < class TFixedImage, class TMovingSpatialObject>
double
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ComputeThirdDerivatives( Vector< InternalComputationValueType, 3> *v,
  const InternalComputationValueType & scale ) const
{
  // We convolve the 1D signal defined by the direction v at point
  // m_CurrentPoint with a second derivative of a gaussian
  InternalComputationValueType result = 0.0;
  InternalComputationValueType wI = 0.0;
  itk::Index<3> index;

  InternalComputationValueType wTotalX = 0.0;
  const InternalComputationValueType scaleSquared = scale * scale;
  const InternalComputationValueType scaleExtentProduct = scale * m_Extent;

  for( double dist = -scaleExtentProduct;
       dist <= scaleExtentProduct;
       dist += 0.1 )
    {
    wI = 2 * dist * exp( -0.5 * std::pow( dist, 2 ) / scaleSquared );

    wTotalX += fabs( wI );

    index[0] =
      static_cast<unsigned int>( m_CurrentPoint[0] + dist * v->GetElement(0) );
    index[1] =
      static_cast<unsigned int>( m_CurrentPoint[1] + dist * v->GetElement(1) );
    index[2] =
      static_cast<unsigned int>( m_CurrentPoint[2] + dist * v->GetElement(2) );

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

/** Test whether the specified point is inside
 * Thsi method overload the one in the ImageMapper class
 * \warning This method cannot be safely used in more than one thread at
 * a time.
 * \sa Evaluate(); */
template < class TFixedImage, class TMovingSpatialObject>
bool
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::IsInside( const InputPointType & point ) const
{
  typename TransformType::MatrixType matrix =  this->GetTransform()->GetMatrix();
  m_CurrentPoint = matrix * point + this->GetTransform()->GetOffset();

  CenterOfRotationType centerOfRotation;
  for( unsigned int ii = 0; ii < TubeDimension; ++ii )
    {
    centerOfRotation[ii] = this->m_CenterOfRotation[ii];
    }

  CenterOfRotationType rotationOffset = matrix * centerOfRotation;

  m_CurrentPoint[0] += centerOfRotation[0] - rotationOffset[0];
  m_CurrentPoint[1] += centerOfRotation[1] - rotationOffset[1];
  m_CurrentPoint[2] += centerOfRotation[2] - rotationOffset[2];

  return ( this->m_Interpolator->IsInsideBuffer( m_CurrentPoint ) );
}

/** Clamp the point bounds to the fixed Image */
template < class TFixedImage, class TMovingSpatialObject>
void
ImageToTubeRigidMetric<TFixedImage, TMovingSpatialObject>
::ClampPointBoundsToImage(int bounds[6])
{
  SizeType size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();

  bounds[0] = vnl_math_max( bounds[0], 0 );
  bounds[1] = vnl_math_min( bounds[1], static_cast<int>( size[0] ) - 1 );
  bounds[2] = vnl_math_max( bounds[2], 0 );
  bounds[3] = vnl_math_min( bounds[3], static_cast<int>( size[1] ) - 1 );
  bounds[4] = vnl_math_max( bounds[4], 0 );
  bounds[5] = vnl_math_min( bounds[5], static_cast<int>( size[2] ) - 1 );
}

} // end namespace itk

#endif
