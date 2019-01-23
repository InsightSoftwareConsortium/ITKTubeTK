/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itktubeRadiusExtractor2_hxx
#define __itktubeRadiusExtractor2_hxx

#include "itktubeRadiusExtractor2.h"

#include "tubeMatrixMath.h"
#include "tubeTubeMath.h"
#include "tubeUserFunction.h"
#include "tubeGoldenMeanOptimizer1D.h"
#include "tubeSplineApproximation1D.h"

#include <itkMinimumMaximumImageFilter.h>

#include <vnl/vnl_math.h>

namespace itk
{

namespace tube
{

template< class TInputImage >
class LocalMedialnessSplineValueFunction
: public ::tube::UserFunction< int, double >
{
public:

  LocalMedialnessSplineValueFunction( RadiusExtractor2< TInputImage > *
    _radiusExtractor )
    {
    m_Value = 0;

    m_RadiusExtractor = _radiusExtractor;
    }

  const double & Value( const int & x )
    {
    m_Value = m_RadiusExtractor->GetKernelMedialness( x *
      m_RadiusExtractor->GetRadiusTolerance() );

    return m_Value;
    }

private:

  double   m_Value;

  typename RadiusExtractor2< TInputImage >::Pointer   m_RadiusExtractor;
};


/** Constructor */
template< class TInputImage >
RadiusExtractor2<TInputImage>
::RadiusExtractor2( void )
{
  m_Image = NULL;

  m_DataMin = 0;
  m_DataMax = -1;

  m_RadiusStart = 1.5;
  m_RadiusMin = 0.33;
  m_RadiusMax = 15.0;
  m_RadiusStep = 0.25;
  m_RadiusTolerance = 0.125;

  m_RadiusCorrectionFunction = RADIUS_CORRECTION_NONE;

  m_MinMedialness = 0.15;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.1;

  m_NumKernelPoints = 7;
  m_KernelTubePoints.resize( m_NumKernelPoints );

  m_KernelPointStep = 14;
  m_KernelStep = 30;
  m_KernelExtent = 1.6;

  m_KernelValues.clear();
  m_KernelDistances.clear();
  m_KernelTangentDistances.clear();

  m_KernelOptimalRadius = 0;
  m_KernelOptimalRadiusMedialness = 0;
  m_KernelOptimalRadiusBranchness = 0;

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
}

/** Destructor */
template< class TInputImage >
RadiusExtractor2<TInputImage>
::~RadiusExtractor2( void )
{
}

/** Set the input image */
template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetInputImage( typename ImageType::Pointer inputImage )
{
  m_Image = inputImage;

  if( m_Image )
    {
    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter =
      MinMaxFilterType::New();
    minMaxFilter->SetInput( m_Image );
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();
    if( this->GetDebug() )
      {
      std::cout << "RadiusExtractor2: SetInputImage: Minimum = "
        << m_DataMin << std::endl;
      std::cout << "RadiusExtractor2: SetInputImage: Maximum = "
        << m_DataMax << std::endl;
      }
    }
}

/** Compute the medialness at a kernel */
template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::GetPointVectorMeasures( std::vector< TubePointType > & points,
  double pntR,
  double & mness,
  double & bness,
  bool doBNess )
{
  ::tube::ComputeVectorTangentsAndNormals( points );

  if( this->GetDebug() )
    {
    std::cout << "Compute values at point" << std::endl;
    }

  unsigned int tempNumPoints = this->GetNumKernelPoints();
  int numPoints = points.size();
  this->SetNumKernelPoints( numPoints );

  this->SetKernelTubePoints( points );

  this->GenerateKernel();

  mness = this->GetKernelMedialness( pntR );
  if( doBNess )
    {
    bness = this->GetKernelBranchness( pntR );
    }

  this->SetNumKernelPoints( tempNumPoints );
}

/** Compute the Optimal scale */
template< class TInputImage >
bool
RadiusExtractor2<TInputImage>
::GetPointVectorOptimalRadius( std::vector< TubePointType > & points,
  double & r0,
  double rMin,
  double rMax,
  double rStep,
  double rTolerance )
{
  unsigned int tempNumPoints = this->GetNumKernelPoints();
  unsigned int numPoints = points.size();
  this->SetNumKernelPoints( numPoints );

  this->SetKernelTubePoints( points );

  double tempXStart = this->GetRadiusStart();
  this->SetRadiusStart( r0 );
  double tempXMin = this->GetRadiusMin();
  this->SetRadiusMin( rMin );
  double tempXMax = this->GetRadiusMax();
  this->SetRadiusMax( rMax );
  double tempXStep = this->GetRadiusStep();
  this->SetRadiusStep( rStep );
  double tempXTolerance = this->GetRadiusTolerance();
  this->SetRadiusTolerance( rTolerance );

  this->GenerateKernel();

  this->UpdateKernelOptimalRadius();

  this->SetRadiusStart( tempXStart );
  this->SetRadiusMin( tempXMin );
  this->SetRadiusMax( tempXMax );
  this->SetRadiusStep( tempXStep );
  this->SetRadiusTolerance( tempXTolerance );

  if( tempNumPoints != numPoints )
    {
    this->SetNumKernelPoints( tempNumPoints );
    }

  r0 = this->GetKernelOptimalRadius();

  return true;
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetNumKernelPoints( unsigned int _numPoints )
{
  m_NumKernelPoints = _numPoints;

  m_KernelTubePoints.resize( m_NumKernelPoints );
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::GenerateKernel( void )
{
  ITKIndexType minX;
  ITKIndexType maxX;

  unsigned int dimension = m_Image->GetImageDimension();

  ITKIndexType buffer;
  for( unsigned int i = 0; i < dimension; ++i )
    {
    buffer[i] = static_cast< unsigned int >(
      this->GetKernelExtent() * this->GetRadiusMax() );
    }

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = m_KernelTubePoints.begin();
  for( unsigned int i = 0; i < dimension; ++i )
    {
    minX[i] = static_cast< int >( m_KernelTubePoints[0].GetPosition()[i]
      - buffer[i] );
    maxX[i] = static_cast< int >( m_KernelTubePoints[0].GetPosition()[i]
      + buffer[i] );
    }
  ++pntIter;
  int tempI;
  while( pntIter != m_KernelTubePoints.end() )
    {
    for( unsigned int i = 0; i < dimension; ++i )
      {
      tempI = static_cast< int >( pntIter->GetPosition()[i] - buffer[i] );
      if( tempI < minX[i] )
        {
        minX[i] = tempI;
        }
      tempI = static_cast< int >( pntIter->GetPosition()[i] + buffer[i] );
      if( tempI > maxX[i] )
        {
        maxX[i] = tempI;
        }
      }
    ++pntIter;
    }
  unsigned int kernelSize = maxX[0] - minX[0] + 1;
  for( unsigned int i = 1; i < dimension; ++i )
    {
    kernelSize *= ( maxX[i] - minX[i] + 1 );
    }
  m_KernelValues.resize( kernelSize );
  m_KernelDistances.resize( kernelSize );
  m_KernelTangentDistances.resize( kernelSize );

  unsigned int count = 0;
  ITKIndexType x = minX;
  bool done = false;
  while( !done )
    {
    if( m_Image->GetLargestPossibleRegion().IsInside( x ) )
      {
      m_KernelValues[ count ] = ( m_Image->GetPixel( x ) - m_DataMin )
        / ( m_DataMax - m_DataMin );
      if( m_KernelValues[ count ] < 0 )
        {
        m_KernelValues[ count ] = 0;
        }
      else if( m_KernelValues[ count ] > 1 )
        {
        m_KernelValues[ count ] = 1;
        }
      ITKPointType p;
      for( unsigned int i = 0; i < dimension; ++i )
        {
        p[ i ] = x[ i ];
        }

      unsigned int pntCount = 0;
      pntIter = m_KernelTubePoints.begin();

      double pntTangentDist = 0;
      double minTangentDist = 0;
      double minNormalDist = 0;
      int minTangentDistCount = -1;
      while( pntIter != m_KernelTubePoints.end() )
        {
        ITKVectorType pDiff = pntIter->GetPosition() - p;
        double d1 = 0;
        for( unsigned int i = 0; i < dimension; ++i )
          {
          d1 += pDiff[i] * pntIter->GetTangent()[i];
          }
        pntTangentDist = d1 * d1;
        if( pntTangentDist < minTangentDist || minTangentDistCount == -1 )
          {
          minTangentDist = pntTangentDist;
          minTangentDistCount = pntCount;
          d1 = 0;
          for( unsigned int i = 0; i < dimension; ++i )
            {
            d1 += pDiff[i] * pntIter->GetNormal1()[i];
            }
          minNormalDist = d1 * d1;
          if( dimension == 3 )
            {
            double d2 = 0;
            for( unsigned int i = 0; i < dimension; ++i )
              {
              d2 += pDiff[i] * pntIter->GetNormal2()[i];
              }
            minNormalDist += d2 * d2;
            }
          }
        ++pntIter;
        ++pntCount;
        }
      m_KernelDistances[ count ] = std::sqrt( minNormalDist );
      m_KernelTangentDistances[ count ] = std::sqrt( minTangentDist );
      }
    else
      {
      m_KernelValues[ count ] = 0;
      m_KernelDistances[ count ] = this->GetKernelExtent() *
        this->GetRadiusMax() + 1;
      m_KernelTangentDistances[ count ] = this->GetKernelExtent() *
        this->GetRadiusMax() + 1;
      }
    ++count;
    unsigned int d = 0;
    while( d < dimension && ++x[d] > maxX[d] )
      {
      x[d] = minX[d];
      ++d;
      }
    if( d >= dimension )
      {
      done = true;
      }
    }
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetKernelTubePoints( const std::vector< TubePointType > & tubePoints )
{
  if( tubePoints.size() != m_NumKernelPoints )
    {
    std::cerr << "Error: number of kernel points not equal to expected."
      << std::endl;
    std::cerr << "   TubePointsSize = " << tubePoints.size() << std::endl;
    std::cerr << "   NumKernelPoints = " << m_NumKernelPoints << std::endl;
    }

  for( unsigned int i=0; i<m_NumKernelPoints; ++i )
    {
    m_KernelTubePoints[ i ] = tubePoints[ i ];
    }

  ::tube::ComputeVectorTangentsAndNormals( m_KernelTubePoints );
}

template< class TInputImage >
double
RadiusExtractor2<TInputImage>
::GetKernelMedialness( double r )
{
  if( r < m_RadiusMin )
    {
    double factor = ( m_RadiusMin - r ) / m_RadiusTolerance;
    double m0 = this->GetKernelMedialness( m_RadiusMin );
    double m1 = this->GetKernelMedialness( m_RadiusMin + m_RadiusTolerance );
    return m0 - factor * std::fabs( m0 - m1 );
    }
  else if( r > m_RadiusMax )
    {
    double factor = ( r - m_RadiusMax ) / m_RadiusTolerance;
    double m0 = this->GetKernelMedialness( m_RadiusMax );
    double m1 = this->GetKernelMedialness( m_RadiusMax - m_RadiusTolerance );
    return m0 - factor * std::fabs( m0 - m1 );
    }

  double pVal = 0;
  double nVal = 0;

  std::vector< double >::iterator iterDist;
  std::vector< double >::iterator iterTanDist;
  std::vector< double >::iterator iterValue;

  double areaR = r * r * vnl_math::pi;
  double distMax = ( r + 1.0 );
  double areaMax = distMax * distMax * vnl_math::pi;
  double areaNeg = areaMax - areaR;
  double distMin2 = ( areaR - areaNeg ) / vnl_math::pi;
  double distMin = 0;
  if( distMin2 > 0 )
    {
    distMin = std::sqrt( distMin2 );
    }
  double areaMin = distMin * distMin * vnl_math::pi;
  double areaPos = areaR - areaMin;
  if( this->GetDebug() )
    {
    std::cout << "R = " << r << std::endl;
    std::cout << "   Dist = " << distMin << " - " << distMax << std::endl;
    std::cout << "   Area = " << areaPos << " - " << areaNeg << std::endl;
    }

  const int histoBins = 500;
  unsigned int histoPos[histoBins];
  unsigned int histoNeg[histoBins];
  unsigned int histoPosCount = 0;
  unsigned int histoNegCount = 0;
  for( int i=0; i<histoBins; ++i )
    {
    histoPos[i] = 0;
    histoNeg[i] = 0;
    }
  int bin = 0;
  iterDist = m_KernelDistances.begin();
  iterTanDist = m_KernelTangentDistances.begin();
  iterValue = m_KernelValues.begin();
  while( iterDist != m_KernelDistances.end() )
    {
    if( ( ( *iterTanDist ) < r || ( *iterTanDist ) < 1 )
      && ( *iterDist ) >= distMin && ( *iterDist ) <= distMax )
      {
      bin = ( *iterValue ) * histoBins;
      if( bin < 0 )
        {
        bin = 0;
        }
      else if( bin > histoBins - 1 )
        {
        bin = histoBins - 1;
        }
      if( ( *iterDist ) <= r )
        {
        ++histoPos[bin];
        ++histoPosCount;
        }
      else
        {
        ++histoNeg[bin];
        ++histoNegCount;
        }
      }
    ++iterValue;
    ++iterDist;
    ++iterTanDist;
    }

  int binCount = 0;
  pVal = 0;
  if( histoPosCount > 2 )
    {
    bin = histoBins - 1;
    while( binCount < 0.5 * histoPosCount && bin > 0 )
      {
      binCount += histoPos[bin];
      --bin;
      }
    pVal = ( bin + 0.5 ) / histoBins;
    }
  nVal = 1;
  if( histoNegCount > 2 )
    {
    binCount = 0;
    bin = 0;
    while( binCount < 0.5 * histoNegCount && bin < histoBins )
      {
      binCount += histoNeg[bin];
      ++bin;
      }
    nVal = ( bin - 0.5 ) / histoBins;
    }

  if( this->GetDebug() )
    {
    std::cout << "   Count = " << histoPosCount << " - " << histoNegCount
      << std::endl;
    std::cout << "   Val = " << pVal << " - " << nVal << " = " << pVal-nVal
      << std::endl;
    }
  double medialness = pVal - nVal;

  return medialness;
}

template< class TInputImage >
double
RadiusExtractor2<TInputImage>
::GetKernelBranchness( double r )
{
  if( r < m_RadiusMin )
    {
    double factor = ( m_RadiusMin - r ) / m_RadiusTolerance;
    double m0 = this->GetKernelBranchness( m_RadiusMin );
    double m1 = this->GetKernelBranchness( m_RadiusMin
      + m_RadiusTolerance );
    return m0 - factor * std::fabs( m0 - m1 );
    }
  else if( r > m_RadiusMax )
    {
    double factor = ( r - m_RadiusMax ) / m_RadiusTolerance;
    double m0 = this->GetKernelBranchness( m_RadiusMax );
    double m1 = this->GetKernelBranchness( m_RadiusMax
      - m_RadiusTolerance );
    return m0 - factor * std::fabs( m0 - m1 );
    }
  double pVal = 0;
  double nVal = 0;

  std::vector< double >::iterator iterDist;
  std::vector< double >::iterator iterTanDist;
  std::vector< double >::iterator iterValue;

  double distMax = r * this->GetKernelExtent();
  double distMin = 0;

  const int histoBins = 500;
  unsigned int histoPos[histoBins];
  unsigned int histoNeg[histoBins];
  unsigned int histoPosCount = 0;
  unsigned int histoNegCount = 0;
  for( unsigned int i=0; i<histoBins; ++i )
    {
    histoPos[i] = 0;
    histoNeg[i] = 0;
    }
  int bin = 0;
  iterDist = m_KernelDistances.begin();
  iterTanDist = m_KernelTangentDistances.begin();
  iterValue = m_KernelValues.begin();
  while( iterDist != m_KernelDistances.end() )
    {
    if( ( *iterTanDist ) < r
      && ( *iterDist ) >= distMin && ( *iterDist ) <= distMax )
      {
      bin = ( *iterValue ) * histoBins;
      if( bin < 0 )
        {
        bin = 0;
        }
      else if( bin > histoBins - 1 )
        {
        bin = histoBins - 1;
        }
      if( ( *iterDist ) <= r )
        {
        ++histoPos[bin];
        ++histoPosCount;
        }
      else
        {
        ++histoNeg[bin];
        ++histoNegCount;
        }
      }
    ++iterValue;
    ++iterDist;
    ++iterTanDist;
    }

  int binCount = 0;
  pVal = 0;
  if( histoPosCount > 1 )
    {
    bin = 0;
    while( binCount < 0.5 * histoPosCount && bin < histoBins )
      {
      binCount += histoPos[bin];
      ++bin;
      }
    pVal = ( bin + 0.5 ) / histoBins;
    }
  nVal = 1;
  if( histoNegCount > 1 )
    {
    binCount = 0;
    bin = histoBins - 1;
    while( binCount < 0.75 * histoNegCount && bin > 1 )
      {
      binCount += histoNeg[bin];
      --bin;
      }
    nVal = ( bin - 0.5 ) / histoBins;
    }

  double branchness = 1.0 - std::fabs( pVal - nVal );

  return branchness;
}

template< class TInputImage >
bool
RadiusExtractor2<TInputImage>
::UpdateKernelOptimalRadius( void )
{
  ::tube::UserFunction< int, double > * myFunc = new
    LocalMedialnessSplineValueFunction< TInputImage >( this );

  ::tube::GoldenMeanOptimizer1D * opt = new
    ::tube::GoldenMeanOptimizer1D();
  opt->SetXStep( m_RadiusStep / m_RadiusTolerance );
  opt->SetTolerance( 1 );
  opt->SetMaxIterations( 20 );
  opt->SetSearchForMin( false );

  int xMin = 0;
  int xMax = 1;
  double x = 0.5;
  switch( m_RadiusCorrectionFunction )
    {
    default:
    case RADIUS_CORRECTION_NONE:
      {
      xMin = ( int )std::ceil( m_RadiusMin / m_RadiusTolerance );
      xMax = ( int )std::floor( m_RadiusMax / m_RadiusTolerance );
      x = m_RadiusStart / m_RadiusTolerance;
      break;
      }
    case RADIUS_CORRECTION_FOR_BINARY_IMAGE:
      {
      xMin = ( int )std::ceil( m_RadiusMin / m_RadiusTolerance );
      xMax = ( int )std::floor( m_RadiusMax / m_RadiusTolerance );
      x = m_RadiusStart / m_RadiusTolerance;
      break;
      }
    case RADIUS_CORRECTION_FOR_CTA:
      {
      xMin = ( int )std::ceil( m_RadiusMin / m_RadiusTolerance );
      xMax = ( int )std::floor( m_RadiusMax / m_RadiusTolerance );
      x = m_RadiusStart / m_RadiusTolerance;
      break;
      }
    case RADIUS_CORRECTION_FOR_MRA:
      {
      xMin = ( int )std::ceil( m_RadiusMin / m_RadiusTolerance );
      xMax = ( int )std::floor( m_RadiusMax / m_RadiusTolerance );
      x = m_RadiusStart / m_RadiusTolerance;
      break;
      }
    };

  opt->SetXMin( xMin );
  opt->SetXMax( xMax );

  ::tube::SplineApproximation1D * spline = new
  ::tube::SplineApproximation1D( myFunc, opt );

  //spline->SetClip( true );
  spline->SetXMin( xMin );
  spline->SetXMax( xMax );

  double xVal = myFunc->Value( x );
  bool result = spline->Extreme( &x, &xVal );

  switch( m_RadiusCorrectionFunction )
    {
    default:
    case RADIUS_CORRECTION_NONE:
      {
      m_KernelOptimalRadius = x * m_RadiusTolerance;
      break;
      }
    case RADIUS_CORRECTION_FOR_BINARY_IMAGE:
      {
      m_KernelOptimalRadius = ( x * m_RadiusTolerance
        * x * m_RadiusTolerance ) / 24 + 0.5;
      break;
      }
    case RADIUS_CORRECTION_FOR_CTA:
      {
      m_KernelOptimalRadius = x * m_RadiusTolerance;
      break;
      }
    case RADIUS_CORRECTION_FOR_MRA:
      {
      m_KernelOptimalRadius = x * m_RadiusTolerance;
      break;
      }
    };

  m_KernelOptimalRadiusMedialness = xVal;

  delete spline;
  delete opt;
  delete myFunc;

  return result;
}

template< class TInputImage >
bool
RadiusExtractor2<TInputImage>
::ExtractRadii( TubeType * tube )
{
  if( tube->GetPoints().size() == 0 )
    {
    return false;
    }

  tube->RemoveDuplicatePoints();
  ::tube::ComputeVectorTangentsAndNormals< TubePointType >(
    tube->GetPoints() );

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = tube->GetPoints().begin();
  while( pntIter != tube->GetPoints().end() )
    {
    pntIter->SetRadius( 0 );
    ++pntIter;
    }

  int pntCount = 0;
  pntIter = tube->GetPoints().begin();
  while( pntIter != tube->GetPoints().end() && pntIter->GetID() != 0 )
    {
    ++pntIter;
    ++pntCount;
    }

  if( pntIter == tube->GetPoints().end() )
    {
    if( this->GetDebug() )
      {
      std::cout << "Warning: PointID 0 not found. Using mid-point of tube."
        << std::endl;
      }
    pntIter = tube->GetPoints().begin();
    unsigned int psize = tube->GetPoints().size();
    for( unsigned int i=0; i<psize/2; i++ )
      {
      ++pntIter;
      }
    }
  else if( this->GetDebug() )
    {
    std::cout << "Found point " << ( *pntIter ).GetID() << std::endl;
    }

  double rStart0 = this->GetRadiusStart();
  double rStart = rStart0;
  for( int p = static_cast< int >( pntCount );
    p < static_cast< int >( tube->GetPoints().size() );
    p += this->GetKernelStep() )
    {
    this->SetRadiusStart( rStart );
    this->GenerateKernelTubePoints( p, tube );
    this->GenerateKernel();
    this->UpdateKernelOptimalRadius();
    this->RecordOptimaAtTubePoints( p, tube );
    rStart = this->GetKernelOptimalRadius();
    }

  rStart = rStart0;
  for( int p = static_cast< int >( pntCount )
    - this->GetKernelPointStep(); p >= 0;
    p -= this->GetKernelStep() )
    {
    this->SetRadiusStart( rStart );
    this->GenerateKernelTubePoints( p, tube );
    this->GenerateKernel();
    this->UpdateKernelOptimalRadius();
    this->RecordOptimaAtTubePoints( p, tube );
    rStart = this->GetKernelOptimalRadius();
    }

  if( this->GetDebug() )
    {
    std::cout << "Radius results:" << std::endl;
    pntIter = tube->GetPoints().begin();
    while( pntIter != tube->GetPoints().end() )
      {
      std::cout << "   " << pntIter->GetID() << " : "
        << pntIter->GetRadius() << std::endl;
      ++pntIter;
      }
    }

  return true;
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::GenerateKernelTubePoints( unsigned int tubePointNum,
  TubeType * tube )
{
  unsigned int tubeSize = tube->GetPoints().size();
  if( tubeSize < m_NumKernelPoints * m_KernelPointStep )
    {
    std::cerr << "RadiusExtractor: Tube length is too short" << std::endl;
    return;
    }

  int startP = tubePointNum - ( m_NumKernelPoints / 2 ) * m_KernelPointStep;

  int endP = startP + ( m_NumKernelPoints - 1 ) * m_KernelPointStep;

  unsigned int count = 0;
  for( int p = startP; p <= endP; p += m_KernelPointStep )
    {
    if( p < 0 )
      {
      typename TubeType::PointType p1 =
        tube->GetPoints()[ 0 ].GetPosition();
      typename TubeType::PointType p2 =
        tube->GetPoints()[ m_NumKernelPoints / 2 *
        m_KernelPointStep ].GetPosition();
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( -p / m_KernelPointStep ) );
        }
      m_KernelTubePoints[ count ].SetPosition( p2 );
      }
    else if( p > static_cast< int >( tubeSize ) - 1 )
      {
      typename TubeType::PointType p1 =
        tube->GetPoints()[ tubeSize - 1 ].GetPosition();
      typename TubeType::PointType p2 =
        tube->GetPoints()[ tubeSize - 1 - m_NumKernelPoints / 2 *
        m_KernelPointStep ].GetPosition();
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( ( p - ( tubeSize-1 ) ) /
            m_KernelPointStep ) );
        }
      m_KernelTubePoints[ count ].SetPosition( p2 );
      }
    else
      {
      m_KernelTubePoints[ count ] = tube->GetPoints()[ p ];
      }
    ++count;
    }

  ::tube::ComputeVectorTangentsAndNormals( m_KernelTubePoints );
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::RecordOptimaAtTubePoints( unsigned int tubePointNum,
  TubeType * tube )
{
  int tubeSize = tube->GetPoints().size();

  double r1 = this->GetKernelOptimalRadius();
  double m1 = this->GetKernelOptimalRadiusMedialness();
  double b1 = this->GetKernelOptimalRadiusBranchness();

  int startP = tubePointNum - ( m_NumKernelPoints / 2 ) * m_KernelPointStep;
  if( startP < 0 )
    {
    startP = 0;
    }
  double r0 = tube->GetPoints()[ startP ].GetRadius();
  double m0 = tube->GetPoints()[ startP ].GetMedialness();
  double b0 = tube->GetPoints()[ startP ].GetBranchness();
  if( r0 == 0 )
    {
    r0 = r1;
    m0 = m1;
    b0 = b1;
    }

  int endP = startP + ( m_NumKernelPoints - 1 ) * m_KernelPointStep;
  if( endP > tubeSize-1 )
    {
    endP = tubeSize-1;
    }
  double r2 = tube->GetPoints()[ endP ].GetRadius();
  double m2 = tube->GetPoints()[ endP ].GetMedialness();
  double b2 = tube->GetPoints()[ endP ].GetBranchness();
  if( r2 == 0 )
    {
    r2 = r1;
    m2 = m1;
    b2 = b1;
    }

  for( int p = startP; p <= endP; ++p )
    {
    if( p < static_cast< int >( tubePointNum ) )
      {
      double d = 1;
      if( static_cast< int >( tubePointNum ) != startP )
        {
        d = static_cast< double >( tubePointNum - p )
          / ( tubePointNum - startP );
        }
      tube->GetPoints()[ p ].SetRadius( d * r0 + ( 1 - d ) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m0 + ( 1 - d ) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b0 + ( 1 - d ) * b1 );
      }
    else
      {
      double d = 1;
      if( static_cast< int >( tubePointNum ) != endP )
        {
        d = static_cast< double >( p - tubePointNum )
          / ( endP - tubePointNum );
        }
      tube->GetPoints()[ p ].SetRadius( d * r2 + ( 1 - d ) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m2 + ( 1 - d ) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b2 + ( 1 - d ) * b1 );
      }
    }
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_Image.IsNotNull() )
    {
    os << indent << "Image = " << m_Image << std::endl;
    }
  else
    {
    os << indent << "Image = NULL" << std::endl;
    }
  os << indent << "DataMin = " << m_DataMin << std::endl;
  os << indent << "DataMax = " << m_DataMax << std::endl;

  os << indent << "RadiusStart = " << m_RadiusStart << std::endl;
  os << indent << "RadiusMin = " << m_RadiusMin << std::endl;
  os << indent << "RadiusMax = " << m_RadiusMax << std::endl;
  os << indent << "RadiusStep = " << m_RadiusStep << std::endl;
  os << indent << "RadiusTolerance = " << m_RadiusTolerance << std::endl;

  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;

  os << indent << "NumKernelPoints = " << m_NumKernelPoints << std::endl;
  os << indent << "KernelTubePoints = " << m_KernelTubePoints.size()
    << std::endl;
  os << indent << "KernelPointStep = " << m_KernelPointStep << std::endl;
  os << indent << "KernelStep = " << m_KernelStep << std::endl;
  os << indent << "KernelExtent = " << m_KernelExtent << std::endl;
  os << indent << "KernelValues = " << m_KernelValues.size() << std::endl;
  os << indent << "KernelDistances = " << m_KernelDistances.size()
    << std::endl;
  os << indent << "KernelTangentDistances = " << m_KernelTangentDistances.size()
    << std::endl;

  os << indent << "KernelOptimalRadius = " << m_KernelOptimalRadius
    << std::endl;
  os << indent << "KernelOptimalRadiusMedialness = "
    << m_KernelOptimalRadiusMedialness << std::endl;
  os << indent << "KernelOptimalRadiusBranchness = "
    << m_KernelOptimalRadiusBranchness
    << std::endl;

  os << indent << "IdleCallBack = " << m_IdleCallBack << std::endl;
  os << indent << "StatusCallBack = " << m_StatusCallBack << std::endl;
}

/**
 * Idle callback */
template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetIdleCallBack( bool ( *idleCallBack )() )
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Status Call back */
template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetStatusCallBack( void ( *statusCallBack )( const char *, const char *,
    int ) )
{
  m_StatusCallBack = statusCallBack;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeRadiusExtractor2_hxx )
