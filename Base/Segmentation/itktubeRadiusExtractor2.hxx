/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itktubeRadiusExtractor2_hxx
#define __itktubeRadiusExtractor2_hxx

#include "itktubeRadiusExtractor2.h"

#include "tubeMatrixMath.h"
#include "tubeTubeMath.h"
#include "tubeUserFunction.h"

#include <itkMinimumMaximumImageFilter.h>

namespace itk
{

namespace tube
{

/** Define the Medialness Function
 * \class RadiusExtractor2MedialnessFunc
 */
template< class TImage >
class RadiusExtractor2MedialnessFunc : public ::tube::UserFunction< int, double >
{
public:

  typedef itk::VesselTubeSpatialObject< TImage::ImageDimension >
                                                           TubeType;
  typedef typename TubeType::TubePointType                 TubePointType;

  RadiusExtractor2MedialnessFunc( RadiusExtractor2< TImage > *
    newRadiusExtractor )
    {
    m_RadiusExtractor = newRadiusExtractor;
    }

  const double & Value( const int & x )
    {
    double r = m_RadiusExtractor->GetRadiusStep() * x;
    m_Value = m_RadiusExtractor->GetKernelMedialness( r );
    return m_Value;
    }

  RadiusExtractor2< TImage > * GetRadiusExtractor( void  )
    {
    return m_RadiusExtractor;
    }

private:

  double                                    m_Value;
  RadiusExtractor2< TImage >              * m_RadiusExtractor;

}; // End class RadiusExtractorMedialnessFunc

/** Constructor */
template< class TInputImage >
RadiusExtractor2<TInputImage>
::RadiusExtractor2( void )
{
  m_Image = NULL;

  m_DataMin = 0;
  m_DataMax = -1;

  m_RadiusStart = 0.5;
  m_RadiusMin = 0.3;
  m_RadiusMax = 6.0;
  m_RadiusStep = 0.1;
  m_RadiusTolerance = 0.01;

  m_MinMedialness = 0.15;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.1;

  m_MedialnessFunc = new RadiusExtractor2MedialnessFunc<TInputImage>(
    this );

  m_MedialnessOpt.SetSearchForMin( false );

  m_MedialnessOptSpline = new SplineType( m_MedialnessFunc,
    &m_MedialnessOpt );

  m_NumKernelPoints = 7;
  m_KernelTubePoints.resize( m_NumKernelPoints );

  m_KernelPointStep = 7;
  m_KernelStep = 20;
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
  if( m_MedialnessFunc != NULL )
    {
    delete m_MedialnessFunc;
    }
  m_MedialnessFunc = NULL;

  if( m_MedialnessOptSpline != NULL )
    {
    delete m_MedialnessOptSpline;
    }
  m_MedialnessOptSpline = NULL;
}

/** Get the medialness operator */
template< class TInputImage >
typename RadiusExtractor2<TInputImage>::OptimizerType &
RadiusExtractor2<TInputImage>
::GetMedialnessOptimizer( void )
{
  return & m_MedialnessOpt;
}

/** Get the medialness operator */
template< class TInputImage >
typename RadiusExtractor2<TInputImage>::SplineType &
RadiusExtractor2<TInputImage>
::GetMedialnessOptimizerSpline( void )
{
  return m_MedialnessOptSpline;
}

/** Set the input image */
template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::SetInputImage( typename ImageType::Pointer inputImage  )
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
  if( pntR < m_MedialnessOptSpline->GetXMin() )
    {
    pntR = m_MedialnessOptSpline->GetXMin();
    }

  if( pntR > m_MedialnessOptSpline->GetXMax() )
    {
    pntR = m_MedialnessOptSpline->GetXMax();
    }

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
      m_KernelDistances[ count ] = vcl_sqrt( minNormalDist );
      m_KernelTangentDistances[ count ] = vcl_sqrt( minTangentDist );
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
    double factor = ( m_RadiusMin - r ) / m_RadiusStep;
    double m0 = this->GetKernelMedialness( m_RadiusMin );
    double m1 = this->GetKernelMedialness( m_RadiusMin + m_RadiusStep );
    return m0 - factor * vnl_math_abs(m0 - m1);
    }
  else if( r > m_RadiusMax )
    {
    double factor = ( r - m_RadiusMax ) / m_RadiusStep;
    double m0 = this->GetKernelMedialness( m_RadiusMax );
    double m1 = this->GetKernelMedialness( m_RadiusMax - m_RadiusStep );
    return m0 - factor * vnl_math_abs(m0 - m1);
    }

  double gfact = -0.5 / ( 0.5 * r * 0.5 * r );

  double pVal = 0;
  double nVal = 0;
  double dwPTot = 0;
  double dwNTot = 0;
  double vTot = 0;
  unsigned int vCount = 0;

  std::vector< double >::iterator iterDist;
  std::vector< double >::iterator iterTanDist;
  std::vector< double >::iterator iterValue;

  double distMax = this->GetKernelExtent() * r;
  double distMin = r - ( distMax - r );
  if( distMin < 0 )
    {
    distMin = 0;
    distMax = 2 * r;
    }

  //std::cout << "Radius = " << r << std::endl;
  iterDist = m_KernelDistances.begin();
  iterTanDist = m_KernelTangentDistances.begin();
  iterValue = m_KernelValues.begin();
  while( iterDist != m_KernelDistances.end() )
    {
    if( ( *iterTanDist ) < 0.5 * r
      && ( *iterDist ) >= distMin && ( *iterDist ) <= distMax )
      {
      vTot += ( *iterValue );
      ++vCount;
      }
    ++iterValue;
    ++iterDist;
    ++iterTanDist;
    }
  double vAvg = 0;
  if( vCount > 0 )
    {
    vAvg = vTot / vCount;
    }

  iterDist = m_KernelDistances.begin();
  iterTanDist = m_KernelTangentDistances.begin();
  iterValue = m_KernelValues.begin();
  while( iterDist != m_KernelDistances.end() )
    {
    if( ( *iterTanDist ) < 0.5 * r
      && ( *iterDist ) >= distMin && ( *iterDist ) <= distMax )
      {
      double w = vcl_exp( gfact * ( r - ( *iterDist ) )
        * ( r - ( *iterDist ) ) );
      double dw = 2 * ( ( *iterDist ) - r ) * gfact * w;
      if( dw > 0 )
        {
        pVal += dw * ( *iterValue - vAvg );
        dwPTot += dw;
        }
      else
        {
        nVal += -dw * ( *iterValue - vAvg );
        dwNTot += -dw;
        }

      //std::cout << *iterDist << ", " << dw << ", " << *iterValue
        //<< std::endl;
      }
    ++iterDist;
    ++iterTanDist;
    ++iterValue;
    }
  if( dwPTot > 0 )
    {
    pVal /= dwPTot;
    }
  if( dwNTot > 0 )
    {
    nVal /= dwNTot;
    }

  double val = ( pVal - nVal );

  //std::cout << "Radius = " << r << "   Value = " << val << std::endl;

  return val;
}

template< class TInputImage >
double
RadiusExtractor2<TInputImage>
::GetKernelBranchness( double r )
{
  double gfact = -0.5 / ( r * r );

  double val = 0;
  double wTot = 0;
  std::vector< double >::iterator iterDist;
  iterDist = m_KernelDistances.begin();
  std::vector< double >::iterator iterValue;
  iterValue = m_KernelValues.begin();
  double distMax = this->GetKernelExtent() * r;
  double distMin = r - ( distMax - r );
  if( distMin < 0 )
    {
    distMin = 0;
    distMax = 2 * r;
    }
  while( iterDist != m_KernelDistances.end() )
    {
    if( ( *iterDist ) > distMin && ( *iterDist ) < distMax )
      {
      double w = vcl_exp( gfact * ( 1.25*r - ( *iterDist ) )
        * ( 1.25*r - ( *iterDist ) ) );
      val += w * ( *iterValue );
      wTot += vnl_math_abs( w );
      }
    ++iterDist;
    ++iterValue;
    }
  if( wTot > 0 )
    {
    val /= wTot;
    }

  return val;
}

template< class TInputImage >
bool
RadiusExtractor2<TInputImage>
::UpdateKernelOptimalRadius( void )
{
  m_MedialnessOptSpline->SetNewData( true );

  m_MedialnessOpt.SetXStep( m_RadiusStep );
  m_MedialnessOpt.SetTolerance( m_RadiusTolerance );
  m_MedialnessOptSpline->SetClip( false );
  m_MedialnessOptSpline->SetXMin( m_RadiusMin / m_RadiusStep );
  m_MedialnessOptSpline->SetXMax( m_RadiusMax / m_RadiusStep );

  double r0 = m_RadiusStart;
  r0 = static_cast<int>( r0 / m_RadiusStep ) * m_RadiusStep;
  int r0Range = static_cast<int>( ( 0.5 * r0 ) / m_RadiusStep );
  if( r0Range < 4 )
    {
    r0Range = 4;
    }
  r0 = r0 - r0Range * m_RadiusStep;
  if( r0 < m_RadiusMin )
    {
    r0 = m_RadiusMin;
    }
  double r0Max = r0;
  double r0MaxMedialness = 0;
  double r0Medialness = this->GetKernelMedialness( r0 );
  double rEnd = r0 + 4 * r0Range * m_RadiusStep;
  while( r0 < rEnd )
    {
    r0 += m_RadiusStep;
    r0Medialness = this->GetKernelMedialness( r0 );
    if( r0Medialness > r0MaxMedialness )
      {
      r0Max = r0;
      r0MaxMedialness = r0Medialness;
      }
    }
  r0 = r0Max;

  /*
  r0 = r0 / m_RadiusStep;
  m_MedialnessOptSpline->Extreme( &r0, &r0Medialness );
  r0 = r0 * m_RadiusStep;
  */

  if( this->GetDebug() )
    {
    std::cout << " cmp: " << r0-m_RadiusStep/2 << " - "
      << m_MedialnessOptSpline->Value( (r0-m_RadiusStep/2) / m_RadiusStep )
      << std::endl;
    std::cout << " cmp: " << r0 << " - "
      << m_MedialnessOptSpline->Value( r0 / m_RadiusStep ) << std::endl;
    std::cout << " cmp: " << r0+m_RadiusStep/2 << " - "
      << m_MedialnessOptSpline->Value( (r0+m_RadiusStep/2) / m_RadiusStep )
      << std::endl;
    }

  if( this->GetDebug() )
    {
    std::cout << "Local extreme at radius r0 = " << r0
      << " with medialness = " << r0Medialness << std::endl;
    std::cout << "  prev radius = " << this->GetRadiusStart()
      << " with medialness = " << m_MedialnessOptSpline->Value(
      this->GetRadiusStart() / m_RadiusStep )
      << std::endl;
    }

  m_KernelOptimalRadius = r0;
  m_KernelOptimalRadiusMedialness = r0Medialness;
  //m_KernelOptimalRadiusBranchness = this->GetKernelBranchness( r0 );

  if( r0Medialness < m_MinMedialnessStart )
    {
    if ( this->GetDebug() )
      {
      std::cout
        << "RadiusExtractor2: calcOptimalScale: kernel fit insufficient"
        << std::endl;
      std::cout << "  Medialness = " << r0Medialness << " < thresh = "
        << m_MinMedialness << std::endl;
      }
    return false;
    }

  return true;
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
    if ( this->GetDebug() )
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
  for( int p = startP; p <= endP; p += m_KernelPointStep  )
    {
    if( p < 0 )
      {
      typename TubeType::PointType p1 =
        tube->GetPoints()[ count ].GetPosition();
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
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( ( p - (tubeSize-1) ) /
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
      tube->GetPoints()[ p ].SetRadius( d * r0 + (1 - d) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m0 + (1 - d) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b0 + (1 - d) * b1 );
      }
    else
      {
      double d = 1;
      if( static_cast< int >( tubePointNum ) != endP )
        {
        d = static_cast< double >( p - tubePointNum )
          / ( endP - tubePointNum );
        }
      tube->GetPoints()[ p ].SetRadius( d * r2 + (1 - d) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m2 + (1 - d) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b2 + (1 - d) * b1 );
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

  os << indent << "MedialnessOpt = " << m_MedialnessOpt << std::endl;
  if( m_MedialnessOptSpline != NULL )
    {
    os << indent << "MedialnessOptSpline = " << m_MedialnessOptSpline
      << std::endl;
    }
  else
    {
    os << indent << "MedialnessOptSpline = NULL" << std::endl;
    }

  os << indent << "RadiusStart = " << m_RadiusStart << std::endl;
  os << indent << "RadiusMin = " << m_RadiusMin << std::endl;
  os << indent << "RadiusMax = " << m_RadiusMax << std::endl;
  os << indent << "RadiusStep = " << m_RadiusStep << std::endl;
  os << indent << "RadiusTolerance = " << m_RadiusTolerance << std::endl;

  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;
  os << indent << "MedialnessFunc = " << m_MedialnessFunc << std::endl;

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

#endif // End !defined(__itktubeRadiusExtractor2_hxx)
