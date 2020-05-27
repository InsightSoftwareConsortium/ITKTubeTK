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
#ifndef __itktubeRadiusExtractor3_hxx
#define __itktubeRadiusExtractor3_hxx

#include "itktubeRadiusExtractor3.h"

#include "tubeMessage.h"
#include "tubeMatrixMath.h"
#include "tubeTubeMathFilters.h"
#include "tubeUserFunction.h"
#include "tubeGoldenMeanOptimizer1D.h"
#include "tubeSplineApproximation1D.h"

#include <itkMinimumMaximumImageFilter.h>


#include <vnl/vnl_math.h>

namespace itk
{

namespace tube
{

/** Constructor */
template< class TInputImage >
RadiusExtractor3<TInputImage>
::RadiusExtractor3( void )
{
  m_InputImage = NULL;

  m_Spacing = 1;
  m_DataMin = 0;
  m_DataMax = -1;

  m_RadiusStartInIndexSpace = 1.0;  // All values are in index space.
  m_RadiusMinInIndexSpace = 0.708/2;
  m_RadiusMaxInIndexSpace = 4.0;
  m_RadiusStepInIndexSpace = 0.708/2;
  m_RadiusToleranceInIndexSpace = 0.708/3;

  m_RadiusCorrectionScale = 1.0;
  m_RadiusCorrectionFunction = RADIUS_CORRECTION_NONE;

  m_MinMedialness = 0.3;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.15;

  m_NumKernelPoints = 7;
  m_KernelTubePoints.resize( m_NumKernelPoints );

  m_KernelPointStep = 7;
  m_KernelStep = 13;
  m_KernelExtent = 1.75;

  m_KernelValue.clear();
  m_KernelCount.clear();

  m_KernelOptimalRadius = 0;
  m_KernelOptimalRadiusMedialness = 0;
  m_KernelOptimalRadiusBranchness = 0;

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
}

/** Destructor */
template< class TInputImage >
RadiusExtractor3<TInputImage>
::~RadiusExtractor3( void )
{
}

/** Set the input image */
template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::SetInputImage( typename InputImageType::Pointer inputImage )
{
  m_InputImage = inputImage;

  if( m_InputImage )
    {
    typedef MinimumMaximumImageFilter<InputImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter =
      MinMaxFilterType::New();
    minMaxFilter->SetInput( m_InputImage );
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();
    for( unsigned int d=1; d<ImageDimension; ++d )
      {
      if( m_InputImage->GetSpacing()[d] != m_InputImage->GetSpacing()[0] )
        {
        ::tube::WarningMessage(
          "Image is not isotropic. Using x-dim spacing as the spacing." );
        std::cout << "  Spacing = " << m_InputImage->GetSpacing() << std::endl;
        break;
        }
      }
    m_Spacing = m_InputImage->GetSpacing()[0];
    if( this->GetDebug() )
      {
      ::tube::DebugMessage( "RadiusExtractor3: SetInputImage: Minimum = "
        + std::to_string(m_DataMin) );
      ::tube::DebugMessage( "RadiusExtractor3: SetInputImage: Maximum = "
        + std::to_string(m_DataMax) );
      }
    }
}

/** Compute the medialness at a kernel */
template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::GetPointVectorMeasures( std::vector< TubePointType > & points,
  double pntR,
  double & mness,
  double & bness,
  bool doBNess )
{
  ::tube::ComputeVectorTangentsAndNormals( points );

  if( this->GetDebug() )
    {
    ::tube::DebugMessage( "Compute values at point" );
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
RadiusExtractor3<TInputImage>
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

  r0 = this->GetKernelOptimalRadius() * m_Spacing;

  return true;
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::SetNumKernelPoints( unsigned int _numPoints )
{
  m_NumKernelPoints = _numPoints;

  m_KernelTubePoints.resize( m_NumKernelPoints );
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::GenerateKernel( void )
{

  IndexType minXI;
  IndexType maxXI;

  double maxKernelR = this->GetRadiusMaxInIndexSpace();
  double maxKernelDist = this->GetKernelExtent() * maxKernelR;

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = m_KernelTubePoints.begin();
  IndexType kernelPointI;
  m_InputImage->TransformPhysicalPointToIndex( 
    m_KernelTubePoints[0].GetPositionInObjectSpace(), kernelPointI );
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    minXI[i] = static_cast< int >( kernelPointI[i] - (maxKernelDist + 0.5) );
    maxXI[i] = static_cast< int >( kernelPointI[i] + (maxKernelDist + 0.5) );
    }
  ++pntIter;
  int tempI;
  while( pntIter != m_KernelTubePoints.end() )
    {
    IndexType kernelPoint;
    m_InputImage->TransformPhysicalPointToIndex(
      pntIter->GetPositionInObjectSpace(), kernelPoint );
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      tempI = static_cast< int >( kernelPoint[i] - (maxKernelDist + 0.5) );
      if( tempI < minXI[i] )
        {
        minXI[i] = tempI;
        }
      tempI = static_cast< int >( kernelPoint[i] + (maxKernelDist + 0.5) );
      if( tempI > maxXI[i] )
        {
        maxXI[i] = tempI;
        }
      }
    ++pntIter;
    }
  unsigned int kernelSize = static_cast<unsigned int>( maxKernelDist*3 );
  m_KernelValue.resize( kernelSize );
  std::fill( m_KernelValue.begin(), m_KernelValue.end(), 0 );
  m_KernelCount.resize( kernelSize );
  std::fill( m_KernelCount.begin(), m_KernelCount.end(), 0 );

  IndexType xI = minXI;
  bool done = false;
  while( !done )
    {
    if( m_InputImage->GetLargestPossibleRegion().IsInside( xI ) )
      {
      double val = ( m_InputImage->GetPixel( xI ) - m_DataMin )
        / ( m_DataMax - m_DataMin );
      if( val < 0 )
        {
        val = 0;
        }
      else if( val > 1 )
        {
        val = 1;
        }
      PointType p;
      m_InputImage->TransformIndexToPhysicalPoint( xI, p );

      unsigned int pntCount = 0;
      pntIter = m_KernelTubePoints.begin();

      double pntTangentDistI = 0;
      double minTangentDistI = maxKernelR;
      double minNormalDist = -1;
      while( pntIter != m_KernelTubePoints.end() )
        {
        VectorType pDiff = p - pntIter->GetPositionInObjectSpace();
        double d1 = 0;
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          double tf = pDiff[i] * pntIter->GetTangentInObjectSpace()[i];
          d1 += tf * tf;
          }
        pntTangentDistI = std::sqrt( d1 ) / m_Spacing;
        if( pntTangentDistI < minTangentDistI )
          {
          minTangentDistI = pntTangentDistI;
          d1 = 0;
          for( unsigned int i = 0; i < ImageDimension; ++i )
            {
            double tf = pDiff[i] * pntIter->GetNormal1InObjectSpace()[i];
            d1 += tf * tf;
            }
          minNormalDist = d1;
          if( ImageDimension == 3 )
            {
            double d2 = 0;
            for( unsigned int i = 0; i < ImageDimension; ++i )
              {
              double tf = pDiff[i] * pntIter->GetNormal2InObjectSpace()[i];
              d2 += tf * tf;
              }
            minNormalDist += d2;
            }
          }
        ++pntIter;
        ++pntCount;
        }
      if( minNormalDist != -1 )
        {
        double distI = std::sqrt( minNormalDist ) / m_Spacing;
        double count = (distI / maxKernelDist) * kernelSize;
        // std::cout << distI << " : " << count << " : " << val << std::endl;
        if( count < 0 )
          {
          count = 0;
          }
        if( count < kernelSize )
          {
          m_KernelValue[ count ] += val;
          m_KernelCount[ count ]++;
          }
        }
      }
    unsigned int d = 0;
    while( d < ImageDimension && ++xI[d] > maxXI[d] )
      {
      xI[d] = minXI[d];
      ++d;
      }
    if( d >= ImageDimension )
      {
      done = true;
      }
    }
  for(unsigned int i=0; i<kernelSize; ++i )
    {
    if( m_KernelCount[i] > 0 )
      {
      m_KernelValue[i] /= m_KernelCount[i];
      }
    std::cout << m_KernelValue[i] << "(" << m_KernelCount[i] << ") ";
    }
  std::cout << std::endl;
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
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
RadiusExtractor3<TInputImage>
::GetKernelMedialness( double r )
{
  int kernelSize = m_KernelValue.size();
  double maxKernelR = this->GetRadiusMaxInIndexSpace();
  double maxKernelDist = this->GetKernelExtent() * maxKernelR;

  double maxV = 0;
  double minV = 1;
  for( unsigned int i = 1; i<kernelSize-1; ++i )
    {
    double val = 0;
    unsigned int count = 0;
    for( unsigned int j = i-1; j<i+1; ++j )
      {
      if( m_KernelCount[i] > 0 )
        {
        val += m_KernelValue[i];
        ++count;
        }
      }
    if( count > 0 )
      {
      val /= count;
      }
    if( val > maxV )
      {
      maxV = val;
      }
    if( val < minV )
      {
      minV = val;
      }
    }

  double medialness = 0;
  if( maxV > minV )
    {
    medialness = minV + (maxV - minV) / 2;
    }

  return medialness;
}

template< class TInputImage >
double
RadiusExtractor3<TInputImage>
::GetKernelBranchness( double r )
{
  return 0;
}

template< class TInputImage >
bool
RadiusExtractor3<TInputImage>
::UpdateKernelOptimalRadius( void )
{
  unsigned int kernelSize = m_KernelValue.size();
  double maxKernelR = this->GetRadiusMaxInIndexSpace();
  double maxKernelDist = this->GetKernelExtent() * maxKernelR;

  double maxV = 0;
  double minV = 1;
  unsigned int maxI = 0;
  unsigned int minI = kernelSize-1;
  for( unsigned int i = 1; i<kernelSize-1; ++i )
    {
    double val = 0;
    unsigned int count = 0;
    for( unsigned int j = i-1; j<i+1; ++j )
      {
      if( m_KernelCount[j] > 0 )
        {
        val += m_KernelValue[j];
        ++count;
        }
      }
    if( count > 0 )
      {
      val /= count;
      if( val < minV)
        {
        minV = val;
        minI = i;
        }
      }
    }
  for( unsigned int i = minI-1; i>0; --i )
    {
    double val = 0;
    unsigned int count = 0;
    for( unsigned int j = i-1; j<i+1; ++j )
      {
      if( m_KernelCount[j] > 0 )
        {
        val += m_KernelValue[j];
        ++count;
        }
      }
    if( count > 0 )
      {
      val /= count;
      if( val > maxV )
        {
        maxV = val;
        maxI = i;
        }
      }
    }

  double thresh = 0;
  if( maxV > minV )
    {
    thresh = minV + (maxV - minV) / 2;
    for(unsigned int i = maxI+1; i<kernelSize-1; ++i )
      {
      int count = 0;
      double val = 0;
      for( unsigned int j = i-1; j<i+1; ++j )
        {
        if( m_KernelCount[j] > 0 )
          {
          val += m_KernelValue[j];
          ++count;
          }
        }
      if( count > 0 )
        {
        val /= count;
        if( val < thresh )
          {
          m_KernelOptimalRadius = ( (i-0.5) / (double)(kernelSize) ) * maxKernelDist;
          break;
          }
        }
      }
    }
  else
    {
    m_KernelOptimalRadius = 1;
    }
  std::cout << "   " << thresh << " : " << m_KernelOptimalRadius << std::endl;

  switch( m_RadiusCorrectionFunction )
    {
    default:
    case RADIUS_CORRECTION_NONE:
      {
      break;
      }
    case RADIUS_CORRECTION_FOR_BINARY_IMAGE:
      {
      m_KernelOptimalRadius = ( m_KernelOptimalRadius * m_KernelOptimalRadius ) / 24 + 0.5;
      break;
      }
    case RADIUS_CORRECTION_FOR_CTA:
      {
      m_KernelOptimalRadius = m_KernelOptimalRadius;
      break;
      }
    case RADIUS_CORRECTION_FOR_MRA:
      {
      m_KernelOptimalRadius = m_KernelOptimalRadius;
      break;
      }
    };

  m_KernelOptimalRadiusMedialness = thresh;

  return true;
}

template< class TInputImage >
bool
RadiusExtractor3<TInputImage>
::ExtractRadii( TubeType * tube, bool verbose )
{
  unsigned int tubeSize = tube->GetPoints().size();
  if( tubeSize < m_NumKernelPoints * m_KernelPointStep )
    {
    return false;
    }

  tube->RemoveDuplicatePointsInObjectSpace();
  ::tube::ComputeVectorTangentsAndNormals< TubePointType >(
    tube->GetPoints() );

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = tube->GetPoints().begin();
  while( pntIter != tube->GetPoints().end() )
    {
    pntIter->SetRadiusInObjectSpace( 0 );
    ++pntIter;
    }

  int pntCount = 0;
  pntIter = tube->GetPoints().begin();
  while( pntIter != tube->GetPoints().end() && pntIter->GetId() != 0 )
    {
    ++pntIter;
    ++pntCount;
    }

  if( pntIter == tube->GetPoints().end() )
    {
    if( this->GetDebug() )
      {
      ::tube::WarningMessage(
        "Warning: PointID 0 not found. Using mid-point of tube." );
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
    ::tube::DebugMessage( "Found point " +
      std::to_string( ( *pntIter ).GetId() ) );
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
    rStart = this->GetKernelOptimalRadius() * m_Spacing;
    if( verbose )
      {
      std::cout << p << " : r = " << rStart << std::endl;
      }
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
    rStart = this->GetKernelOptimalRadius() * m_Spacing;
    if( verbose )
      {
      std::cout << p << " : r = " << rStart << std::endl;
      }
    }

  //if( this->GetDebug() )
    {
    ::tube::DebugMessage( "Radius results:" );
    pntIter = tube->GetPoints().begin();
    while( pntIter != tube->GetPoints().end() )
      {
      ::tube::DebugMessage( "   " + std::to_string(pntIter->GetId()) + " : "
        + std::to_string(pntIter->GetRadiusInObjectSpace()) ); 
      ++pntIter;
      }
    }

  return true;
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
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
        tube->GetPoints()[ 0 ].GetPositionInObjectSpace();
      typename TubeType::PointType p2 =
        tube->GetPoints()[ m_NumKernelPoints / 2 *
        m_KernelPointStep ].GetPositionInObjectSpace();
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( -p / m_KernelPointStep ) );
        }
      m_KernelTubePoints[ count ].SetPositionInObjectSpace( p2 );
      }
    else if( p > static_cast< int >( tubeSize ) - 1 )
      {
      typename TubeType::PointType p1 =
        tube->GetPoints()[ tubeSize - 1 ].GetPositionInObjectSpace();
      typename TubeType::PointType p2 =
        tube->GetPoints()[ tubeSize - 1 - m_NumKernelPoints / 2 *
        m_KernelPointStep ].GetPositionInObjectSpace();
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( ( p - ( tubeSize-1 ) ) /
            m_KernelPointStep ) );
        }
      m_KernelTubePoints[ count ].SetPositionInObjectSpace( p2 );
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
RadiusExtractor3<TInputImage>
::RecordOptimaAtTubePoints( unsigned int tubePointNum,
  TubeType * tube )
{
  int tubeSize = tube->GetPoints().size();

  int midNum = static_cast<int>(tubePointNum);

  double r1 = this->GetKernelOptimalRadius() * m_Spacing;
  double m1 = this->GetKernelOptimalRadiusMedialness();
  double b1 = this->GetKernelOptimalRadiusBranchness();

  int startP = midNum
    - ( m_NumKernelPoints / 2 ) * m_KernelPointStep - 1;
  if( startP < 0 )
    {
    startP = 0;
    }
  double r0 = tube->GetPoints()[ startP ].GetRadiusInObjectSpace();
  double m0 = tube->GetPoints()[ startP ].GetMedialness();
  double b0 = tube->GetPoints()[ startP ].GetBranchness();
  if( r0 == 0 )
    {
    r0 = r1;
    m0 = m1;
    b0 = b1;
    }

  int endP = startP + ( m_NumKernelPoints ) * m_KernelPointStep + 1;
  if( endP > tubeSize-1 )
    {
    endP = tubeSize-1;
    }
  double r2 = tube->GetPoints()[ endP ].GetRadiusInObjectSpace();
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
    if( p < midNum )
      {
      double d = 0;
      if( midNum != startP )
        {
        d = static_cast< double >( midNum - p ) / ( midNum - startP );
        }
      tube->GetPoints()[ p ].SetRadiusInObjectSpace( d * r0 + ( 1 - d ) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m0 + ( 1 - d ) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b0 + ( 1 - d ) * b1 );
      }
    else
      {
      double d = 0;
      if( midNum != endP )
        {
        d = static_cast< double >( p - midNum ) / ( endP - midNum );
        }
      tube->GetPoints()[ p ].SetRadiusInObjectSpace( d * r2 + ( 1 - d ) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m2 + ( 1 - d ) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b2 + ( 1 - d ) * b1 );
      }
    }
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_InputImage.IsNotNull() )
    {
    os << indent << "InputImage = " << m_InputImage << std::endl;
    }
  else
    {
    os << indent << "InputImage = NULL" << std::endl;
    }
  os << indent << "Spacing = " << m_Spacing << std::endl;
  os << indent << "DataMin = " << m_DataMin << std::endl;
  os << indent << "DataMax = " << m_DataMax << std::endl;

  os << indent << "RadiusStartInIndexSpace = "
    << m_RadiusStartInIndexSpace << std::endl;
  os << indent << "RadiusMinInIndexSpace = "
    << m_RadiusMinInIndexSpace << std::endl;
  os << indent << "RadiusMaxInIndexSpace = "
    << m_RadiusMaxInIndexSpace << std::endl;
  os << indent << "RadiusStepInIndexSpace = "
    << m_RadiusStepInIndexSpace << std::endl;
  os << indent << "RadiusToleranceInIndexSpace = "
    << m_RadiusToleranceInIndexSpace << std::endl;

  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;

  os << indent << "NumKernelPoints = " << m_NumKernelPoints << std::endl;
  os << indent << "KernelTubePoints = " << m_KernelTubePoints.size()
    << std::endl;
  os << indent << "KernelPointStep = " << m_KernelPointStep << std::endl;
  os << indent << "KernelStep = " << m_KernelStep << std::endl;
  os << indent << "KernelExtent = " << m_KernelExtent << std::endl;
  os << indent << "KernelValue = " << m_KernelValue.size() << std::endl;
  os << indent << "KernelCount = " << m_KernelCount.size()
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
RadiusExtractor3<TInputImage>
::SetIdleCallBack( bool ( *idleCallBack )() )
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Status Call back */
template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::SetStatusCallBack( void ( *statusCallBack )( const char *, const char *,
    int ) )
{
  m_StatusCallBack = statusCallBack;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeRadiusExtractor3_hxx )
