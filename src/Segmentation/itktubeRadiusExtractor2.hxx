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

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itktubeRadiusExtractor2_hxx
#define __itktubeRadiusExtractor2_hxx


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

template< class TInputImage >
class LocalMedialnessSplineValueFunction
: public ::tube::UserFunction< int, double >
{
public:

  LocalMedialnessSplineValueFunction( RadiusExtractor2< TInputImage > *
    _radiusExtractor )
    {
    m_Value = 0;

    m_Spacing = 1;

    m_RadiusExtractor = _radiusExtractor;
    }

  void SetSpacing( double spacing )
    {
    m_Spacing = spacing;
    }

  double GetVariableMultiplier( void )
    {
    return 0.708 * m_Spacing;   // srt(2) / 2 = max distance within a voxel to
                                //              a corner.
    }

  const double & Value( const int & x )
    {
    m_Value = m_RadiusExtractor->GetKernelMedialness(
      x * this->GetVariableMultiplier() );

    return m_Value;
    }

private:

  double   m_Spacing;
  double   m_Value;

  typename RadiusExtractor2< TInputImage >::Pointer   m_RadiusExtractor;
};


/** Constructor */
template< class TInputImage >
RadiusExtractor2<TInputImage>
::RadiusExtractor2( void )
{
  m_InputImage = NULL;

  m_Spacing = 1;
  m_DataMin = 0;
  m_DataMax = -1;

  m_RadiusStart = 0.708;  // All values are in physical space (e.g., mm).
  m_RadiusMin = 0.708/2;
  m_RadiusMax = 15.0;
  m_RadiusStep = 0.708/2;
  m_RadiusTolerance = 0.708/4;

  m_RadiusCorrectionFunction = RADIUS_CORRECTION_NONE;

  m_MinMedialness = 0.15;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.1;

  m_NumKernelPoints = 7;
  m_KernelTube = TubeType::New();
  m_KernelTube->GetPoints().resize(7);

  m_KernelPointStep = 14;
  m_KernelStep = 30;
  m_KernelExtent = 1.25;

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
      ::tube::DebugMessage( "RadiusExtractor2: SetInputImage: Minimum = "
        + std::to_string(m_DataMin) );
      ::tube::DebugMessage( "RadiusExtractor2: SetInputImage: Maximum = "
        + std::to_string(m_DataMax) );
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
  typename TubeType::Pointer tmpTube = TubeType::New();
  tmpTube->SetPoints( points );
  tmpTube->ComputeTangentsAndNormals();

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
  m_KernelTube->GetPoints().resize( m_NumKernelPoints );
}

template< class TInputImage >
void
RadiusExtractor2<TInputImage>
::GenerateKernel( void )
{
  IndexType minXI;
  IndexType maxXI;

  IndexType bufferI;
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    bufferI[i] = static_cast< unsigned int >(
      this->GetKernelExtent() * this->GetRadiusMax() / m_Spacing );
    }

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = m_KernelTube->GetPoints().begin();
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    IndexType kernelPointI;
    m_InputImage->TransformPhysicalPointToIndex( 
      m_KernelTube->GetPoints()[0].GetPositionInObjectSpace(), kernelPointI );
    minXI[i] = static_cast< int >( kernelPointI[i] - bufferI[i] );
    maxXI[i] = static_cast< int >( kernelPointI[i] + bufferI[i] );
    }
  ++pntIter;
  int tempI;
  while( pntIter != m_KernelTube->GetPoints().end() )
    {
    IndexType kernelPointI;
    m_InputImage->TransformPhysicalPointToIndex(
      pntIter->GetPositionInObjectSpace(), kernelPointI );
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      tempI = static_cast< int >( kernelPointI[i] - bufferI[i] );
      if( tempI < minXI[i] )
        {
        minXI[i] = tempI;
        }
      tempI = static_cast< int >( kernelPointI[i] + bufferI[i] );
      if( tempI > maxXI[i] )
        {
        maxXI[i] = tempI;
        }
      }
    ++pntIter;
    }
  unsigned int kernelSize = maxXI[0] - minXI[0] + 1;
  for( unsigned int i = 1; i < ImageDimension; ++i )
    {
    kernelSize *= ( maxXI[i] - minXI[i] + 1 );
    }
  m_KernelValues.resize( kernelSize );
  std::fill( m_KernelValues.begin(), m_KernelValues.end(), 0 );
  m_KernelDistances.resize( kernelSize );
  std::fill( m_KernelDistances.begin(), m_KernelDistances.end(), 0 );
  m_KernelTangentDistances.resize( kernelSize );
  std::fill( m_KernelTangentDistances.begin(), m_KernelTangentDistances.end(),
    0 );

  unsigned int count = 0;
  IndexType xI = minXI;
  bool done = false;
  double thresh = m_DataMin + 0.05 * ( m_DataMax - m_DataMin );
  while( !done )
    {
    if( m_InputImage->GetLargestPossibleRegion().IsInside( xI ) &&
      m_InputImage->GetPixel( xI ) > thresh )
      {
      m_KernelValues[ count ] = ( m_InputImage->GetPixel( xI ) - m_DataMin )
        / ( m_DataMax - m_DataMin );
      if( m_KernelValues[ count ] < 0 )
        {
        m_KernelValues[ count ] = 0;
        }
      else if( m_KernelValues[ count ] > 1 )
        {
        m_KernelValues[ count ] = 1;
        }
      PointType p;
      m_InputImage->TransformIndexToPhysicalPoint( xI, p );

      unsigned int pntCount = 0;
      pntIter = m_KernelTube->GetPoints().begin();

      double pntTangentDist = 0;
      double minTangentDist = 0;
      double minNormalDist = 0;
      int minTangentDistCount = -1;
      while( pntIter != m_KernelTube->GetPoints().end() )
        {
        VectorType pDiff = pntIter->GetPositionInObjectSpace() - p;
        double d1 = 0;
        for( unsigned int i = 0; i < ImageDimension; ++i )
          {
          double tf = pDiff[i] * pntIter->GetTangentInObjectSpace()[i];
          d1 += tf * tf;
          }
        pntTangentDist = std::sqrt( d1 ) * m_Spacing;
        if( pntTangentDist < minTangentDist || minTangentDistCount == -1 )
          {
          minTangentDist = pntTangentDist;
          minTangentDistCount = pntCount;
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
      m_KernelDistances[ count ] = std::sqrt( minNormalDist ) * m_Spacing;
      m_KernelTangentDistances[ count ] = minTangentDist;
      }
    else
      {
      m_KernelValues[ count ] = 0;
      m_KernelDistances[ count ] = this->GetRadiusMax() + 1;
      m_KernelTangentDistances[ count ] = this->GetRadiusMax() + 1;
      }
    ++count;
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

  m_KernelTube->SetPoints( tubePoints );
  m_KernelTube->ComputeTangentsAndNormals();
}

template< class TInputImage >
double
RadiusExtractor2<TInputImage>
::GetKernelMedialness( double r )
{
  // If radius is outside of the legal range, calculate the slope
  // just inside of the legal range and project that rate of decline
  // in medialness outside of the legal range
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
  if( r * this->GetKernelExtent() > distMax )
    {
    distMax = r * this->GetKernelExtent();
    }

  double areaMax = distMax * distMax * vnl_math::pi;
  double areaNeg = areaMax - areaR;

  // Set the areaPos to be equal to areaNeg
  double distMin2 = ( areaR - areaNeg ) / vnl_math::pi;
  double distMin = 0;  // Where areaPos starts
  if( distMin2 > 0 )
    {
    distMin = std::sqrt( distMin2 );  // areaPos should be same size
    }
  double areaMin = distMin * distMin * vnl_math::pi;
  double areaPos = areaR - areaMin;
  if( this->GetDebug() )
    {
    ::tube::InformationMessage( "R = " + std::to_string(r)  );
    ::tube::InformationMessage( "   Dist = " + std::to_string(distMin) + " - "
     + std::to_string(distMax) );
    ::tube::InformationMessage( "   Area = " + std::to_string(areaPos) + " - "
     + std::to_string(areaNeg) );
    }

  const int histoBins = 500;
  double histoPos[histoBins];
  double histoNeg[histoBins];
  double histoPosCount = 0;
  double histoNegCount = 0;
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
    if( ( ( *iterTanDist ) < distMax || ( *iterTanDist ) < 1 )
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
        double distFromR = (((*iterDist) - distMin)/(r - distMin)) / 0.5 + 0.5;
        histoPos[bin] += distFromR;
        histoPosCount += distFromR;
        }
      else
        {
        double distFromR = ((distMax - (*iterDist))/(distMax - r)) / 0.5 + 0.5;
        histoNeg[bin] += distFromR;
        histoNegCount += distFromR;
        }
      }
    ++iterValue;
    ++iterDist;
    ++iterTanDist;
    }

  pVal = 0;
  //std::cout << "HistoCount = " << histoPosCount << " - " << histoNegCount
    //<< std::endl;
  if( histoPosCount > 0 )
    {
    bin = histoBins - 1;
    double binCount = 0;
    double binSum = 0;
    while( binCount < 0.75 * histoPosCount && bin > 0 )
      {
      binCount += histoPos[bin];
      binSum += histoPos[bin] * (bin / (double)histoBins);
      --bin;
      }
    //std::cout << "Pos binSum = " << binSum << "  binCount = " << binCount
      //<< std::endl;
    pVal = binSum / binCount; // ( bin + 0.5 ) / histoBins;
    }
  nVal = 1;
  if( histoNegCount > 0 )
    {
    bin = 0;
    double binCount = 0;
    double binSum = 0;
    while( binCount < 0.75 * histoNegCount && bin < histoBins )
      {
      binCount += histoNeg[bin];
      binSum += histoNeg[bin] * (bin / (double)histoBins);
      ++bin;
      }
    //std::cout << "Neg binSum = " << binSum << "  binCount = " << binCount
      //<< std::endl;
    nVal = binSum / binCount; //( bin - 0.5 ) / histoBins;
    }

  if( this->GetDebug() )
    {
    std::cout <<  "   Count = " << histoPosCount << " - " << histoNegCount
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
  // if r is within RadiusMin, apply a decreasing linear extrapolation
  // from RadiusMin using the slope from just inside RadiusMin
  if( r < m_RadiusMin )
    {
    double factor = ( m_RadiusMin - r ) / m_RadiusTolerance;
    double m0 = this->GetKernelBranchness( m_RadiusMin );
    double m1 = this->GetKernelBranchness( m_RadiusMin
      + m_RadiusTolerance );
    return m0 - factor * std::fabs( m0 - m1 );
    }
  // if r is outside RadiusMax, apply a decreasing linear extrapolation
  // from RadiusMax using the slope from just inside RadiusMin
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
  LocalMedialnessSplineValueFunction< TInputImage > * myFunc = new
    LocalMedialnessSplineValueFunction< TInputImage >( this );
  ::tube::UserFunction< int, double > * myUserFunc = myFunc;

  ::tube::GoldenMeanOptimizer1D * opt = new
    ::tube::GoldenMeanOptimizer1D();
  opt->SetXStep( m_RadiusStep );
  opt->SetTolerance( m_RadiusTolerance );
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
      xMin = ( int )std::ceil( m_RadiusMin );
      xMax = ( int )std::floor( m_RadiusMax );
      x = m_RadiusStart;
      break;
      }
    case RADIUS_CORRECTION_FOR_BINARY_IMAGE:
      {
      xMin = ( int )std::ceil( m_RadiusMin );
      xMax = ( int )std::floor( m_RadiusMax );
      x = m_RadiusStart;
      break;
      }
    case RADIUS_CORRECTION_FOR_CTA:
      {
      xMin = ( int )std::ceil( m_RadiusMin );
      xMax = ( int )std::floor( m_RadiusMax );
      x = m_RadiusStart;
      break;
      }
    case RADIUS_CORRECTION_FOR_MRA:
      {
      xMin = ( int )std::ceil( m_RadiusMin );
      xMax = ( int )std::floor( m_RadiusMax );
      x = m_RadiusStart;
      break;
      }
    };

  myFunc->SetSpacing( m_Spacing );

  xMin /= myFunc->GetVariableMultiplier();
  xMax /= myFunc->GetVariableMultiplier();

  opt->SetXMin( xMin );
  opt->SetXMax( xMax );


  ::tube::SplineApproximation1D * spline = new
  ::tube::SplineApproximation1D( myUserFunc, opt );

  //spline->SetClip( true );
  spline->SetXMin( xMin );
  spline->SetXMax( xMax );

  x /= myFunc->GetVariableMultiplier();
  double xVal = myFunc->Value( x );
  bool result = spline->Extreme( &x, &xVal );
  x *= myFunc->GetVariableMultiplier();

  switch( m_RadiusCorrectionFunction )
    {
    default:
    case RADIUS_CORRECTION_NONE:
      {
      m_KernelOptimalRadius = x;
      break;
      }
    case RADIUS_CORRECTION_FOR_BINARY_IMAGE:
      {
      m_KernelOptimalRadius = ( x * x ) / 24 + 0.5;
      break;
      }
    case RADIUS_CORRECTION_FOR_CTA:
      {
      m_KernelOptimalRadius = x;
      break;
      }
    case RADIUS_CORRECTION_FOR_MRA:
      {
      m_KernelOptimalRadius = x;
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
::ExtractRadii( TubeType * tube, bool verbose )
{
  if( tube->GetPoints().size() == 0 )
    {
    return false;
    }

  tube->RemoveDuplicatePointsInObjectSpace();
  tube->ComputeTangentsAndNormals();

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
    rStart = this->GetKernelOptimalRadius();
    if( verbose )
      {
      std::cout << p << " : r = " << rStart << " -> "
        << tube->GetPoint(p)->GetRadiusInObjectSpace() << std::endl;
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
    rStart = this->GetKernelOptimalRadius();
    if( verbose )
      {
      std::cout << p << " : r = " << rStart << " -> "
        << tube->GetPoint(p)->GetRadiusInObjectSpace() << std::endl;
      }
    }

  if( this->GetDebug() )
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
        tube->GetPoints()[ 0 ].GetPositionInObjectSpace();
      typename TubeType::PointType p2 =
        tube->GetPoints()[ m_NumKernelPoints / 2 *
        m_KernelPointStep ].GetPositionInObjectSpace();
      for( unsigned int i = 0; i < ImageDimension; ++i )
        {
        p2[i] = p1[i] - ( ( p2[i] - p1[i] ) * ( -p / m_KernelPointStep ) );
        }
      m_KernelTube->GetPoints()[ count ].SetPositionInObjectSpace( p2 );
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
      m_KernelTube->GetPoints()[ count ].SetPositionInObjectSpace( p2 );
      }
    else
      {
      m_KernelTube->GetPoints()[ count ] = tube->GetPoints()[ p ];
      }
    ++count;
    }

  m_KernelTube->ComputeTangentsAndNormals();
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
  double r0 = tube->GetPoints()[ startP ].GetRadiusInObjectSpace();
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
    if( p < static_cast< int >( tubePointNum ) )
      {
      double d = 1;
      if( static_cast< int >( tubePointNum ) != startP )
        {
        d = static_cast< double >( tubePointNum - p )
          / ( tubePointNum - startP );
        }
      tube->GetPoints()[ p ].SetRadiusInObjectSpace( d * r0 + ( 1 - d ) * r1 );
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
      tube->GetPoints()[ p ].SetRadiusInObjectSpace( d * r2 + ( 1 - d ) * r1 );
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

  os << indent << "RadiusStart = " << m_RadiusStart << std::endl;
  os << indent << "RadiusMin = " << m_RadiusMin << std::endl;
  os << indent << "RadiusMax = " << m_RadiusMax << std::endl;
  os << indent << "RadiusStep = " << m_RadiusStep << std::endl;
  os << indent << "RadiusTolerance = " << m_RadiusTolerance << std::endl;

  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;

  os << indent << "NumKernelPoints = " << m_NumKernelPoints << std::endl;
  os << indent << "KernelTube = " << m_KernelTube << std::endl;
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
