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

#include <algorithm>
#include <math.h>

#include "itktubeRadiusExtractor3.h"

#include "tubeMessage.h"
#include "tubeMatrixMath.h"
#include "tubeTubeMathFilters.h"

#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkFRPROptimizer.h"

#include "itkMinimumMaximumImageFilter.h"

#include <vnl/vnl_math.h>

namespace itk
{

class ProfileCurve : public SingleValuedCostFunction
{
public:
  ITK_DISALLOW_COPY_AND_MOVE(ProfileCurve);

  using Self = ProfileCurve;
  using Superclass = SingleValuedCostFunction;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  itkNewMacro(Self);

  itkTypeMacro(ProfileCurve, SingleValuedCostFunction);

  typedef SingleValuedCostFunction::ParametersType ParametersType;
  typedef SingleValuedCostFunction::DerivativeType DerivativeType;

  void SetData( std::vector<double> * data )
    {
    m_Data = data;
    }

  double GetValue( const ParametersType & p ) const override
    {
    double err = 0;
    for( unsigned int i=0; i<m_Data->size(); ++i )
      {
      double tf = (*m_Data)[i] - (p[0]-(p[1]/(1+exp(-p[2]*(i-p[3])))));
      if( isnan(tf) )
        {
        tf = 1;
        }
      err += tf * tf;
      }
    err = std::sqrt(err/m_Data->size());
    return err;
    }

  void GetDerivative( const ParametersType & p, DerivativeType & d ) const override
    {
    d[0] = 0;
    d[1] = 0;
    d[2] = 0;
    d[3] = 0;
    double v = this->GetValue(p);
    for( unsigned int i=0; i<m_Data->size(); ++i )
      {
      double expV = exp(-p[2]*(i-p[3]));
      double denom = 1 + expV;
      d[0] += - 2 * p[0] + (2 * p[1])/denom;
      d[1] +=  2 * (p[0]-(p[1]/denom)) / denom;
      d[2] += -(2*p[1]*(p[0]-(p[1]/denom))*(p[3]-i)*expV)/(denom*denom);
      d[3] += -(2*p[1]*p[2]*(p[0]-(p[1]/denom))*expV)/(denom*denom);
      }
    d[0] /= m_Data->size();
    d[1] /= m_Data->size();
    d[2] /= m_Data->size();
    d[3] /= m_Data->size();
    for( unsigned int i=0; i<4; ++i )
      {
      if( isnan(d[i]) )
        {
        d[i] = 0;
        }
      }
    }

  unsigned int GetNumberOfParameters(void) const override
    {
      return 4;
    }

protected:
  ProfileCurve() {};
  ~ProfileCurve() override = default;

  std::vector<double> * m_Data;
};

namespace tube
{

/** Constructor */
template< class TInputImage >
RadiusExtractor3<TInputImage>
::RadiusExtractor3( void )
{
  m_InputImage = nullptr;

  m_Spacing = 1;
  m_DataMin = 0;
  m_DataMax = -1;

  m_RadiusStartInIndexSpace = 0.75;  // All values are in index space.
  m_RadiusMinInIndexSpace = 0.708/2.0;
  m_RadiusMaxInIndexSpace = 8.0;

  m_MinMedialness = 0.10;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.10;

  m_KernelNumberOfPoints = 5;
  m_KernelPointStep = 12;
  m_KernelStep = 17;

  m_KernelTube = TubeType::New();
  m_KernelTube->GetPoints().resize(m_KernelNumberOfPoints);

  m_ProfileNumberOfBins = 20;
  m_ProfileBinValue.resize( m_ProfileNumberOfBins );
  m_ProfileBinCount.resize( m_ProfileNumberOfBins );
  std::fill( m_ProfileBinValue.begin(), m_ProfileBinValue.end(), 0 );
  std::fill( m_ProfileBinCount.begin(), m_ProfileBinCount.end(), 0 );

  m_KernelOptimalRadius = 0;
  m_KernelOptimalRadiusMedialness = 0;
  m_KernelOptimalRadiusBranchness = 0;

  m_IdleCallBack = nullptr;
  m_StatusCallBack = nullptr;
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

  if( this->GetDebug() )
    {
    ::tube::DebugMessage( "Compute values at point" );
    }

  unsigned int tempNumPoints = this->GetKernelNumberOfPoints();

  unsigned int numPoints = points.size();
  this->SetKernelNumberOfPoints( numPoints );
  this->SetKernelTubePoints( points );

  this->GenerateKernelProfile();

  mness = this->GetKernelMedialness( pntR );
  if( doBNess )
    {
    bness = this->GetKernelBranchness( pntR );
    }

  this->SetKernelNumberOfPoints( tempNumPoints );
}

/** Compute the Optimal scale */
template< class TInputImage >
bool
RadiusExtractor3<TInputImage>
::GetPointVectorOptimalRadius( std::vector< TubePointType > & points,
  double & r0,
  double rMin,
  double rMax )
{
  unsigned int tempNumPoints = this->GetKernelNumberOfPoints();

  unsigned int numPoints = points.size();
  this->SetKernelNumberOfPoints( numPoints );
  this->SetKernelTubePoints( points );

  double tempXStart = this->GetRadiusStart();
  this->SetRadiusStart( r0 );
  double tempXMin = this->GetRadiusMin();
  this->SetRadiusMin( rMin );
  double tempXMax = this->GetRadiusMax();
  this->SetRadiusMax( rMax );

  this->GenerateKernelProfile();

  this->OptimizeKernelRadius();

  this->SetRadiusStart( tempXStart );
  this->SetRadiusMin( tempXMin );
  this->SetRadiusMax( tempXMax );

  this->SetKernelNumberOfPoints( tempNumPoints );

  r0 = this->GetKernelOptimalRadius();

  return true;
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::SetKernelNumberOfPoints( unsigned int _numPoints )
{
  m_KernelNumberOfPoints = _numPoints;
  m_KernelTube->GetPoints().resize( m_KernelNumberOfPoints );
}

template< class TInputImage >
double
RadiusExtractor3<TInputImage>
::GetProfileMaxDistance()
{
  double maxR = this->GetRadiusMax() - this->GetRadiusMin();

  double profileMaxDistance = maxR * pow(m_ProfileNumberOfBins, 1.6)
    / pow(m_ProfileNumberOfBins-2, 1.6);

  profileMaxDistance += this->GetRadiusMin();

  return profileMaxDistance;
}

template< class TInputImage >
double
RadiusExtractor3<TInputImage>
::GetProfileBinRadius( double i )
{
  i = fabs(i);

  double maxR = this->GetRadiusMax() - this->GetRadiusMin();
  double profileMaxDistance = this->GetProfileMaxDistance();

  double x = pow( i, 1.6 ) / pow( m_ProfileNumberOfBins, 1.6 )
    * profileMaxDistance;

  x += this->GetRadiusMin();

  return x;
}

template< class TInputImage >
double
RadiusExtractor3<TInputImage>
::GetProfileBinNumber( double x )
{
  x = fabs(x);

  double maxR = this->GetRadiusMax() - this->GetRadiusMin();
  double profileMaxDistance = this->GetProfileMaxDistance();

  x -= this->GetRadiusMin();

  double i = pow( ( x * pow(m_ProfileNumberOfBins, 1.6 ) / profileMaxDistance ),
    ( 1.0 / 1.6) ); 

  return i;
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::GenerateKernelProfile( void )
{
  IndexType minXIndex;
  IndexType maxXIndex;

  double maxR = this->GetRadiusMax();
  double profileMaxDistance = this->GetProfileMaxDistance();
  double profileMaxIndex = profileMaxDistance/m_Spacing;

  typename std::vector< TubePointType >::iterator pntIter;
  pntIter = m_KernelTube->GetPoints().begin();
  PointType p = pntIter->GetPositionInObjectSpace();
  IndexType kernelPointIndex;
  m_InputImage->TransformPhysicalPointToIndex( p, kernelPointIndex );
  for( unsigned int i = 0; i < ImageDimension; ++i )
    {
    minXIndex[i] = static_cast< int >( kernelPointIndex[i] -
      (profileMaxIndex + 0.5) );
    maxXIndex[i] = static_cast< int >( kernelPointIndex[i] +
      (profileMaxIndex + 0.5) );
    }
  ++pntIter;
  int tempIndex;
  while( pntIter != m_KernelTube->GetPoints().end() )
    {
    p = pntIter->GetPositionInObjectSpace();
    m_InputImage->TransformPhysicalPointToIndex( p, kernelPointIndex );
    for( unsigned int i = 0; i < ImageDimension; ++i )
      {
      tempIndex = static_cast< int >( kernelPointIndex[i] -
        (profileMaxIndex + 0.5) );
      if( tempIndex < minXIndex[i] )
        {
        minXIndex[i] = tempIndex;
        }
      tempIndex = static_cast< int >( kernelPointIndex[i] +
        (profileMaxIndex + 0.5) );
      if( tempIndex > maxXIndex[i] )
        {
        maxXIndex[i] = tempIndex;
        }
      }
    ++pntIter;
    }

  std::fill( m_ProfileBinValue.begin(), m_ProfileBinValue.end(), 0 );
  std::fill( m_ProfileBinCount.begin(), m_ProfileBinCount.end(), 0 );
  IndexType xIndex = minXIndex;
  bool done = false;
  while( !done )
    {
    if( m_InputImage->GetLargestPossibleRegion().IsInside( xIndex ) )
      {
      double val = ( m_InputImage->GetPixel( xIndex ) - m_DataMin )
        / ( m_DataMax - m_DataMin );
      if( val >= 0 && val < 1 )
        {
        PointType p;
        m_InputImage->TransformIndexToPhysicalPoint( xIndex, p );
  
        unsigned int pntCount = 0;
        pntIter = m_KernelTube->GetPoints().begin();
  
        double pntTangentDistance = 0;
        double minTangentDistance = 2*m_Spacing;
        double minNormalDistance = -1;
        typename std::vector< TubePointType >::iterator minTangentPnt;
        minTangentPnt = m_KernelTube->GetPoints().end();
        while( pntIter != m_KernelTube->GetPoints().end() )
          {
          VectorType pDiff = p - pntIter->GetPositionInObjectSpace();
          double d1 = 0;
          for( unsigned int i = 0; i < ImageDimension; ++i )
            {
            double tf = pDiff[i] * pntIter->GetTangentInObjectSpace()[i];
            d1 += tf * tf;
            }
          pntTangentDistance = std::sqrt( d1 );
          if( pntTangentDistance < minTangentDistance )
            {
            minTangentDistance = pntTangentDistance;
            minTangentPnt = pntIter;
            }
          ++pntIter;
          ++pntCount;
          }
        if( minTangentPnt != m_KernelTube->GetPoints().end())
          {
          double d1 = 0;
          VectorType pDiff = p - minTangentPnt->GetPositionInObjectSpace();
          for( unsigned int i = 0; i < ImageDimension; ++i )
            {
            double tf = pDiff[i] * minTangentPnt->GetNormal1InObjectSpace()[i];
            d1 += tf * tf;
            }
          minNormalDistance = d1;
          if( ImageDimension == 3 )
            {
            double d2 = 0;
            for( unsigned int i = 0; i < ImageDimension; ++i )
              {
              double tf = pDiff[i] * minTangentPnt->GetNormal2InObjectSpace()[i];
              d2 += tf * tf;
              }
            minNormalDistance += d2;
            }
  
          double dist = std::sqrt( minNormalDistance );
          double bin = this->GetProfileBinNumber( dist );
          if( bin >= 0 && bin < static_cast<int>(m_ProfileNumberOfBins) )
            {
            m_ProfileBinValue[ (int)bin ] += val;
            m_ProfileBinCount[ (int)bin ]++;
            if( bin > 0 )
              {
              m_ProfileBinValue[ bin-1 ] += 0.5 * val * (1-(bin-(int)bin));
              m_ProfileBinCount[ bin-1 ] += 0.5 * (1-(bin-(int)bin));
              }
            if( bin < static_cast<int>(m_ProfileNumberOfBins)-1 )
              {
              m_ProfileBinValue[ bin+1 ] += 0.5 * val * (bin-(int)bin);
              m_ProfileBinCount[ bin+1 ] += 0.5 * (bin-(int)bin);
              }
            }
          }
        }
      }
    unsigned int d = 0;
    while( d < ImageDimension && ++xIndex[d] > maxXIndex[d] )
      {
      xIndex[d] = minXIndex[d];
      ++d;
      }
    if( d >= ImageDimension )
      {
      done = true;
      }
    }
  for(unsigned int i=0; i<m_ProfileNumberOfBins; ++i )
    {
    if( m_ProfileBinCount[i] > 0 && m_ProfileBinValue[i] > 0 )
      {
      m_ProfileBinValue[i] /= m_ProfileBinCount[i];
      }
    else
      {
      if( i>0 )
        {
        m_ProfileBinValue[i] = m_ProfileBinValue[i-1];
        }
      }
    }
  int i=0;
  while(i<static_cast<int>(m_ProfileNumberOfBins) &&
    m_ProfileBinValue[i]<=m_ProfileBinValue[i+1])
    {
    ++i;
    }
  while(i>0)
    {
    m_ProfileBinValue[i-1] = m_ProfileBinValue[i];
    --i;
    }
  i = m_ProfileNumberOfBins-1;
  while(i>0 && m_ProfileBinValue[i]>=m_ProfileBinValue[i-1])
    {
    --i;
    }
  while(i<static_cast<int>(m_ProfileNumberOfBins)-1)
    {
    m_ProfileBinValue[i+1] = m_ProfileBinValue[i];
    ++i;
    }
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::SetKernelTubePoints( const std::vector< TubePointType > & tubePoints )
{
  if( tubePoints.size() != m_KernelNumberOfPoints )
    {
    std::cerr << "Error: number of kernel points not equal to expected."
      << std::endl;
    std::cerr << "   TubePointsSize = " << tubePoints.size() << std::endl;
    std::cerr << "   KernelNumberOfPoints = " << m_KernelNumberOfPoints
      << std::endl;
    }

  m_KernelTube->SetPoints( tubePoints );
  if( tubePoints.size() > 1 )
    {
    m_KernelTube->ComputeTangentsAndNormals();
    }

  if( tubePoints.size() == 1 )
    {
    auto kernPnt = m_KernelTube->GetPoints().begin();
    VectorType v;
    v.Fill(0);
    v[0] = 1;
    CovariantVector<double, ImageDimension> cv;
    cv.Fill(0);
    cv[1] = 1;
    double sum = 0;
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      sum += std::abs(kernPnt->GetTangentInObjectSpace()[i]);
      }
    if( sum == 0 )
      {
      sum = 0;
      for( unsigned int i=0; i<ImageDimension; ++i )
        {
        sum += std::abs(kernPnt->GetNormal1InObjectSpace()[i]);
        }
      if( sum == 0 )
        {
        std::cout << "WARNING: Single point kernel, setting tangent and normals."
          << std::endl;
        kernPnt->SetTangentInObjectSpace(v);
        kernPnt->SetNormal1InObjectSpace(cv);
        if( ImageDimension > 2 )
          {
          cv.Fill(0);
          cv[2] = 1;
          kernPnt->SetNormal2InObjectSpace(cv);
          }
        }
      else
        {
        std::cout << "WARNING: Single point kernel, setting tangent."
          << std::endl;
        kernPnt->SetTangentInObjectSpace(v);
        }
      }
    sum = 0;
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      sum += std::abs(kernPnt->GetNormal1InObjectSpace()[i]);
      }
    if( sum == 0 )
      {
      std::cout << "WARNING: Single point kernel, resetting normal 1"
        << std::endl;
      kernPnt->SetNormal1InObjectSpace(cv);
      }
    if( ImageDimension > 2 )
      {
      sum = 0;
      for( unsigned int i=0; i<ImageDimension; ++i )
        {
        sum += std::abs(kernPnt->GetNormal2InObjectSpace()[i]);
        }
      if( sum == 0 )
        {
        std::cout << "WARNING: Single point kernel, resetting normal 2"
          << std::endl;
        kernPnt->SetNormal2InObjectSpace(cv);
        }
      }
    }
}

template< class TInputImage >
bool
RadiusExtractor3<TInputImage>
::OptimizeKernelRadius( void )
{
  double rMin = this->GetRadiusMin();
  double rMax = this->GetRadiusMax();

  m_KernelOptimalRadius = this->GetRadiusStart();

  typedef FRPROptimizer OptimizerType;
  OptimizerType::Pointer opt = OptimizerType::New();

  itk::ProfileCurve::Pointer curvFunc = itk::ProfileCurve::New();
  curvFunc->SetData( &m_ProfileBinValue );

  OptimizerType::ParametersType params;
  params.SetSize(4);
  params[0] = ( m_ProfileBinValue[0] + m_ProfileBinValue[1] ) / 2;
  params[1] = params[0] - (m_ProfileBinValue[m_ProfileNumberOfBins-2]
    + m_ProfileBinValue[m_ProfileNumberOfBins-1])/2;
  params[2] = 1;
  params[3] = this->GetProfileBinNumber(m_KernelOptimalRadius);

  OptimizerType::ScalesType scales;
  scales.SetSize(4);
  scales[0] = 10.0;
  scales[1] = 10.0;
  scales[2] = 0.8;
  scales[3] = 0.001;

  opt->SetCostFunction( curvFunc.GetPointer() );
  opt->SetScales( scales );
  opt->SetInitialPosition( params );
  opt->SetUseUnitLengthGradient(true);
  opt->SetStepLength(1.0);
  opt->SetCatchGetValueException( true );
  opt->SetMaximumIteration(200);
  opt->SetMaximumLineIteration(100);
  opt->SetStepTolerance(0.01);
  opt->StartOptimization();

  params = opt->GetCurrentPosition();

  m_KernelOptimalRadius = this->GetProfileBinRadius(params[3]);
  m_KernelOptimalRadiusMedialness = params[1];
  m_KernelOptimalRadiusBranchness = params[2];

  if( this->GetKernelOptimalRadiusMedialness() < m_MinMedialness )
    {
    m_KernelOptimalRadius = ( m_KernelOptimalRadius + this->GetRadiusStart() )
      / 2.0;
    if(this->GetDebug())
      {
      std::cout << "r = " << m_KernelOptimalRadius << " : Medialness Limit = "
        << m_KernelOptimalRadiusMedialness << std::endl;
      }
    }

  if(m_KernelOptimalRadius<this->GetRadiusMin())
  {
    m_KernelOptimalRadius = this->GetRadiusMin();
  }
  else if(m_KernelOptimalRadius>this->GetRadiusMax())
  {
    m_KernelOptimalRadius = this->GetRadiusMax();
  }

  if( this->GetDebug() )
    {
    std::cout << "Params = " << params << std::endl;
    std::cout << "............ Kernel = ";
    for(unsigned int i=0; i<m_ProfileNumberOfBins; ++i )
      {
      std::cout << "   " << m_ProfileBinValue[i] << " (" << m_ProfileBinCount[i]
        << ")" << std::endl;
      }
    std::cout << std::endl;
    }

  return true;
}

template< class TInputImage >
bool
RadiusExtractor3<TInputImage>
::ExtractRadii( TubeType * tube, bool verbose )
{
  unsigned int tubeSize = tube->GetPoints().size();
  if( tubeSize < m_KernelNumberOfPoints * m_KernelPointStep )
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
    pntCount = 0;
    pntIter = tube->GetPoints().begin();
    unsigned int psize = tube->GetPoints().size();
    for( unsigned int i=0; i<psize/2; i++ )
      {
      ++pntIter;
      ++pntCount;
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
    this->GenerateKernelProfile();
    this->OptimizeKernelRadius();
    this->RecordOptimaAtTubePoints( p, tube );
    rStart = (rStart0 + m_KernelOptimalRadius) / 2;
    if( verbose )
      {
      std::cout << p << " : x = "
        << tube->GetPoints()[p].GetPositionInObjectSpace()
        << " : r = " << tube->GetPoints()[p].GetRadiusInObjectSpace()
        << std::endl;
      if( this->GetDebug() )
        {
        std::cout << "............ Kernel = ";
        for(unsigned int i=0; i<m_ProfileNumberOfBins; ++i )
          {
          std::cout << this->GetProfileBinRadius(i)
            << "   " << m_ProfileBinValue[i]
            << " (" << m_ProfileBinCount[i] << ")" << std::endl;
          }
        std::cout << std::endl;
        }
      }
    }

  rStart = rStart0;
  for( int p = static_cast< int >( pntCount ) + this->GetKernelStep()/2; p >= 0;
    p -= this->GetKernelStep() )
    {
    this->SetRadiusStart( rStart );
    this->GenerateKernelTubePoints( p, tube );
    this->GenerateKernelProfile();
    this->OptimizeKernelRadius();
    this->RecordOptimaAtTubePoints( p, tube );
    rStart = (rStart0 + m_KernelOptimalRadius) / 2;
    if( verbose )
      {
      std::cout << p << " : x = "
        << tube->GetPoints()[p].GetPositionInObjectSpace()
        << " : r = " << tube->GetPoints()[p].GetRadiusInObjectSpace()
        << std::endl;
      if(this->GetDebug())
        {
        std::cout << "............ Kernel = ";
        for(unsigned int i=0; i<m_ProfileNumberOfBins; ++i )
          {
          std::cout << "   " << m_ProfileBinValue[i]
            << " (" << m_ProfileBinCount[i] << ")" << std::endl;
          }
        std::cout << std::endl;
        }
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
RadiusExtractor3<TInputImage>
::GenerateKernelTubePoints( unsigned int tubePointNum,
  TubeType * tube )
{
  unsigned int tubeSize = tube->GetPoints().size();
  if( tubeSize < m_KernelNumberOfPoints * m_KernelPointStep )
    {
    std::cerr << "RadiusExtractor: Tube length is too short" << std::endl;
    return;
    }

  int startP = tubePointNum - ( (m_KernelNumberOfPoints-1) / 2 ) * m_KernelPointStep;
  int endP = startP + ( m_KernelNumberOfPoints - 1 ) * m_KernelPointStep;

  if( startP < 0 || endP >= static_cast<int>(tubeSize))
  {
    if( startP < 0)
    {
      startP = 0;
      endP = startP + ( m_KernelNumberOfPoints - 1 ) * m_KernelPointStep;
    }
    else
    {
      endP = tubeSize-1;
      startP = endP - ( m_KernelNumberOfPoints - 1 ) * m_KernelPointStep;
    }
  }

  unsigned int count = 0;
  for( int p = startP; p <= endP; p += m_KernelPointStep )
    {
    m_KernelTube->GetPoints()[ count ] = tube->GetPoints()[ p ];
    ++count;
    }

  m_KernelTube->ComputeTangentsAndNormals();
}

template< class TInputImage >
void
RadiusExtractor3<TInputImage>
::RecordOptimaAtTubePoints( unsigned int tubePointNum,
  TubeType * tube )
{
  int tubeSize = tube->GetPoints().size();

  int midNum = static_cast<int>(tubePointNum);

  double r1 = this->GetKernelOptimalRadius();
  double m1 = this->GetKernelOptimalRadiusMedialness();
  double b1 = this->GetKernelOptimalRadiusBranchness();
  if( tube->GetPoints()[midNum].GetRadiusInObjectSpace() > 0 )
  {
    r1 = (r1 + tube->GetPoints()[midNum].GetRadiusInObjectSpace())/2.0;
    m1 = (m1 + tube->GetPoints()[midNum].GetMedialness())/2.0;
    b1 = (b1 + tube->GetPoints()[midNum].GetBranchness())/2.0;
  }

  int startP = midNum
    - ( m_KernelNumberOfPoints / 2 ) * m_KernelPointStep - 1;
  if( startP < 0 )
    {
    startP = 0;
    }
  double r0 = tube->GetPoints()[ startP ].GetRadiusInObjectSpace();
  double m0 = tube->GetPoints()[ startP ].GetMedialness();
  double b0 = tube->GetPoints()[ startP ].GetBranchness();
  if( r0 <= 0 )
    {
    r0 = r1;
    m0 = m1;
    b0 = b1;
    }

  int endP = startP + ( m_KernelNumberOfPoints ) * m_KernelPointStep + 1;
  if( endP > tubeSize-1 )
    {
    endP = tubeSize-1;
    }
  double r2 = tube->GetPoints()[ endP ].GetRadiusInObjectSpace();
  double m2 = tube->GetPoints()[ endP ].GetMedialness();
  double b2 = tube->GetPoints()[ endP ].GetBranchness();
  if( r2 <= 0 )
    {
    r2 = r1;
    m2 = m1;
    b2 = b1;
    }

  double rMin = this->GetRadiusMin();
  double rMax = this->GetRadiusMax();
  if( r0 < rMin || r1 < rMin || r2 < rMin )
    {
    std::cout << "ERROR: Max r exceeded." << r0 << ", " << r1 << ", " << r2
      << std::endl;
    }
  if( r0 > rMax || r1 > rMax || r2 > rMax )
    {
    std::cout << "ERROR: Max r exceeded." << r0 << ", " << r1 << ", " << r2
      << std::endl;
    }
  for( int p = startP; p <= endP; ++p )
    {
    if( p < midNum )
      {
      double d = 0;
      if( midNum != startP )
        {
        d = static_cast< double >( midNum - p ) / ( midNum - startP );
        if( d < 0 )
          {
          d = 0;
          }
        if( d > 1 )
          {
          d = 1;
          }
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
        if( d < 0 )
          {
          d = 0;
          }
        if( d > 1 )
          {
          d = 1;
          }
        }
      tube->GetPoints()[ p ].SetRadiusInObjectSpace( d * r2 + ( 1 - d ) * r1 );
      tube->GetPoints()[ p ].SetMedialness( d * m2 + ( 1 - d ) * m1 );
      tube->GetPoints()[ p ].SetBranchness( d * b2 + ( 1 - d ) * b1 );
      }
    if( tube->GetPoints()[p].GetRadiusInObjectSpace() > rMax )
      {
      std::cout << "ERROR: Max r exceeded."
        << tube->GetPoints()[p].GetRadiusInObjectSpace() << std::endl;
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

  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;

  os << indent << "KernelNumberOfPoints = " << m_KernelNumberOfPoints << std::endl;
  os << indent << "KernelPointStep = " << m_KernelPointStep << std::endl;
  os << indent << "KernelStep = " << m_KernelStep << std::endl;

  os << indent << "ProfileNumberOfBins = " << m_ProfileNumberOfBins
    << std::endl;
  os << indent << "ProfileBinValue = " << m_ProfileBinValue.size() << std::endl;
  os << indent << "ProfileBinCount = " << m_ProfileBinCount.size()
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
