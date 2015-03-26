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

#ifndef __itktubeRadiusExtractor_hxx
#define __itktubeRadiusExtractor_hxx

#include "itktubeRadiusExtractor.h"

#include "tubeMatrixMath.h"
#include "tubeTubeMath.h"
#include "tubeUserFunction.h"

#include <itkMinimumMaximumImageFilter.h>

namespace itk
{

namespace tube
{

/** Define the Medialness Function
 * \class RadiusExtractorMedialnessFunc
 */
template< class TImage >
class RadiusExtractorMedialnessFunc : public ::tube::UserFunction< int, double >
{
public:

  typedef itk::VesselTubeSpatialObject< TImage::ImageDimension >
                                                           TubeType;
  typedef typename TubeType::TubePointType                 TubePointType;

  RadiusExtractorMedialnessFunc( RadiusExtractor< TImage > *
    newRadiusExtractor,
    double newMedialnessRadiusStep )
    : m_KernelArray(0),
      m_Value(0.0)
    {
    m_RadiusExtractor = newRadiusExtractor;
    m_MedialnessRadiusStep = newMedialnessRadiusStep;
    }

  void SetKernelArray( std::vector<TubePointType> * newKernelArray )
    {
    m_KernelArray = newKernelArray;
    }

  void SetMedialnessRadiusStep( double newMedialnessRadiusStep )
    {
    m_MedialnessRadiusStep = newMedialnessRadiusStep;
    }

  const double & Value( const int & x )
    {
    double bness = 0;

    m_RadiusExtractor->MeasuresInKernelArray(
      *m_KernelArray, x*m_MedialnessRadiusStep, m_Value, bness, false );

    return m_Value;
    }

  std::vector< TubePointType > * GetKernelArray( void  )
    {
    return m_KernelArray;
    }

  RadiusExtractor< TImage > * GetRadiusExtractor( void  )
    {
    return m_RadiusExtractor;
    }

private:

  std::vector<TubePointType>             * m_KernelArray;

  RadiusExtractor< TImage >              * m_RadiusExtractor;

  double                                   m_Value;

  double                                   m_MedialnessRadiusStep;

}; // End class RadiusExtractorMedialnessFunc

/** Constructor */
template< class TInputImage >
RadiusExtractor<TInputImage>
::RadiusExtractor( void )
{
  m_Image = NULL;
  m_ImageXMin.Fill( 0 );
  m_ImageXMax.Fill( -1 );

  m_DataOp = BlurImageFunction<ImageType>::New();
  m_DataOp->SetScale( 1.0 );
  m_DataOp->SetExtent( 1.1 );
  m_DataMin = 0;
  m_DataMax = -1;

  m_NumKernelPoints = 5;
  m_KernelPointSpacing = 10;

  m_RadiusStart = 1.5;
  m_RadiusMin = 0.5;
  m_RadiusMax = 10.0;

  m_MinMedialness = 0.0015;       // 0.015; larger = harder
  m_MinMedialnessStart = 0.001;

  if( ImageDimension == 2 )
    {
    m_KernNumDirs = 2;
    m_KernX.set_size( ImageDimension, m_KernNumDirs );
    m_KernX(0, 0) = 1;
    m_KernX(1, 0) = 0;
    m_KernX(0, 1) = -1;
    m_KernX(1, 1) = 0;
    }
  else if( ImageDimension == 3 )
    {
    m_KernNumDirs = 8;
    m_KernX.set_size( ImageDimension, m_KernNumDirs );
    int dir = 0;
    for( double theta=0; theta<vnl_math::pi-vnl_math::pi/8;
      theta+=( double )( vnl_math::pi/4 ) )
      {
      m_KernX(0, dir) = vcl_cos( theta );
      m_KernX(1, dir) = vcl_sin( theta );
      m_KernX(2, dir) = 0;
      ++dir;
      m_KernX(0, dir) = -vcl_cos( theta );
      m_KernX(1, dir) = -vcl_sin( theta );
      m_KernX(2, dir) = 0;
      ++dir;
      }
    }
  else
    {
    std::cerr
      << "Error: Radius estimation only supports 2 & 3 dimensions."
      << std::endl;
    throw "Error: Radius estimation only supports 2 & 3 dimensions.";
    }

  m_MedialnessRadiusStep = 0.25;
  m_MedialnessRadiusTolerance = 0.01;

  m_MedialnessFunc = new RadiusExtractorMedialnessFunc<TInputImage>(
    this, m_MedialnessRadiusStep );

  m_MedialnessOpt.SetTolerance( m_MedialnessRadiusTolerance /
    m_MedialnessRadiusStep );

  m_MedialnessOpt.SetXStep( 0.5 / m_MedialnessRadiusStep );

  m_MedialnessOpt.SetSearchForMin( false );

  m_MedialnessOptSpline = new SplineType( m_MedialnessFunc,
    &m_MedialnessOpt );

  m_MedialnessOptSpline->SetClip( true );
  m_MedialnessOptSpline->SetXMin( m_RadiusMin / m_MedialnessRadiusStep );
  m_MedialnessOptSpline->SetXMax( m_RadiusMax / m_MedialnessRadiusStep );

  m_OptimalRadiusMetric = MAXIMIZE_NORMALIZED_DIFFERENCE;

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
}

/** Destructor */
template< class TInputImage >
RadiusExtractor<TInputImage>
::~RadiusExtractor( void )
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


/** Set Radius Min */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetRadiusStart( double radius0 )
{
  m_RadiusStart = radius0;

  if( m_RadiusStart < m_RadiusMin )
    {
    m_RadiusStart = m_RadiusMin;
    }

  if( m_RadiusStart > m_RadiusMax )
    {
    m_RadiusStart = m_RadiusMax;
    }
}

template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetRadiusStep( double radiusStep )
{
  m_MedialnessRadiusStep = radiusStep;
}

template< class TInputImage >
double
RadiusExtractor<TInputImage>
::GetRadiusStep( void )
{
  return m_MedialnessRadiusStep;
}

template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetRadiusTolerance( double radiusTolerance )
{
  m_MedialnessRadiusTolerance = radiusTolerance;
}

template< class TInputImage >
double
RadiusExtractor<TInputImage>
::GetRadiusTolerance( void )
{
  return m_MedialnessRadiusTolerance;
}

/** Set Radius Min */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetRadiusMin( double radiusMin )
{
  if( radiusMin < 0.3 )
    {
    radiusMin = 0.3;
    }
  m_RadiusMin = radiusMin;
  m_MedialnessOptSpline->SetXMin( m_RadiusMin / m_MedialnessRadiusStep );

  if( m_RadiusStart < m_RadiusMin )
    {
    m_RadiusStart = m_RadiusMin;
    }
}

/** Set Radius Max */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetRadiusMax( double radiusMax )
{
  this->m_RadiusMax = radiusMax;
  m_MedialnessOptSpline->SetXMax( this->m_RadiusMax / m_MedialnessRadiusStep );
  if( m_RadiusStart > m_RadiusMax )
    {
    m_RadiusStart = m_RadiusMax;
    }
}

/** Get the medialness operator */
template< class TInputImage >
typename RadiusExtractor<TInputImage>::OptimizerType &
RadiusExtractor<TInputImage>
::GetMedialnessOptimizer( void )
{
  return & m_MedialnessOpt;
}

/** Get the medialness operator */
template< class TInputImage >
typename RadiusExtractor<TInputImage>::SplineType &
RadiusExtractor<TInputImage>
::GetMedialnessOptimizerSpline( void )
{
  return m_MedialnessOptSpline;
}

/** Set the input image */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetInputImage( typename ImageType::Pointer inputImage  )
{
  m_Image = inputImage;

  if( m_Image )
    {
    m_DataOp->SetUseRelativeSpacing( true );
    m_DataOp->SetInputImage( m_Image );
    m_DataOp->SetScale( 3.0 );
    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter =
      MinMaxFilterType::New();
    minMaxFilter->SetInput( m_Image );
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();

    m_ImageXMin = m_Image->GetLargestPossibleRegion().GetIndex();
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      m_ImageXMax[i] = m_ImageXMin[i] +
        m_Image->GetLargestPossibleRegion().GetSize()[i] - 1;
      }

    if( this->GetDebug() )
      {
      std::cout << "RadiusExtractor: SetInputImage: Minimum = "
        << m_DataMin << std::endl;
      std::cout << "RadiusExtractor: SetInputImage: Maximum = "
        << m_DataMax << std::endl;
      }
    }
}

/**
 * Compute the medialness at a point
 */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::MeasuresAtPoint( TubePointType & pnt, double pntR,
  double & mness,
  double & bness,
  bool doBNess )
{
  if( pntR < m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep;
    }

  if( pntR > m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep;
    }

  MatrixType n( ImageDimension, ImageDimension-1 );

  // Verify the normal directions stored in the point are actually
  //   normal to the point's tangent direction.
  double dotP = vnl_math_abs( dot_product( pnt.GetNormal1().GetVnlVector(),
    pnt.GetTangent().GetVnlVector() ) );
  if( ImageDimension == 3 )
    {
    dotP += vnl_math_abs( dot_product( pnt.GetNormal2().GetVnlVector(),
      pnt.GetTangent().GetVnlVector() ) );
    }
  double len = pnt.GetNormal1().GetNorm();
  if( dotP < 0.001 && vnl_math_abs( 1 - len ) < 0.01 )
    {
    n.set_column( 0, pnt.GetNormal1().GetVnlVector() );
    if( ImageDimension == 3 )
      {
      n.set_column( 1, pnt.GetNormal2().GetVnlVector() );
      }
    }
  else
    {
    // If the point's normals, aren't normal, then create new normals.
    //   However, given only one point, we cannot compute the tube's
    //   local Frenet frame.   Ideally the user should call
    //   ComputeTangents on the tube prior to calling this function
    //   to avoid this situation.  If inconsistent normals are used,
    //   branchness computations suffer due to normal flipping.
    if( this->GetDebug() )
      {
      std::cout
        << "Warning: Point normals invalid. Recomputing. Frenet frame lost."
        << std::endl;
      std::cout << "   DotProd = " << dotP << " and Norm = " << len
        << std::endl;
      std::cout << "   pos = " << pnt.GetPosition() << std::endl;
      std::cout << "   t = " << pnt.GetTangent() << std::endl;
      std::cout << "   n1 = " << pnt.GetNormal1() << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   n2 = " << pnt.GetNormal2() << std::endl;
        }
      }
    n.set_column( 0,
      ::tube::ComputeOrthogonalVector( pnt.GetTangent().GetVnlVector() ) );
    n.set_column( 0, n.get_column( 0 ).normalize() );
    if( ImageDimension == 3 )
      {
      n.set_column( 1,
        ::tube::ComputeCrossVector( pnt.GetTangent().GetVnlVector(),
          n.get_column( 0 ) ) );
      n.set_column( 1, n.get_column( 1 ).normalize() );
      pnt.SetNormal1( n.get_column(0)[0], n.get_column(0)[1], n.get_column(0)[2] );
      pnt.SetNormal2( n.get_column(1)[0], n.get_column(1)[1], n.get_column(1)[2] );
      }
    else
      {
      pnt.SetNormal1( n.get_column(0)[0], n.get_column(0)[1] );
      }
    if( this->GetDebug() )
      {
      std::cout << "   new n1 = " << pnt.GetNormal1() << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   new n2 = " << pnt.GetNormal2() << std::endl;
        }
      }
    }

  VectorType kernPos( m_KernNumDirs );
  VectorType kernNeg( m_KernNumDirs );
  VectorType kernBrn( m_KernNumDirs );
  kernPos.fill( 0 );
  kernNeg.fill( 0 );
  kernBrn.fill( 0 );

  if( this->GetDebug() )
    {
    std::cout << "Compute values at point" << std::endl;
    }
  this->ValuesInKernel( pnt, pntR, n, kernPos, kernNeg, kernBrn,
    doBNess );

  this->MeasuresInKernel( pntR, kernPos,
    kernNeg, kernBrn, mness, bness, doBNess );
  pnt.SetMedialness( mness );
  std::cout << "   " << pntR << " : " << mness << std::endl;
  if( doBNess )
    {
    pnt.SetBranchness( bness );
    }
}

/** Compute the medialness at a kernel */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::MeasuresInKernelArray( KernArrayType & kernArray,
  double pntR,
  double & mness,
  double & bness,
  bool doBNess )
{
  unsigned int len = kernArray.size();

  if( len == 0 )
    {
    mness = 0;
    bness = 0;
    return;
    }
  else if( len <= 2 )
    {
    this->MeasuresAtPoint( kernArray[0], pntR, mness, bness,
      doBNess );
    return;
    }

  unsigned int mid = ( len - 1 ) / 2;

  double wTot = 0;
  VectorType w( len );
  for( unsigned int i=0; i<len; i++ )
    {
    w[i] = 1.0 - vnl_math_abs( (double)i - (double)mid ) / ( 2.0 * mid );
    wTot += w[i];
    }
  for( unsigned int i=0; i<len; i++ )
    {
    w[i] /= wTot;
    }

  if( pntR < m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep;
    }

  if( pntR > m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep;
    }

  if( this->GetDebug() )
    {
    std::cout << "RadiusExtractor: MedialnessInKern: pntR = "
      << pntR << " : size = " << kernArray.size()
      << std::endl;
    }

  VectorType kernPosTot( m_KernNumDirs );
  VectorType kernNegTot( m_KernNumDirs );
  VectorType kernBrnTot( m_KernNumDirs );
  kernPosTot.fill( 0 );
  kernNegTot.fill( 0 );
  kernBrnTot.fill( 0 );

  VectorType kernPos( m_KernNumDirs );
  VectorType kernNeg( m_KernNumDirs );
  VectorType kernBrn( m_KernNumDirs );
  kernPos.fill( 0 );
  kernNeg.fill( 0 );
  kernBrn.fill( 0 );

  MatrixType norms( ImageDimension, ImageDimension-1 );
  double dotP = vnl_math_abs( dot_product(
    kernArray[mid].GetNormal1().GetVnlVector(),
    kernArray[mid].GetTangent().GetVnlVector() ) );
  if( ImageDimension == 3 )
    {
    dotP += vnl_math_abs( dot_product(
      kernArray[mid].GetNormal2().GetVnlVector(),
      kernArray[mid].GetTangent().GetVnlVector() ) );
    }
  double sum = kernArray[mid].GetNormal1().GetNorm();
  if( dotP < 0.001 && vnl_math_abs( 1 - sum ) < 0.01 )
    {
    norms.set_column(0, kernArray[mid].GetNormal1().GetVnlVector() );
    if( ImageDimension > 2 )
      {
      norms.set_column(1, kernArray[mid].GetNormal2().GetVnlVector() );
      }
    }
  else
    {
    // If the point's normals, aren't normal, then create new normals.
    //   However, given only one point, we cannot compute the tube's
    //   local Frenet frame.   Ideally the user should call
    //   ComputeTangents on the tube prior to calling this function
    //   to avoid this situation.  If inconsistent normals are used,
    //   branchness computations suffer due to normal flipping.
    //if( this->GetDebug() )
      {
      std::cout
        << "Warning: Point normals invalid. Recomputing. Frenet frame lost."
        << std::endl;
      std::cout << "   DotProd = " << dotP << " and Norm = " << sum
        << std::endl;
      std::cout << "   pos = " << kernArray[mid].GetPosition()
        << std::endl;
      std::cout << "   t = " << kernArray[mid].GetTangent() << std::endl;
      std::cout << "   n1 = " << kernArray[mid].GetNormal1() << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   n2 = " << kernArray[mid].GetNormal2()
          << std::endl;
        }
      }
    norms.set_column( 0,
      ::tube::ComputeOrthogonalVector(
      kernArray[mid].GetTangent().GetVnlVector() ) );
    norms.set_column( 0, norms.get_column( 0 ).normalize() );
    if( ImageDimension == 3 )
      {
      norms.set_column( 1,
        ::tube::ComputeCrossVector(
        kernArray[mid].GetTangent().GetVnlVector(),
        norms.get_column( 0 ) ) );
      norms.set_column( 1, norms.get_column( 1 ).normalize() );
      kernArray[mid].SetNormal1( norms.get_column(0)[0],
        norms.get_column(0)[1], norms.get_column(0)[2] );
      kernArray[mid].SetNormal2( norms.get_column(1)[0],
        norms.get_column(1)[1], norms.get_column(1)[2] );
      }
    else
      {
      kernArray[mid].SetNormal1( norms.get_column(0)[0],
        norms.get_column(0)[1] );
      }
    if( this->GetDebug() )
      {
      std::cout << "   new n1 = " << kernArray[mid].GetNormal1()
        << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   new n2 = " << kernArray[mid].GetNormal2()
          << std::endl;
        }
      }
    }

  // With the coordinate frame defined, compute medialness for other points
  typename std::vector<TubePointType>::iterator pnt = kernArray.begin();
  for( unsigned int i=0; i<len; ++i )
    {
    /*if( this->GetDebug() )
      {
      std::cout << "***___ Kern point = " << i << " ___***" << std::endl;
      }*/
    this->ValuesInKernel( *pnt, pntR, norms, kernPos, kernNeg,
      kernBrn, doBNess );
    /*if( this->GetDebug() )
      {
      std::cout << "  w[" << i << "] = " << w[i] << std::endl;
      }*/
    for( unsigned int d=0; d<m_KernNumDirs; d++ )
      {
      if( this->GetDebug() )
        {
        std::cout << "  Pos[" << d << "] = " << kernPos[d] << std::endl;
        std::cout << "  Neg[" << d << "] = " << kernNeg[d] << std::endl;
        std::cout << "  Brn[" << d << "] = " << kernBrn[d] << std::endl;
        }
      kernPosTot[d] += w[i] * kernPos[d];
      kernNegTot[d] += w[i] * kernNeg[d];
      kernBrnTot[d] += w[i] * kernBrn[d];
      }
    ++pnt;
    }

  mness = 0;
  bness = 0;
  this->MeasuresInKernel( pntR, kernPosTot, kernNegTot, kernBrnTot,
    mness, bness, doBNess );
  if( this->GetDebug() )
    {
    std::cout << "      At radius = " << pntR << " medialness = " << mness
      << std::endl;
    }
}

/** Compute the Optimal scale */
template< class TInputImage >
bool
RadiusExtractor<TInputImage>
::OptimalRadiusAtPoint( TubePointType & pnt,
  double & r0,
  double rMin,
  double rMax,
  double rStep,
  double rTolerance )
{
  TubePointType tmpPnt;
  ITKPointType x;

  KernArrayType kernArray;

  if( rMin < 0.5 )
    {
    rMin = 0.5;
    }

  if( r0 < rMin )
    {
    r0 = rMin;
    }
  else if( r0 > rMax )
    {
    r0 = rMax;
    }

  double xStep = 0.5 * r0;
  if( xStep < 0.5 )
    {
    xStep = 0.5;
    }
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    x[i] = pnt.GetPosition()[i] - pnt.GetTangent()[i] * xStep;
    }
  tmpPnt.SetPosition( x );
  tmpPnt.SetRadius( r0 );
  tmpPnt.SetTangent( pnt.GetTangent() );
  tmpPnt.SetNormal1( pnt.GetNormal1() );
  tmpPnt.SetNormal2( pnt.GetNormal2() );
  kernArray.push_back( tmpPnt );

  ITKPointType x1 = pnt.GetPosition();
  tmpPnt.SetPosition( x1 );
  tmpPnt.SetRadius( r0 );
  tmpPnt.SetTangent( pnt.GetTangent() );
  tmpPnt.SetNormal1( pnt.GetNormal1() );
  tmpPnt.SetNormal2( pnt.GetNormal2() );
  kernArray.push_back( tmpPnt );

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    x1[i] = pnt.GetPosition()[i] + pnt.GetTangent()[i] * xStep;
    }
  tmpPnt.SetPosition( x1 );
  tmpPnt.SetRadius( r0 );
  tmpPnt.SetTangent( pnt.GetTangent() );
  tmpPnt.SetNormal1( pnt.GetNormal1() );
  tmpPnt.SetNormal2( pnt.GetNormal2() );
  kernArray.push_back( tmpPnt );

  double mness = 0.0;

  double tempXMin = m_MedialnessOptSpline->GetXMin();
  m_MedialnessOptSpline->SetXMin( rMin / m_MedialnessRadiusStep );

  double tempXMax = m_MedialnessOptSpline->GetXMax();
  m_MedialnessOptSpline->SetXMax( rMax / m_MedialnessRadiusStep );

  double tempXStep = m_MedialnessOpt.GetXStep();
  m_MedialnessOpt.SetXStep( rStep / m_MedialnessRadiusStep );

  double tempTol = m_MedialnessOpt.GetTolerance();
  m_MedialnessOpt.SetTolerance( rTolerance / m_MedialnessRadiusStep );

  /*
  if( this->GetDebug() )
    {
    std::cout << "kern pnt = " << kernArray.begin()->GetPosition()
    << std::endl;
    } */
  static_cast< RadiusExtractorMedialnessFunc< TInputImage > *>(
    m_MedialnessFunc )->SetKernelArray( & kernArray );
  m_MedialnessOptSpline->SetNewData( true );
  double oldR0 = r0;
  r0 /= m_MedialnessRadiusStep;
  m_MedialnessOptSpline->Extreme( &r0, &mness );
  //if( this->GetDebug() )
    {
    std::cout << " cmp: " << r0-0.1/m_MedialnessRadiusStep << " - "
      << m_MedialnessOptSpline->Value( r0-0.1/m_MedialnessRadiusStep )
      << std::endl;
    std::cout << " cmp: " << r0 << " - "
      << m_MedialnessOptSpline->Value( r0 ) << std::endl;
    std::cout << " cmp: " << r0+0.1/m_MedialnessRadiusStep << " - "
      << m_MedialnessOptSpline->Value( r0+0.1/m_MedialnessRadiusStep )
      << std::endl;
    }
  r0 *= m_MedialnessRadiusStep;

  if( this->GetDebug() )
    {
    std::cout << "Local extreme at radius r0 = " << r0
      << " with medialness = " << mness << std::endl;
    std::cout << "  prev radius = " << oldR0 << std::endl;
    std::cout << std::endl;
    }

  m_MedialnessOptSpline->SetXMin( tempXMin );
  m_MedialnessOptSpline->SetXMax( tempXMax );
  m_MedialnessOpt.SetXStep( tempXStep );
  m_MedialnessOpt.SetTolerance( tempTol );

  pnt.SetRadius( r0 );

  if( mness > m_MinMedialnessStart )
    {
    return true;
    }
  else
    {
    if ( this->GetDebug() )
      {
      std::cout
        << "RadiusExtractor: calcOptimalScale: kernel fit insufficient"
        << std::endl;
      std::cout << "  Medialness = " << mness << " < thresh = "
        << m_MinMedialness << std::endl;
      }
    return false;
    }

  return true;
}

template< class TInputImage >
void
RadiusExtractor<TInputImage>
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
  os << indent << "ImageXMin = " << m_ImageXMin << std::endl;
  os << indent << "ImageXMax = " << m_ImageXMax << std::endl;
  if( m_DataOp.IsNotNull() )
    {
    os << indent << "DataOp = " << m_DataOp << std::endl;
    }
  else
    {
    os << indent << "DataOp = NULL" << std::endl;
    }
  os << indent << "MedialnessRadiusStep = " << m_MedialnessRadiusStep
    << std::endl;
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
  os << indent << "NumKernelPoints = " << m_NumKernelPoints << std::endl;
  os << indent << "KernelPointSpacing = " << m_KernelPointSpacing
    << std::endl;
  os << indent << "DataMin = " << m_DataMin << std::endl;
  os << indent << "DataMax = " << m_DataMax << std::endl;
  os << indent << "RadiusStart = " << m_RadiusStart << std::endl;
  os << indent << "RadiusMin = " << m_RadiusMin << std::endl;
  os << indent << "RadiusMax = " << m_RadiusMax << std::endl;
  os << indent << "MinMedialness = " << m_MinMedialness << std::endl;
  os << indent << "MinMedialnessStart = " << m_MinMedialnessStart
    << std::endl;
  os << indent << "MedialnessFunc = " << m_MedialnessFunc << std::endl;
  os << indent << "KernNumDirs = " << m_KernNumDirs << std::endl;
  os << indent << "KernX = " << m_KernX << std::endl;
  os << indent << "IdleCallBack = " << m_IdleCallBack << std::endl;
  os << indent << "StatusCallBack = " << m_StatusCallBack << std::endl;
}

template< class TInputImage >
void
RadiusExtractor<TInputImage>
::ValuesInSubKernel( TubePointType pnt, double pntR,
  MatrixType & kernN, VectorType & kern, double & kernCnt )
{
  if( pntR < m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep;
    }

  if( pntR > m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep;
    }

  /*
  if( this->GetDebug() )
    {
    std::cout << "Kern pnt = " << pnt.GetPosition() << std::endl;
    for( unsigned int i=0; i<ImageDimension-1; i++ )
      {
      std::cout << "Kern N[" << i << "] = " << kernN.get_column( i )
        << std::endl;
      }
    std::cout << "Kern r = " << pntR << std::endl;
    std::cout << "Kern X = " << m_KernX << std::endl;
    }
  */

  VectorType nodePnt;
  for( unsigned int dir=0; dir<m_KernNumDirs; dir++ )
    {
    nodePnt = pnt.GetPosition().GetVnlVector();
    for( unsigned int i=0; i<ImageDimension-1; i++ )
      {
      for( unsigned int d=0; d<ImageDimension; d++ )
        {
        nodePnt[d] += pntR * m_KernX( i, dir ) * kernN( d, i );
        }
      }

    if( this->GetDebug() )
      {
      VectorType xV = nodePnt - pnt.GetPosition().GetVnlVector();
      std::cout << "  dirX = " << m_KernX.get_column( dir ) << std::endl;
      std::cout << "  node pnt = " << nodePnt << std::endl;
      std::cout << "  dist(kernPnt, nodePnt) = " << xV.magnitude()
        << std::endl;
      xV = xV.normalize();
      double dotP = dot_product( pnt.GetTangent().GetVnlVector(), xV );
      std::cout << "  dotp(kernPnt.tan, nodePnt.dist) = " << dotP
        << std::endl;
      }

    bool inBounds = true;
    ContinuousIndex<double, ImageDimension> nodeCIndx;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      if( nodePnt[i] < m_ImageXMin[i] || nodePnt[i] > m_ImageXMax[i] )
        {
        inBounds = false;
        break;
        }
      nodeCIndx[i] = nodePnt[i];
      }

    if( inBounds )
      {
      double val = ( m_DataOp->EvaluateAtContinuousIndex( nodeCIndx )
        - m_DataMin ) / ( m_DataMax - m_DataMin );

      if( val < 0 )
        {
        val = 0;
        }

      if( val > 1 )
        {
        val = 1;
        }

      kern[dir] = val;
      ++kernCnt;
      }
    else
      {
      kern[dir] = 0;
      }
    }
}

template< class TInputImage >
void
RadiusExtractor<TInputImage>
::ValuesInKernel( TubePointType pnt, double pntR,
  MatrixType & kernN, VectorType & kernPos, VectorType & kernNeg,
  VectorType & kernBrn, bool doBNess )
{
  if( pntR < m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMin() * m_MedialnessRadiusStep;
    }

  if( pntR > m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep )
    {
    pntR = m_MedialnessOptSpline->GetXMax() * m_MedialnessRadiusStep;
    }

  MatrixType n( ImageDimension, ImageDimension-1 );

  // Verify the normal directions stored in the point are actually
  //   normal to the point's tangent direction.
  double dotP = vnl_math_abs( dot_product( pnt.GetNormal1().GetVnlVector(),
    pnt.GetTangent().GetVnlVector() ) );
  if( ImageDimension == 3 )
    {
    dotP += vnl_math_abs( dot_product( pnt.GetNormal2().GetVnlVector(),
      pnt.GetTangent().GetVnlVector() ) );
    }
  double sum = pnt.GetNormal1().GetNorm();
  if( dotP < 0.001 && vnl_math_abs( 1 - sum ) < 0.01 )
    {
    n.set_column( 0, pnt.GetNormal1().GetVnlVector() );
    if( ImageDimension == 3 )
      {
      n.set_column( 1, pnt.GetNormal2().GetVnlVector() );
      }
    }
  else
    {
    // If the point's normals, aren't normal, then create new normals.
    //   However, given only one point, we cannot compute the tube's
    //   local Frenet frame.   Ideally the user should call
    //   ComputeTangents on the tube prior to calling this function
    //   to avoid this situation.  If inconsistent normals are used,
    //   branchness computations suffer due to normal flipping.
    if( this->GetDebug() )
      {
      std::cout
        << "Warning: Point normals invalid. Recomputing. Frenet frame lost."
        << std::endl;
      std::cout << "   DotProd = " << dotP << " and Norm = " << sum
        << std::endl;
      std::cout << "   pos = " << pnt.GetPosition() << std::endl;
      std::cout << "   t = " << pnt.GetTangent() << std::endl;
      std::cout << "   n1 = " << pnt.GetNormal1() << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   n2 = " << pnt.GetNormal2() << std::endl;
        }
      }
    n.set_column( 0,
      ::tube::ComputeOrthogonalVector( pnt.GetTangent().GetVnlVector() ) );
    n.set_column( 0, n.get_column( 0 ).normalize() );
    if( ImageDimension == 3 )
      {
      n.set_column( 1,
        ::tube::ComputeCrossVector( pnt.GetTangent().GetVnlVector(),
        n.get_column( 0 ) ) );
      n.set_column( 1, n.get_column( 1 ).normalize() );
      pnt.SetNormal1( n.get_column(0)[0], n.get_column(0)[1],
        n.get_column(0)[2] );
      pnt.SetNormal2( n.get_column(1)[0], n.get_column(1)[1],
        n.get_column(1)[2] );
      }
    else
      {
      pnt.SetNormal1( n.get_column(0)[0], n.get_column(0)[1] );
      }
    if( this->GetDebug() )
      {
      std::cout << "   new n1 = " << pnt.GetNormal1() << std::endl;
      if( ImageDimension == 3 )
        {
        std::cout << "   new n2 = " << pnt.GetNormal2() << std::endl;
        }
      }
    }

  VectorType n0( ImageDimension );
  n0.fill( 0 );
  for( unsigned int j=0; j<ImageDimension-1; j++ )
    {
    double dp = dot_product( kernN.get_column( 0 ),
      n.get_column( j ) );
    n0 += dp * n.get_column( j );
    }
  n0 = n0.normalize();
  if( this->GetDebug() )
    {
    std::cout << "n0 = " << n0 << std::endl;
    std::cout << "kernN0 = " << kernN.get_column(0) << std::endl;
    }
  n.set_column( 0, n0 );
  if( ImageDimension == 3 )
    {
    n.set_column( 1,
      ::tube::ComputeCrossVector( pnt.GetTangent().GetVnlVector(),
      n.get_column( 0 ) ) );
    n.set_column( 1, n.get_column( 1 ).normalize() );
    if( this->GetDebug() )
      {
      std::cout << "n1 = " << n.get_column( 1 ) << std::endl;
      std::cout << "kernN1 = " << kernN.get_column(1) << std::endl;
      }
    double tf = dot_product( kernN.get_column( 1 ),
      n.get_column( 1 ) );
    if( tf < 0 )
      {
      n.set_column( 1, n.get_column( 1 ) * -1 );
      }
    }

  kernPos.fill( 0 );
  kernNeg.fill( 0 );
  kernBrn.fill( 0 );
  double kernPosCnt = 0;
  double kernNegCnt = 0;
  double kernBrnCnt = 0;

  double e = 1.0;
  double f = 2.0;
  if( ( pntR / f ) * e < 0.71 )  // std::sqrt( 2 )/2 = 0.7071 approx = 0.71
    {
    e = ( 0.71 * f ) / pntR;
    }
  double r = pntR - (pntR/f) * e;
  if( r < 0 )
    {
    r = 0;
    e = 1;
    f = 1;
    }
  m_DataOp->SetScale( pntR / f  );
  m_DataOp->SetExtent( e );
  if( this->GetDebug() )
    {
    std::cout << "Pos: opR = " << pntR << " r = " << r << " opE = " << e
      << " dist = " << r + pntR/f * e << std::endl;
    std::cout << "   n = " << n << std::endl;
    }
  this->ValuesInSubKernel( pnt, r, n, kernPos, kernPosCnt );
  if( this->GetDebug() )
    {
    std::cout << "   kernPos = " << kernPos << std::endl;
    }

  r = pntR + (pntR - r);
  if( this->GetDebug() )
    {
    std::cout << "Neg: opR = " << pntR << " r = " << r << " opE = " << e
      << " dist = " << r - pntR/f * e << std::endl;
    }
  this->ValuesInSubKernel( pnt, r, n, kernNeg, kernNegCnt );
  if( this->GetDebug() )
    {
    std::cout << "   kernNeg = " << kernNeg << std::endl;
    }

  if( doBNess )
    {
    e = 1.1;
    f = 3.0;
    if( ( pntR / f ) * e < 1.1 )  // std::sqrt( 2 )/2 = 0.7071 approx = 0.71
      {
      f = ( ( pntR * e ) / 1.1 + f ) / 2;
      e = 1.1 / ( pntR / f );
      }
    if( ( pntR / f ) * e > 3.1 )
      {
      f = ( ( pntR * e ) / 3.1 + f ) / 2;
      e = 3.1 / ( pntR / f );
      }
    m_DataOp->SetScale( pntR / f );
    m_DataOp->SetExtent( e );
    r = f * pntR;
    if( this->GetDebug() )
      {
      std::cout << "Brn: opR = " << pntR/f << " opE = " << e
        << " dist = " << r << std::endl;
      }
    this->ValuesInSubKernel( pnt, r, n, kernBrn, kernBrnCnt );
    }

  int kernCnt = 0;
  for( unsigned int i=0; i<m_KernNumDirs; i++ )
    {
    if( kernNeg[i]!=0 && kernPos[i]!=0 )
      {
      kernCnt++;
      }
    }
  if( kernCnt < 2 )
    {
    if( this->GetDebug() )
      {
      std::cout << "Warning: Medialness kernel does not intersect image."
        << std::endl;
      std::cout << "   Warn pos = " << pnt.GetPosition() << std::endl;
      std::cout << "   Warn tan = " << pnt.GetTangent() << std::endl;
      std::cout << "   Warn kern = " << std::endl;
      for( unsigned int i=0; i<m_KernNumDirs; i++ )
        {
        std::cout << "    " << i << " : " << kernPos[i] << ", "
          << kernNeg[i] << ", " << kernBrn[i] << std::endl;
        }
      }
    }
}

/** Compute kernel array */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::ValuesInFullKernelArray( TubeType * tube,
  KernArrayType & kernArray,
  KernArrayTubePointIndexType & kernArrayTubePointIndex )
{
  typename KernArrayType::iterator tubeFromPnt = tube->GetPoints().begin();
  typename KernArrayType::iterator tubeToPnt = tube->GetPoints().end();

  typename KernArrayType::iterator iterPnt;
  iterPnt = tubeFromPnt;

  kernArray.clear();
  kernArrayTubePointIndex.clear();

  unsigned int kernPointCount = 0;
  unsigned int tubePointCount = 0;
  while( iterPnt != tubeToPnt )
    {
    int avgCount = 0;
    TubePointType kernPnt;
    for( unsigned int j = 0;
      iterPnt != tubeToPnt && j < m_KernelPointSpacing;
      j++ )
      {
      ITKPointType tmpPoint;

      for( unsigned int id = 0; id < ImageDimension; id++ )
        {
        tmpPoint[id] = kernPnt.GetPosition()[id]
          + ( *iterPnt ).GetPosition()[id];
        }
      kernPnt.SetPosition( tmpPoint );

      double tmpRadius = kernPnt.GetRadius() + (*iterPnt).GetRadius();
      kernPnt.SetRadius( tmpRadius );

      double dotP = dot_product( kernPnt.GetTangent().GetVnlVector(),
        ( *iterPnt ).GetTangent().GetVnlVector() );
      if( dotP > 0 )
        {
        kernPnt.SetTangent( kernPnt.GetTangent() +
          ( *iterPnt ).GetTangent() );
        }
      else
        {
        kernPnt.SetTangent( kernPnt.GetTangent() -
          ( *iterPnt ).GetTangent() );
        }

      avgCount++;
      iterPnt++;
      tubePointCount++;
      }

    ITKPointType tmpPoint;
    for( unsigned int id = 0; id < ImageDimension; id++ )
      {
      tmpPoint[id] = kernPnt.GetPosition()[id] / avgCount;
      }
    kernPnt.SetPosition( tmpPoint );

    ITKVectorType tempVect = kernPnt.GetTangent();
    tempVect.Normalize();
    kernPnt.SetTangent( tempVect );

    kernPnt.SetRadius( kernPnt.GetRadius() / avgCount );

    if( iterPnt == tubeToPnt )
      {
      if( kernPointCount == 0 )
        {
        kernArray.push_back( kernPnt );
        kernArrayTubePointIndex.push_back( tubePointCount );
        kernPointCount++;
        }
      return;
      }

    kernArray.push_back( kernPnt );
    kernArrayTubePointIndex.push_back( tubePointCount );
    kernPointCount++;
    }
}


/**
 * Compute the medialness and the branchness */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::MeasuresInKernel( double pntR,
  VectorType & kernPos,
  VectorType & kernNeg,
  VectorType & kernBrn,
  double & mness,
  double & bness,
  bool doBNess )
{
  int kernAvgCnt = 0;
  for( unsigned int i=0; i<m_KernNumDirs; i++ )
    {
    if( kernNeg[i] != 0 && kernPos[i] != 0 )
      {
      kernAvgCnt++;
      }
    }

  if( kernAvgCnt<0.25*m_KernNumDirs )
    {
    if ( this->GetDebug() )
      {
      std::cout
        << "Error: insufficient neg/pos kernel pairs intersect image."
        << std::endl;
      }
    mness = 0;
    bness = 0;
    return;
    }

  int kernNegMaxI = -1;
  for( unsigned int i=0; i<m_KernNumDirs; i++ )
    {
    if( kernNeg[i] != 0 && kernPos[i] != 0 )
      {
      kernNegMaxI = i;
      break;
      }
    }
  double kernNegMax = kernNeg[kernNegMaxI];
  double kernNegAvg = kernNegMax;

  int kernPosMinI = kernNegMaxI;
  double kernPosMin = kernPos[kernPosMinI];
  double kernPosAvg = kernPosMin;

  for( unsigned int i = kernNegMaxI+1; i < m_KernNumDirs; i++ )
    {
    if( kernNeg[i] != 0 && kernPos[i] != 0 )
      {
      if( kernNeg[i]>kernNegMax )
        {
        kernNegMax = kernNeg[i];
        kernNegMaxI = i;
        }
      if( kernPos[i]<kernPosMin )
        {
        kernPosMin = kernPos[i];
        kernPosMinI = i;
        }
      kernNegAvg += kernNeg[i];
      kernPosAvg += kernPos[i];
      }
    }

  if( kernAvgCnt > 5 )
    {
    int kernI = kernNegMaxI;  // first dampen most positive negative-node
    kernPosAvg -= kernPos[kernI];
    kernNegAvg -= kernNeg[kernI];
    /* Dampens the effect of the extreme, by replacing it with a point
     * 1/2 to the extreme. */
    kernPosAvg += ( kernPos[kernI] + kernPosAvg / ( kernAvgCnt-1 ) )
      / 2;
    kernNegAvg += ( kernNeg[kernI] + kernNegAvg / ( kernAvgCnt-1 ) )
      / 2;

    if( kernPosMinI == kernNegMaxI )
      {
      /* Dampens the effect of the largest adjacent node */
      int l = kernI + 1;
      if( l >= (int)m_KernNumDirs )
        {
        l = 0;
        }
      int m = kernI - 1;
      if( m < 0 )
        {
        m = m_KernNumDirs - 1;
        }

      if( kernNeg[l] != 0 && kernPos[l] != 0 &&
          kernNeg[m] != 0 && kernPos[m] != 0 )
        {
        int kernAdjI;
        if( ( kernPos[l] - kernNeg[l] ) < ( kernPos[m] - kernNeg[m] ) )
          {
          kernAdjI = l;
          }
        else
          {
          kernAdjI = m;
          }
        kernPosAvg -= kernPos[kernAdjI];
        kernNegAvg -= kernNeg[kernAdjI];
        kernPosAvg += ( kernPos[kernAdjI] + kernPosAvg / (kernAvgCnt-1) )
          / 2;
        kernNegAvg += ( kernNeg[kernAdjI] + kernNegAvg / (kernAvgCnt-1) )
          / 2;
        }
      }
    else
      {
      kernI = kernPosMinI; // second dampen most negative positive-node
      kernPosAvg -= kernPos[kernI];
      kernNegAvg -= kernNeg[kernI];
      /* Dampens the effect of the extreme, by replacing it with a point
       * 1/2 to the extreme. */
      kernPosAvg += ( kernPos[kernI] + kernPosAvg / ( kernAvgCnt-1 ) )
        / 2;
      kernNegAvg += ( kernNeg[kernI] + kernNegAvg / ( kernAvgCnt-1 ) )
        / 2;
      }
    }

  if( kernAvgCnt == 0 )
    {
    kernPosAvg = 0;
    kernNegAvg = 0;
    mness = 0;
    }
  else
    {
    kernPosAvg /= kernAvgCnt;
    kernNegAvg /= kernAvgCnt;
    switch( m_OptimalRadiusMetric )
      {
      default:
      case MAXIMIZE_DIFFERENCE:
        {
        mness = ( kernPosAvg - kernNegAvg );
        break;
        }
      case MAXIMIZE_RATIO:
        {
        if( kernPosAvg != 0 )
          {
          // We put kernPosAvg in the denominator since kernNegAvg could
          // go towards zero resulting in INF results.   This way the
          // result varies 0 to 1, in most cases;
          mness = kernPosAvg / kernNegAvg;
          }
        else
          {
          mness = 0;
          }
        break;
        }
      case MAXIMIZE_NORMALIZED_DIFFERENCE:
        {
        if( kernPosAvg != 0 )
          {
          // See comment on MAXIMIZE_RATIO to understand this formulation.
          mness = ( kernPosAvg - kernNegAvg ) / kernNegAvg;
          }
        else
          {
          mness = 0;
          }
        break;
        }
      }
    }


  // if( this->GetDebug() )
    {
    std::cout << "  PntR = " << pntR << " MNESS = " << mness
      << " kernPosAvg = " << kernPosAvg
      << " kernNegAvg = " << kernNegAvg << std::endl;
    }

  if( doBNess )
    {
    unsigned int kernBrnMaxI = m_KernNumDirs;
    double kernBrnMax;
    double kernBrnAvg;
    double kernBrnAvgCnt = 1;
    for( unsigned int i=0; i<m_KernNumDirs; i++ )
      {
      if( kernPos[i] != 0 && kernNeg[i] != 0 && kernBrn[i] != 0 )
        {
        kernBrnMaxI = i;
        break;
        }
      }

    if( kernBrnMaxI < m_KernNumDirs )
      {
      kernBrnMax = kernBrn[kernBrnMaxI];
      kernBrnAvg = kernBrnMax;

      for( unsigned int i = kernBrnMaxI+1; i<m_KernNumDirs; i++ )
        {
        if( kernPos[i] != 0 && kernNeg[i] != 0 && kernBrn[i] != 0 )
          {
          kernBrnAvg += kernBrn[i];
          kernBrnAvgCnt++;
          if( kernBrn[i] > kernBrnMax )
            {
            kernBrnMax = kernBrn[i];
            }
          }
        }

      kernBrnAvg -= kernBrnMax;
      if( ( kernBrnAvgCnt-1 ) != 0 )
        {
        kernBrnAvg /= ( kernBrnAvgCnt-1 );
        }
      else
        {
        kernBrnAvg=0;
        }

      if( ( kernPosAvg - kernBrnAvg ) != 0 )
        {
        bness = ( 3 * ( kernBrnMax - kernBrnAvg )
          / ( kernPosAvg - kernBrnAvg )
          + ( kernNegMax - kernNegAvg )
          / ( kernPosAvg - kernBrnAvg ) ) / 4;
        }
      else
        {
        bness=0;
        }

      if( bness>2 )
        {
        bness = 2;
        }
      if( bness<0 )
        {
        bness = 0;
        }
      }
    else
      {
      bness = 0;
      }
    }
  else
    {
    bness = 0;
    }
}


/** Compute Kernel radii one way */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::MeasuresInFullKernelArray( KernArrayType & kernArray,
  unsigned int kernPntStart,
  unsigned int kernPntEnd,
  double radiusStart )
{
  double pntR = radiusStart;
  double prevPntR = pntR;
  double mness = 0;

  unsigned int kernMid = ( m_NumKernelPoints - 1 ) / 2;

  KernArrayType pntKernArray;
  int step = 1;
  if( kernPntStart > kernPntEnd )
    {
    step = -1;
    }

  unsigned int kernArraySize = kernArray.size();
  unsigned int count = 0;
  for( int kernPnt = (int)(kernPntStart);
    kernPnt != (int)(kernPntEnd)+step;
    kernPnt += step )
    {
    pntKernArray.clear();
    for( int j = (int)(kernPnt) - (int)(kernMid);
      j <= (int)(kernPnt) + (int)(kernMid); j++ )
      {
      if( j >= 0 && j < (int)(kernArraySize) )
        {
        pntKernArray.push_back( kernArray[j] );
        }
      }

    static_cast< RadiusExtractorMedialnessFunc< TInputImage > *>(
      m_MedialnessFunc )->SetKernelArray( & pntKernArray );
    m_MedialnessOptSpline->SetNewData( true );
    double oldPntR = pntR;
    pntR /= m_MedialnessRadiusStep;
    m_MedialnessOptSpline->Extreme( &pntR, &mness );
    pntR *= m_MedialnessRadiusStep;

    //if( this->GetDebug() )
      {
      std::cout << "Local extreme at radius pntR = " << pntR
        << " with medialness = " << mness << std::endl;
      std::cout << "  prev radius = " << oldPntR << std::endl;
      /*
      std::cout << "  cmp: " << pntR-0.1 << " - "
        << m_MedialnessOptSpline->Value( (pntR-0.1)/m_MedialnessRadiusStep )
        << std::endl;
      std::cout << "  cmp: " << pntR << " - "
        << m_MedialnessOptSpline->Value( pntR / m_MedialnessRadiusStep )
        << std::endl;
      std::cout << "  cmp: " << pntR+0.1 << " - "
        << m_MedialnessOptSpline->Value( (pntR+0.1)/m_MedialnessRadiusStep )
        << std::endl;
        */
      std::cout << std::endl;
      }

    if( mness < m_MinMedialness )
      {
      //if( this->GetDebug() )
        {
        std::cout << "Bad mnessVal( " << pntR << " ) = " << mness << std::endl;
        }
      pntR = prevPntR;
      oldPntR = pntR;
      pntR /= m_MedialnessRadiusStep;
      m_MedialnessOptSpline->Extreme( &pntR, &mness );
      pntR *= m_MedialnessRadiusStep;
      if( this->GetDebug() )
        {
        std::cout << "NEW local extreme at radius pntR = " << pntR
          << " with medialness = " << mness << std::endl;
        std::cout << "  prev radius = " << oldPntR << std::endl;
        std::cout << "  cmp: " << pntR-0.1 << " - "
          << m_MedialnessOptSpline->Value( (pntR-0.1)/m_MedialnessRadiusStep )
          << std::endl;
        std::cout << "  cmp: " << pntR << " - "
          << m_MedialnessOptSpline->Value( pntR / m_MedialnessRadiusStep )
          << std::endl;
        std::cout << "  cmp: " << pntR+0.1 << " - "
          << m_MedialnessOptSpline->Value( (pntR+0.1)/m_MedialnessRadiusStep )
          << std::endl;
        std::cout << std::endl;
        }

      if( mness >= m_MinMedialness )
        {
        if( this->GetDebug() )
          {
          std::cout << "   *** new mnessVal( " << pntR << " ) = " << mness
            << std::endl;
          }
        if( mness > 2 * m_MinMedialness )
          {
          prevPntR = pntR;
          }
        }
      else
        {
        pntR = prevPntR;
        mness = 2 * m_MinMedialness;
        if( this->GetDebug() )
          {
          std::cout << "   using old mnessVal( " << pntR << " ) = "
            << mness << std::endl;
          }
        }
      }
    else
      {
      if( mness > 2 * m_MinMedialness )
        {
        prevPntR = pntR;
        }
      }

    kernArray[kernPnt].SetRadius( pntR );

    count++;
    if( count/5 == count/5.0 && m_StatusCallBack )
      {
      char mesg[80];
      std::sprintf( mesg, "Radius at %d = %f", count*m_KernelPointSpacing,
        pntR );
      char loc[80];
      std::sprintf( loc, "Extract:Widths" );
      m_StatusCallBack( loc, mesg, 0 );
      }
    }
}

/** Compute Kernel Measures */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SmoothMeasuresInFullKernelArray( KernArrayType & kernArray )
{
  unsigned int len = kernArray.size();

  if( this->GetDebug() )
    {
    std::cout << "Smoothing..." << std::endl;
    for( unsigned int kernPnt=0; kernPnt<len-1; kernPnt++ )
      {
      std::cout << kernPnt << " : r = " << kernArray[kernPnt].GetRadius()
        << std::endl;
      }
    }

  for( unsigned int iter=0; iter<2; iter++ )
    {
    for( unsigned int kernPnt=0; kernPnt<len-1; kernPnt++ )
      {
      kernArray[kernPnt].SetRadius(
        ( 1.0 * kernArray[kernPnt].GetRadius()
        + 1.0 * kernArray[kernPnt+1].GetRadius() ) / 2.0 );
      }
    for( int kernPnt=(int)(len)-1; kernPnt>0; kernPnt-- )
      {
      kernArray[kernPnt].SetRadius(
        ( 1.0 * kernArray[kernPnt].GetRadius()
        + 1.0 * kernArray[kernPnt-1].GetRadius() ) / 2.0 );
      }
    }

  unsigned int kernMid = ( m_NumKernelPoints - 1 ) / 2;

  KernArrayType kernTemp;
  for( unsigned int kernPnt=0; kernPnt<len; kernPnt++ )
    {
    kernTemp.clear();
    for( int j = (int)( kernPnt )-(int)( kernMid );
      j <= (int)( kernPnt )+(int)( kernMid ); j++ )
      {
      if( j >= 0 && j < (int)( len ) )
        {
        kernTemp.push_back( kernArray[j] );
        }
      else
        {
        kernTemp.push_back( kernArray[kernPnt] );
        }
      }

    double mness = 0;
    double bness = 0;
    this->MeasuresInKernelArray( kernTemp,
      kernArray[kernPnt].GetRadius(), mness, bness, true );
    kernArray[kernPnt].SetMedialness( mness );
    kernArray[kernPnt].SetBranchness( bness );
    }
}

/**
 * Apply kernel measures */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::ApplyMeasuresInFullKernelArray( TubeType * tube,
  KernArrayType & kernArray,
  KernArrayTubePointIndexType & kernArrayTubePointIndex )
{

  if( this->GetDebug() )
    {
    std::cout << "Applying measures." << std::endl;
    }

  unsigned int len = kernArray.size();
  if( this->GetDebug() )
    {
    std::cout << "Array len = " << len << std::endl;
    std::cout << "Point = " << kernArrayTubePointIndex[0] << std::endl;
    }

  if( len == 0 )
    {
    return;
    }

  if( kernArrayTubePointIndex[0] >= tube->GetPoints().size() )
    {
    return;
    }

  double r0 = kernArray[0].GetRadius();
  double m0 = kernArray[0].GetMedialness();
  double b0 = kernArray[0].GetBranchness();

  KernArrayType & pnts = tube->GetPoints();
  typename KernArrayType::iterator pntIter = pnts.begin();
  unsigned int pntCnt = 0;

  while( pntCnt < kernArrayTubePointIndex[0] )
    {
    if( this->GetDebug() )
      {
      std::cout << "Lead: " << pntCnt << " : r = " << r0 << std::endl;
      }
    (*pntIter).SetRadius( r0 );
    (*pntIter).SetMedialness( m0 );
    (*pntIter).SetBranchness( b0 );
    ++pntIter;
    ++pntCnt;
    }

  double w;
  double r;
  double m;
  double b;

  double r1 = kernArray[1].GetRadius();
  double m1 = kernArray[1].GetMedialness();
  double b1 = kernArray[1].GetBranchness();
  for( unsigned int arrayCnt=0; arrayCnt<len-1; arrayCnt++ )
    {
    r0 = kernArray[arrayCnt].GetRadius();
    m0 = kernArray[arrayCnt].GetMedialness();
    b0 = kernArray[arrayCnt].GetBranchness();
    r1 = kernArray[arrayCnt+1].GetRadius();
    m1 = kernArray[arrayCnt+1].GetMedialness();
    b1 = kernArray[arrayCnt+1].GetBranchness();
    double kernCnt = 0;
    double kernCntMax = kernArrayTubePointIndex[arrayCnt+1] - pntCnt;
    while( pntCnt < kernArrayTubePointIndex[arrayCnt+1] )
      {
      w = kernCnt / kernCntMax;
      r = r0 + w * ( r1 - r0 );
      m = m0 + w * ( m1 - m0 );
      b = b0 + w * ( b1 - b0 );
      if( this->GetDebug() )
        {
        std::cout << "Mid: " << arrayCnt << " : " << pntCnt
          << " : r = " << r << std::endl;
        }
      (*pntIter).SetRadius( r );
      (*pntIter).SetMedialness( m );
      (*pntIter).SetBranchness( b );
      ++pntIter;
      ++pntCnt;
      ++kernCnt;
      }
    }

  r0 = kernArray[len-1].GetRadius();
  m0 = kernArray[len-1].GetMedialness();
  b0 = kernArray[len-1].GetBranchness();
  while( pntIter != tube->GetPoints().end() )
    {
    (*pntIter).SetRadius( r0 );
    (*pntIter).SetMedialness( m0 );
    (*pntIter).SetBranchness( b0 );
    if( this->GetDebug() )
      {
      std::cout << "End: " << pntCnt << " : r = " << r0 << std::endl;
      }
    ++pntIter;
    ++pntCnt;
    }

  if( m_StatusCallBack )
    {
    char mesg[80];
    std::sprintf( mesg, "Applied to %ld tube points.",
      (long int)(tube->GetPoints().size()) );
    char loc[80];
    std::sprintf( loc, "Extract:Widths" );
    m_StatusCallBack( loc, mesg, 0 );
    }
}

/** Specify function to be optimized when selecting a radius */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetOptimalRadiusMetric( const OptimalRadiusMetricEnum & metric )
{
  m_OptimalRadiusMetric = metric;
}

/** Get function to be optimized when selecting a radius */
template< class TInputImage >
const typename RadiusExtractor<TInputImage>::OptimalRadiusMetricEnum &
RadiusExtractor<TInputImage>
::GetOptimalRadiusMetric( void ) const
{
  return m_OptimalRadiusMetric;
}

/** Compute Radii */
template< class TInputImage >
bool
RadiusExtractor<TInputImage>
::ExtractRadii( TubeType * tube )
{
  if( tube->GetPoints().size() == 0 )
    {
    return false;
    }

  KernArrayType               kernArray;
  KernArrayTubePointIndexType kernArrayTubePointIndex;

  this->ValuesInFullKernelArray( tube, kernArray,
    kernArrayTubePointIndex );

  if( this->GetDebug() )
    {
    std::cout << "Computing tangents and normals..." << std::endl;
    }

  ::tube::ComputeVectorTangentsAndNormals< TubePointType >( kernArray );

  unsigned int len = kernArray.size();

  typename KernArrayType::iterator pntIter;
  pntIter = tube->GetPoints().begin();
  while( pntIter != tube->GetPoints().end() && (*pntIter).GetID() != 0 )
    {
    pntIter++;
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

  int minDistI = 0;
  double minDist = ::tube::ComputeEuclideanDistance(
    kernArray[0].GetPosition(), (*pntIter).GetPosition() );
  for( unsigned int kPnt=1; kPnt<len; kPnt++ )
    {
    double tf = ::tube::ComputeEuclideanDistance(
      kernArray[kPnt].GetPosition(), (*pntIter).GetPosition() );
    if( tf < minDist )
      {
      minDist = tf;
      minDistI = (int)(kPnt);
      }
    }
  if( this->GetDebug() )
    {
    std::cout << "Found point i = " << minDistI << std::endl;
    }

  int kernPnt = minDistI;
  if( this->GetDebug() )
    {
    std::cout << "Calling MeasuresInFullKernelArray" << std::endl;
    }
  this->MeasuresInFullKernelArray( kernArray, kernPnt, len-1, m_RadiusStart );

  if( this->GetDebug() )
    {
    std::cout << "Calling MeasuresInFullKernelArray2" << std::endl;
    }
  this->MeasuresInFullKernelArray( kernArray, kernPnt, 0,
    kernArray[kernPnt].GetRadius() );

  if( this->GetDebug() )
    {
    std::cout << "Smoothing" << std::endl;
    }
  this->SmoothMeasuresInFullKernelArray( kernArray );

  this->ApplyMeasuresInFullKernelArray( tube, kernArray,
    kernArrayTubePointIndex );

  if( this->GetDebug() )
    {
    std::cout << "Radius results:" << std::endl;
    pntIter = tube->GetPoints().begin();
    while( pntIter != tube->GetPoints().end() )
      {
      std::cout << "   " << pntIter->GetID() << " : " << pntIter->GetRadius()
        << std::endl;
      ++pntIter;
      }
    }

  return true;
}

/**
 * Idle callback */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetIdleCallBack( bool ( *idleCallBack )() )
{
  m_IdleCallBack = idleCallBack;
}

/**
 * Status Call back */
template< class TInputImage >
void
RadiusExtractor<TInputImage>
::SetStatusCallBack( void ( *statusCallBack )( const char *, const char *,
    int ) )
{
  m_StatusCallBack = statusCallBack;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeRadiusExtractor_hxx)
