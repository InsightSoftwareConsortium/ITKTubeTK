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

#ifndef __itktubeRidgeExtractor_hxx
#define __itktubeRidgeExtractor_hxx

#ifdef WIN32
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include "itktubeRidgeExtractor.h"
#include "tubeMatrixMath.h"

#include <itkImageRegionIterator.h>
#include <itkMinimumMaximumImageFilter.h>
#include <itkNeighborhoodIterator.h>

#include <list>

namespace itk
{

namespace tube
{

template< class TInputImage >
class RidgeExtractorSplineValue
  : public ::tube::UserFunction< vnl_vector< int >, double >
{
public:

  RidgeExtractorSplineValue(
    RidgeExtractor<TInputImage> * newRidgeExtractor )
    {
    m_Ridge = newRidgeExtractor;
    m_XVal = 0;
    m_XIndx.Fill( 0 );
    }

  const double & Value( const vnl_vector<int> & x )
    {
    for( unsigned int i=0; i<x.size(); i++ )
      {
      m_XIndx[i] = x[i];
      }

    m_XVal = m_Ridge->Intensity( m_XIndx );

    return m_XVal;
    }

protected:

  RidgeExtractor<TInputImage> *   m_Ridge;

  typename TInputImage::IndexType m_XIndx;
  double                          m_XVal;

}; // End class RidgeExtractorSplineValue


/**
 * Constructor */
template< class TInputImage >
RidgeExtractor<TInputImage>
::RidgeExtractor( void )
{
  m_DataFunc = BlurImageFunction<ImageType>::New();
  m_DataFunc->SetScale( 3 ); // 1.5
  m_DataFunc->SetExtent( 1.5 ); // 3
  m_DataMin = 0;
  m_DataMax = 1;
  m_DataRange = 1;

  m_StepX = 0.1;
  m_X.set_size( ImageDimension );
  m_X.fill( 0.0 );
  m_XP.set_size( ImageDimension );
  m_XP.fill( 0.0 );
  m_XVal = 0.0;
  m_XD.set_size( ImageDimension );
  m_XD.fill( 0.0 );
  m_XH.set_size( ImageDimension, ImageDimension );
  m_XH.fill( 0.0 );
  m_XHEVal.set_size( ImageDimension );
  m_XHEVal.fill( 0.0 );
  m_XHEVect.set_size( ImageDimension, ImageDimension );
  m_XHEVect.fill( 0.0 );
  m_XRidgeness = 0;
  m_XRoundness = 0;
  m_XCurvature = 0;
  m_XLevelness = 0;

  m_DynamicScale = true;
  m_DynamicScaleUsed = 3;
  m_DynamicStepSize = false;
  m_RadiusExtractor = NULL;

  m_ExtractBoundMin.Fill( 0 );
  m_ExtractBoundMax.Fill( 0 );

  m_MaxTangentChange = 0.5;
  m_MaxXChange = 3.0;
  m_MinRidgeness = 0.70;    // near 1 = harder
  m_MinRidgenessStart = 0.69;
  m_MinRoundness = 0.15;    // near 1 = harder
  m_MinRoundnessStart = 0.1;
  m_MinCurvature = 0.1;
  m_MinCurvatureStart = 0.05;
  m_MinLevelness = 0.5;
  m_MinLevelnessStart = 0.45;
  m_MaxRecoveryAttempts = 4;

  m_SplineValueFunc = new RidgeExtractorSplineValue<TInputImage>( this );
  m_DataSpline = new ::tube::SplineND( ImageDimension, m_SplineValueFunc,
    &m_DataSpline1D, &m_DataSplineOpt );

  m_DataSpline->SetClip( true );

  m_DataSpline->GetOptimizerND()->SetSearchForMin( false );
  m_DataSpline->GetOptimizerND()->SetTolerance( 0.001 );
  m_DataSpline->GetOptimizerND()->SetMaxIterations( 200 );
  m_DataSpline->GetOptimizerND()->SetMaxLineSearches( 10 );
  vnl_vector< double > xStep( ImageDimension, 0.1 );
  m_DataSpline->GetOptimizerND()->SetXStep( xStep );

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;

  m_CurrentFailureCode = SUCCESS;
  m_FailureCodeCount.set_size( this->GetNumberOfFailureCodes() );
  m_FailureCodeCount.fill( 0 );

  m_Tube = NULL;
}

/**
 * Destructor */
template< class TInputImage >
RidgeExtractor<TInputImage>
::~RidgeExtractor( void )
{
  if( m_SplineValueFunc != NULL )
    {
    delete m_SplineValueFunc;
    }
  m_SplineValueFunc = NULL;

  if( m_DataSpline != NULL )
    {
    delete m_DataSpline;
    }
  m_DataSpline = NULL;

}

/**
 * Set the input image */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetInputImage( typename ImageType::Pointer inputImage )
{
  if( this->GetDebug() )
    {
    std::cout << std::endl << "Ridge::SetInputImage" << std::endl;
    }

  m_InputImage = inputImage;

  if( m_InputImage.IsNotNull() )
    {
    m_DataFunc->SetUseRelativeSpacing( true );
    m_DataFunc->SetInputImage( inputImage );

    typedef MinimumMaximumImageFilter<ImageType> MinMaxFilterType;
    typename MinMaxFilterType::Pointer minMaxFilter =
      MinMaxFilterType::New();
    minMaxFilter->SetInput( m_InputImage );
    minMaxFilter->Update();
    m_DataMin = minMaxFilter->GetMinimum();
    m_DataMax = minMaxFilter->GetMaximum();
    m_DataRange = m_DataMax-m_DataMin;

    if( this->GetDebug() )
      {
      std::cout << "  Data Minimum = " << m_DataMin << std::endl;
      std::cout << "  Data Maximum = " << m_DataMax << std::endl;
      std::cout << "  Data Range = " << m_DataRange << std::endl;
      }

    typename ImageType::RegionType region;
    region = m_InputImage->GetLargestPossibleRegion();
    vnl_vector<int> vMin( ImageDimension );
    vnl_vector<int> vMax( ImageDimension );
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      m_ExtractBoundMin[i] = region.GetIndex()[i];
      m_ExtractBoundMax[i] = m_ExtractBoundMin[i] + region.GetSize()[i]-1;
      vMin[i] = m_ExtractBoundMin[i];
      vMax[i] = m_ExtractBoundMax[i];
      }
    m_DataSpline->SetXMin( vMin );
    m_DataSpline->SetXMax( vMax );

    if( this->GetDebug() )
      {
      std::cout << "  Origin = " << m_InputImage->GetOrigin() << std::endl;
      std::cout << "  Dim Minimum = " << m_ExtractBoundMin << std::endl;
      std::cout << "  Dim Maximum = " << m_ExtractBoundMax << std::endl;
      }

    /** Allocate the mask image */
    m_TubeMaskImage = TubeMaskImageType::New();
    m_TubeMaskImage->SetRegions( region );
    m_TubeMaskImage->CopyInformation( m_InputImage );
    m_TubeMaskImage->Allocate();
    m_TubeMaskImage->FillBuffer( 0 );

    } // end Image == NULL
}

/**
 * Get the input image */
template< class TInputImage >
typename TInputImage::Pointer
RidgeExtractor<TInputImage>
::GetInputImage( void )
{
  return m_InputImage;
}

/**
 * Set Data Min value */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetDataMin( double dataMin )
{
  m_DataMin = dataMin;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set Data Min value */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetDataMax( double dataMax )
{
  m_DataMax = dataMax;
  m_DataRange = m_DataMax-m_DataMin;
}

/**
 * Set the scale */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetScale( double scale )
{
  if( this->GetDebug() )
    {
    std::cout << "Ridge::SetScale = " << scale << std::endl;
    }
  m_DataSpline->SetNewData( true );
  m_DataFunc->SetScale( scale );
}

/**
 * Get the scale */
template< class TInputImage >
double
RidgeExtractor<TInputImage>
::GetScale( void )
{
  return m_DataFunc->GetScale();
}

/**
 * Set the extent */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetScaleKernelExtent( double extent )
{
  m_DataSpline->SetNewData( true );
  m_DataFunc->SetExtent( extent );
}


/**
 * Get the extent */
template< class TInputImage >
double
RidgeExtractor<TInputImage>
::GetScaleKernelExtent( void )
{
  return m_DataFunc->GetExtent();
}

/**
 * Get the data spline */
template< class TInputImage >
::tube::SplineND *
RidgeExtractor<TInputImage>
::GetDataSpline( void )
{
  return m_DataSpline;
}

/**
 * Get the data spline 1D */
template< class TInputImage >
::tube::Spline1D *
RidgeExtractor<TInputImage>
::GetDataSpline1D( void )
{
  return & m_DataSpline1D;
}

/**
 * Get the data spline optimizer */
template< class TInputImage >
::tube::Optimizer1D *
RidgeExtractor<TInputImage>
::GetDataSplineOptimizer( void )
{
  return & m_DataSplineOpt;
}

/**
 * Set the dynamic scale */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetDynamicScale( bool dynamicScale )
{
  this->m_DynamicScale = dynamicScale;
}

/**
 * Set the dynamic step size */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetDynamicStepSize( bool dynamicStepSize )
{
  this->m_DynamicStepSize = dynamicStepSize;
}


/**
 * Set the radius extractor */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::SetRadiusExtractor( RadiusExtractor2<TInputImage> * radiusExtractor )
{
  m_RadiusExtractor = radiusExtractor;
}

/**
 * Return the intensity */
template< class TInputImage >
double
RidgeExtractor<TInputImage>
::Intensity( const IndexType & x )
{
  double tf = ( m_DataFunc->EvaluateAtIndex( x )-m_DataMin )/m_DataRange;

  if( tf<0 )
    {
    tf = 0;
    }
  else if( tf>1 )
    {
    tf = 1;
    }

  return tf;
}

/**
 * Ridgeness
 */
template< class TInputImage >
double
RidgeExtractor<TInputImage>
::Ridgeness( const ContinuousIndexType & x,
  double & intensity,
  double & roundness,
  double & curvature,
  double & levelness,
  const vnl_vector<double> & prevTangent )
{
  if( this->GetDebug() )
    {
    std::cout << "Ridge::Ridgeness" << std::endl;
    }

  // update current location - m_X
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    m_X[i] = x[i];
    }

  // compute and update the intensity value, first-derivative,
  // and hessian at m_X
  m_XVal = m_DataSpline->ValueJet( m_X, m_XD, m_XH );

  // test for nan
  if( m_XVal != m_XVal || m_XD[0] != m_XD[0] || m_XH( 0, 0 ) != m_XH( 0, 0 ) )
    {
    std::cerr << "NAN at " << m_X << std::endl;

    intensity = 0;
    roundness = 0;
    curvature = 0;
    levelness = 0;

    m_XRoundness = 0;
    m_XCurvature = 0;
    m_XLevelness = 0;
    m_XRidgeness = 0;

    return m_XRidgeness;
    }

  if( this->GetDebug() )
    {
    std::cout << "  Scale = " << m_DataFunc->GetScale() << std::endl;
    std::cout << "  X = " << m_X << std::endl;
    std::cout << "  XD = " << m_XD << std::endl;
    std::cout << "  XH = " << m_XH << std::endl;
    }

  ::tube::ComputeRidgeness<double>( m_XH, m_XD, prevTangent,
    m_XRidgeness, m_XRoundness, m_XCurvature, m_XLevelness, m_XHEVect,
    m_XHEVal );

  intensity = m_XVal;
  roundness = m_XRoundness;
  curvature = m_XCurvature;
  levelness = m_XLevelness;

  return m_XRidgeness;
}

/**
 * Get Current Location
 */
template< class TInputImage >
const typename RidgeExtractor<TInputImage>::VectorType &
RidgeExtractor<TInputImage>::GetCurrentLocation() const
{
  return m_X;
}

/**
 * Get the Hessian Eigen Basis at the Current Location
 */
template< class TInputImage >
const typename RidgeExtractor<TInputImage>::MatrixType &
RidgeExtractor<TInputImage>::GetCurrentBasis() const
{
  return m_XHEVect;
}

/**
 * Get Intensity at the Current Location
 */
template< class TInputImage >
double RidgeExtractor<TInputImage>::GetCurrentIntensity() const
{
  return m_XVal;
}

/**
 * Get Ridgeness at the Current Location
 */
template< class TInputImage >
double RidgeExtractor<TInputImage>::GetCurrentRidgeness() const
{
  return m_XRidgeness;
}

/**
 * Get Roundess at the Current Location
 */
template< class TInputImage >
double RidgeExtractor<TInputImage>::GetCurrentRoundness() const
{
  return m_XRoundness;
}

/**
 * Get Curvature at the Current Location
 */
template< class TInputImage >
double RidgeExtractor<TInputImage>::GetCurrentCurvature() const
{
  return m_XCurvature;
}

/**
 * Get Levelness at the Current Location
 */
template< class TInputImage >
double RidgeExtractor<TInputImage>::GetCurrentLevelness() const
{
  return m_XLevelness;
}

/**
 * PrintSelf
 */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_InputImage.IsNotNull() )
    {
    os << indent << "Image = " << m_InputImage << std::endl;
    }
  else
    {
    os << indent << "Image = NULL" << std::endl;
    }
  if( m_TubeMaskImage.IsNotNull() )
    {
    os << indent << "DataMask = " << m_TubeMaskImage << std::endl;
    }
  else
    {
    os << indent << "DataMask = NULL" << std::endl;
    }
  if( m_DataFunc.IsNotNull() )
    {
    os << indent << "DataFunc = " << m_DataFunc << std::endl;
    }
  else
    {
    os << indent << "DataFunc = NULL" << std::endl;
    }
  os << std::endl;

  if( m_DynamicScale )
    {
    os << indent << "DynamicScale = True" << std::endl;
    }
  else
    {
    os << indent << "DynamicScale = False" << std::endl;
    }
  os << indent << "DynamicScaleUsed = " << m_DynamicScaleUsed << std::endl;
  if( m_DynamicStepSize )
    {
    os << indent << "DynamicStepSize = True" << std::endl;
    }
  else
    {
    os << indent << "DynamicStepSize = False" << std::endl;
    }
  os << indent << "RadiusExtractor = " << m_RadiusExtractor << std::endl;
  os << indent << "MaxRecoveryAttempts = " << m_MaxRecoveryAttempts
    << std::endl;
  os << indent << "DataMin = " << m_DataMin << std::endl;
  os << indent << "DataMax = " << m_DataMax << std::endl;
  os << indent << "DataRange = " << m_DataRange << std::endl;
  os << indent << "StepX = " << m_StepX << std::endl;
  os << indent << "MaxTangentChange = " << m_MaxTangentChange << std::endl;
  os << indent << "MaxXChange = " << m_MaxXChange << std::endl;
  os << indent << "ExtractBoundMin = " << m_ExtractBoundMin << std::endl;
  os << indent << "ExtractBoundMax = " << m_ExtractBoundMax << std::endl;
  os << indent << "DataSpline1D = " << m_DataSpline1D << std::endl;
  os << indent << "DataSplineOpt = " << m_DataSplineOpt << std::endl;
  os << indent << "DataSpline = " << m_DataSpline << std::endl;
  os << indent << "SplineValueFunc = " << m_SplineValueFunc << std::endl;
  os << indent << "MinRidgeness = " << m_MinRidgeness << std::endl;
  os << indent << "MinRidgenessStart = " << m_MinRidgenessStart
    << std::endl;
  os << indent << "MinRoundness = " << m_MinRoundness << std::endl;
  os << indent << "MinRoundnessStart = " << m_MinRoundnessStart
    << std::endl;
  os << indent << "MinCurvature = " << m_MinCurvature << std::endl;
  os << indent << "MinCurvatureStart = " << m_MinCurvatureStart
    << std::endl;
  os << indent << "MinLevelness = " << m_MinLevelness << std::endl;
  os << indent << "MinLevelnessStart = " << m_MinLevelnessStart
    << std::endl;
  os << std::endl;

  os << indent << "X = " << m_X << std::endl;
  os << indent << "XP = " << m_XP << std::endl;
  os << indent << "XVal = " << m_XVal << std::endl;
  os << indent << "XD = " << m_XD << std::endl;
  os << indent << "XH = " << m_XH << std::endl;
  os << indent << "XHEVal = " << m_XHEVal << std::endl;
  os << indent << "XHEVect = " << m_XHEVect << std::endl;
  os << indent << "XRidgeness = " << m_XRidgeness << std::endl;
  os << indent << "XRoundness = " << m_XRoundness << std::endl;
  os << indent << "XCurvature = " << m_XCurvature << std::endl;
  os << indent << "XLevelness = " << m_XLevelness << std::endl;
  os << std::endl;

  if( m_Tube.IsNotNull() )
    {
    os << indent << "Tube = " << m_Tube << std::endl;
    }
  else
    {
    os << indent << "Tube = NULL" << std::endl;
    }

  os << indent << "IdleCallBack = " << m_IdleCallBack << std::endl;
  os << indent << "StatusCallBack = " << m_StatusCallBack << std::endl;
}

/**
 * Traverse one way */
template< class TInputImage >
bool
RidgeExtractor<TInputImage>
::TraverseOneWay( ContinuousIndexType & newX, VectorType & newT,
                  MatrixType & newN, int dir, bool verbose )
{
  if( this->GetDebug() )
    {
    std::cout << "Ridge::TraverseOneWay" << std::endl;
    }

  ContinuousIndexType indxX;

  double      lVal;
  VectorType  lX( ImageDimension );
  VectorType  lT( ImageDimension );
  MatrixType  lN( ImageDimension, ImageDimension-1 );
  VectorType  lNTEVal( ImageDimension );
  VectorType  lStepDir( ImageDimension );
  MatrixType  lSearchDir( ImageDimension, ImageDimension-1 );

  VectorType  pX( ImageDimension );
  VectorType  pT( ImageDimension );
  MatrixType  pN( ImageDimension, ImageDimension-1 );
  VectorType  pStepDir( ImageDimension );
  MatrixType  pSearchDir( ImageDimension, ImageDimension-1 );

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    lX[i] = newX[i];
    lT[i] = newT[i];
    lStepDir[i] = newT[i];
    }
  for( unsigned int j=0; j<ImageDimension-1; j++ )
    {
    lN.set_column( j, newN.get_column( j ) );
    lSearchDir.set_column( j, newN.get_column( j ) );
    }

  pX = lX;
  pT = lT;
  pN = lN;
  pStepDir = lStepDir;
  pSearchDir = lSearchDir;

  int tubeId = m_Tube->GetId();
  int tubePointCount = m_Tube->GetPoints().size();
  int tubePointCountStart = tubePointCount;

  VectorType prod( ImageDimension );
  VectorType tV( ImageDimension );

  typename ImageType::IndexType indx;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indx[i] = ( int )( lX[i]+0.5 );
    if( lX[i] < ( double )m_ExtractBoundMin[i]
      || lX[i]+0.5 > ( double )m_ExtractBoundMax[i] )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "Ridge: TraverseOneWay: Exited boundary" << std::endl;
        }
      m_CurrentFailureCode = EXITED_IMAGE;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      return false;
      }
    }

  double intensity;
  double ridgeness;
  double roundness;
  double curvature;
  double levelness;

  std::vector< TubePointType > pnts;
  pnts.clear();

  typename TubeMaskImageType::PixelType value =
    m_TubeMaskImage->GetPixel( indx );
  if( value != 0 && ( int )value != tubeId )
    {
    if( verbose || this->GetDebug() )
      {
      std::cout << "Ridge: TraverseOneWay: Encountered another tube"
        << std::endl;
      }
    m_CurrentFailureCode = REVISITED_VOXEL;
    ++m_FailureCodeCount[ m_CurrentFailureCode ];
    return false;
    }
  else
    {
    m_TubeMaskImage->SetPixel( indx, ( float )( tubeId
      + ( tubePointCount/10000.0 ) ) );
    if( dir == 1 )
      {
      if( this->GetDebug() )
        {
        std::cout << "Ridge: dir = 1" << std::endl;
        }
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        indxX[i] = lX[i];
        }
      if( this->GetDebug() )
        {
        std::cout << "Initial point ridgeness..." << std::endl;
        }
      ridgeness = Ridgeness( indxX, intensity, roundness, curvature,
        levelness );

      if( this->GetDebug() )
        {
        std::cout << "Adding point..." << std::endl;
        }
      TubePointType pnt;
      pnt.SetID( tubePointCount );
      typename TubePointType::PointType tubeX;
      typename TubePointType::CovariantVectorType tubeN;
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        tubeX[i] = lX[i];
        if( i < ImageDimension-1 )
          {
          for( unsigned int j=0; j<ImageDimension; j++ )
            {
            tubeN[j] = m_XHEVect( j, i );
            }
          if( i == 0 )
            {
            pnt.SetNormal1( tubeN );
            pnt.SetAlpha1( m_XHEVal[0] );
            if( this->GetDebug() )
              {
              std::cout << " n1 = " << tubeN << std::endl;
              }
            }
          else if( i == 1 )
            {
            pnt.SetNormal2( tubeN );
            pnt.SetAlpha2( m_XHEVal[1] );
            if( this->GetDebug() )
              {
              std::cout << " n2 = " << tubeN << std::endl;
              }
            }
          }
        }
      pnt.SetPosition( tubeX );
      if( m_XHEVal[0] != 0 )
        {
        pnt.SetRidgeness( ridgeness );
        }
      else
        {
        pnt.SetRidgeness( 0.0 );
        }
      pnts.push_back( pnt );
      tubePointCount++;

      if( verbose || this->GetDebug() )
        {
        std::cout << "Ridge: Added initial tube point." << std::endl;
        }
      }
    }

  double iScale0 = GetScale();

  double stepFactor0 = 1;
  double stepFactor = stepFactor0;

  int recovery = 0;
  int prevRecoveryPoint = tubePointCount;
  while( recovery < m_MaxRecoveryAttempts &&
    prevRecoveryPoint+( 2.0/m_StepX ) > tubePointCount )
    {
    if( recovery > 0 )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "Attempting recovery : " << recovery;
        std::cout << " : Scale = " << GetScale() << std::endl;
        std::cout << "   x = " << pX << std::endl;
        }
      switch( recovery )
        {
        default:
        case 1:
          prevRecoveryPoint = tubePointCount;
          stepFactor = 1.5 * stepFactor0;
          break;
        case 2:
          stepFactor = 0.5 * stepFactor0;
          break;
        case 3:
          stepFactor = 2.0 * stepFactor0;
          SetScale( GetScale() * 1.25 );
          break;
        }
      if( this->GetDebug() )
        {
        std::cout << "   Point = " << tubePointCount
          << ": Recovery: new scale = " << GetScale() << std::endl;
        }
      lX = pX;
      lT = pT;
      lN = pN;
      lStepDir = pStepDir;
      lSearchDir = pSearchDir;
      }
    else
      {
      if( this->GetDebug() )
        {
        std::cout << "Ridge: No recovery needed" << std::endl;
        }
      if( !m_DynamicScale && GetScale() != iScale0 )
        {
        SetScale( iScale0 + 0.5 * ( GetScale() - iScale0 ) );
        }

      if( dot_product( lStepDir, pStepDir ) <
        m_MaxTangentChange )
        {
        if( verbose || this->GetDebug() )
          {
          std::cout << "Curving" << std::endl;
          }
        stepFactor *= ( double )0.75;
        }
      else
        {
        if( stepFactor<stepFactor0 )
          {
          stepFactor *= 1.25;
          }
        }
      if( stepFactor>stepFactor0 )
        {
        stepFactor = stepFactor0;
        }
      pX = lX;
      pT = lT;
      pN = lN;
      pStepDir = lStepDir;
      pSearchDir = lSearchDir;
      }

    double currentStepX = m_StepX * stepFactor;
    if( m_DynamicStepSize )
      {
      currentStepX *= GetScale();
      }
    if( currentStepX < 0.5 * m_StepX )
      {
      currentStepX = 0.5 * m_StepX;
      }
    else if( currentStepX > 4.0 * m_StepX )
      {
      currentStepX = 4.0 * m_StepX;
      }
    vnl_vector<double> v = ::tube::ComputeLineStep( lX, currentStepX,
      lStepDir );
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      lX[i] = v[i];
      }
    if( this->GetDebug() )
      {
      std::cout << "Ridge: Computed line step = " << v << std::endl;
      std::cout << "Ridge: TraverseOW: lStepDir = "
                << lStepDir << std::endl;
      std::cout << "Ridge: TraverseOW: lX0 = "
                << lX << std::endl;
      std::cout << "Ridge: TraverseOW: lSearchDir1 = "
                << lSearchDir.get_column( 0 ) << std::endl;
      if( ImageDimension > 2 )
        {
        std::cout << "Ridge: TraverseOW: lSearchDir2 = "
                  << lSearchDir.get_column( 1 ) << std::endl;
        }
      std::cout << "Ridge: TraverseOW: lSearchDir1*StepDir = "
                << dot_product( lSearchDir.get_column( 0 ), lStepDir )
                << std::endl;
      }

    if( !m_DataSpline->Extreme( lX, &lVal, ImageDimension-1, lSearchDir ) )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge terminated: Local max not found"
          << std::endl;
        }
      recovery++;
      m_CurrentFailureCode = OTHER_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      continue;
      }

    if( this->GetDebug() )
      {
      std::cout << "Ridge: TraverseOW: lXExtreme = " << lX << std::endl;
      }

    bool inbounds = true;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      if( lX[i] < ( double )m_ExtractBoundMin[i]
        || lX[i] + 0.5 > ( double )m_ExtractBoundMax[i] )
        {
        inbounds = false;
        if( verbose || this->GetDebug() )
          {
          std::cout << "*** Ridge term: Exited extraction bounds"
            << std::endl;
          }
        m_CurrentFailureCode = EXITED_IMAGE;
        ++m_FailureCodeCount[ m_CurrentFailureCode ];
        break;
        }
      }

    if( !inbounds )
      {
      break;
      }

    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      indxX[i] = lX[i];
      }

    ridgeness = Ridgeness( indxX, intensity, roundness, curvature,
      levelness, lStepDir );

    for( unsigned int i=0; i<ImageDimension-1; i++ )
      {
      lN.set_column( i, m_XHEVect.get_column( i ) );
      lNTEVal[i] = m_XHEVal[i];
      }
    lSearchDir = lN;

    lT = m_XHEVect.get_column( ImageDimension-1 );
    lNTEVal[ImageDimension-1] = m_XHEVal[ImageDimension-1];

    double dProd = dot_product( pT, lT );
    if( dProd<0 )
      {
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        lT[i] = lT[i] * -1;
        }
      }

    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      lStepDir[i] = ( pStepDir[i] * 0.25 ) + ( lT[i] * 0.75 );
      }
    lStepDir.normalize();

    dProd = dot_product( lStepDir, pStepDir );
    if( dProd < m_MaxTangentChange )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge term: Rapid change in step direction "
          << "( " << dProd << " )" << std::endl;
        std::cout << "       pStepDir = " << pStepDir << std::endl;
        std::cout << "       StepDir = " << lStepDir << std::endl;
        std::cout << "       Tangent = " << lT << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        }
      m_CurrentFailureCode = TANGENT_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      recovery++;
      continue;
      }

    double diffX = std::sqrt(
      ::tube::ComputeEuclideanDistanceVector( lX, pX ) );
    if( diffX > m_MaxXChange * GetScale() * stepFactor )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge term: Rapid change in spatial location "
          << "( " << diffX << " )" << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        std::cout << "       maxDiffX = " << m_MaxXChange * GetScale()
          * stepFactor << std::endl;
        }
      m_CurrentFailureCode = DISTANCE_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      recovery++;
      continue;
      }

    if( ridgeness < m_MinRidgeness )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge terminated: Local max not a ridge point "
          << "( ridgeness = " << ridgeness << " )" << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        }
      m_CurrentFailureCode = RIDGE_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      if( ridgeness != 0 && curvature != 0 )
        {
        recovery++;
        }
      else
        {
        recovery = m_MaxRecoveryAttempts;
        }
      continue;
      }

    if( curvature < m_MinCurvature )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge terminated: Low curvature "
          << "( " << curvature << " )" << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        }
      m_CurrentFailureCode = CURVE_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      recovery++;
      continue;
      }

    if( levelness < m_MinLevelness )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge terminated: Low levelness "
          << "( " << levelness << " )" << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        }
      m_CurrentFailureCode = LEVEL_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      recovery++;
      continue;
      }

    if( roundness < m_MinRoundness )
      {
      if( verbose || this->GetDebug() )
        {
        std::cout << "*** Ridge terminated: Roundness : Planar point "
                  << "( " << roundness << " )" << std::endl;
        std::cout << "       Intensity = " << intensity << std::endl;
        std::cout << "       Ridgeness = " << ridgeness << std::endl;
        std::cout << "       Roundness = " << roundness << std::endl;
        std::cout << "       Curvature = " << curvature << std::endl;
        std::cout << "       Levelness = " << levelness << std::endl;
        }
      m_CurrentFailureCode = ROUND_FAIL;
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      if( vnl_math_abs( lNTEVal[0] ) )
        {
        recovery++;
        }
      else
        {
        recovery = m_MaxRecoveryAttempts;
        }
      continue;
      }

    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      indx[i] = ( int )( lX[i]+0.5 );
      }
    double maskVal = m_TubeMaskImage->GetPixel( indx );

    if( maskVal != 0 )
      {
      int oldPoint = ( maskVal - ( int )maskVal ) * 10000;
      if( ( int )maskVal != tubeId ||
        ( ( tubePointCount - oldPoint ) > ( 20 / m_StepX )
        && ( tubePointCount - tubePointCountStart ) > ( 20 / m_StepX ) ) )
        {
        m_CurrentFailureCode = REVISITED_VOXEL;
        ++m_FailureCodeCount[ m_CurrentFailureCode ];
        if( verbose || this->GetDebug() )
          {
          std::cout << "*** Ridge terminated: Revisited voxel" << std::endl;
          std::cout << "  indx = " << indx << std::endl;
          std::cout << "  maskVal = " << maskVal << std::endl;
          std::cout << "  tubeId = " << tubeId << std::endl;
          std::cout << "  tubePointCount = " << tubePointCount << std::endl;
          std::cout << "  StepX = " << m_StepX << std::endl;
          }
        break;
        }
      }
    else
      {
      m_TubeMaskImage->SetPixel( indx, ( float )( tubeId
        + ( tubePointCount/10000.0 ) ) );
      }

    /** Show the satus every 50 points */
    if( tubePointCount%50==0 )
      {
      if( m_StatusCallBack )
        {
        char st[80];
        std::sprintf( st, "Point #%d", tubePointCount );
        m_StatusCallBack( NULL, st, 0 );
        }
      }
    if( this->GetDebug() )
      {
      std::cout << "Ridge: TraverseOW: Adding point " << tubePointCount
        << " = " << lX << std::endl;
      std::cout << "       Intensity = " << intensity << std::endl;
      std::cout << "       Ridgeness = " << ridgeness << std::endl;
      std::cout << "       Roundness = " << roundness << std::endl;
      std::cout << "       Curvature = " << curvature << std::endl;
      std::cout << "       Levelness = " << levelness << std::endl;
      std::cout << "       StepX = " << currentStepX << std::endl;
      std::cout << "       Scale = " << this->GetScale()
        << std::endl;
      }

    TubePointType pnt;
    pnt.SetID( tubePointCount );
    typename TubePointType::PointType tubeX;
    typename TubePointType::VectorType tubeV;
    typename TubePointType::CovariantVectorType tubeN;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      tubeX[i] = lX[i];
      if( i < ImageDimension-1 )
        {
        for( unsigned int j=0; j<ImageDimension; j++ )
          {
          tubeN[j] = m_XHEVect( j, i );
          }
        if( i == 0 )
          {
          pnt.SetNormal1( tubeN );
          pnt.SetAlpha1( m_XHEVal[0] );
          }
        else if( i == 1 )
          {
          pnt.SetNormal2( tubeN );
          pnt.SetAlpha2( m_XHEVal[1] );
          }
        }
      }
    pnt.SetPosition( tubeX );
    if( m_XHEVal[0] != 0 )
      {
      pnt.SetRidgeness( roundness );
      }
    else
      {
      pnt.SetRidgeness( 0.0 );
      }
    pnts.push_back( pnt );
    ++tubePointCount;

    recovery = 0;
    prevRecoveryPoint = tubePointCount;

    if( tubePointCount % static_cast< int >( 2 / currentStepX ) == 0 )
      {
      if( m_IdleCallBack )
        {
        m_IdleCallBack();
        }
      if( m_DynamicScale
        && ( tubePointCount % static_cast< int >( 4 / currentStepX ) ) == 0
        && m_RadiusExtractor )
        {
        if( this->GetDebug() )
          {
          std::cout << "Ridge: TraverseOW: DynamicScale" << std::endl;
          }

        TubePointType tmpPoint;
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          tubeV[i] = lStepDir[i];
          }
        tmpPoint.SetTangent( tubeV );
        tmpPoint.SetRadius( this->GetScale() );
        ::tube::ComputeNormalsFromTangent( tmpPoint, tubeV );
        m_RadiusExtractor->SetRadiusStart( this->GetScale() );
        double radiusMin = m_RadiusExtractor->GetRadiusMin();
        double radiusMax = m_RadiusExtractor->GetRadiusMax();
        double radiusStep = m_RadiusExtractor->GetRadiusStep();
        double radiusTolerance = m_RadiusExtractor->GetRadiusTolerance();
        std::vector< TubePointType > points;
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          tubeX[i] = pX[i];
          }
        tmpPoint.SetPosition( tubeX );
        points.push_back( tmpPoint );
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          tubeX[i] = ( pX[i] + lX[i] ) / 2;
          }
        tmpPoint.SetPosition( tubeX );
        points.push_back( tmpPoint );
        points.push_back( pnt );
        if( m_RadiusExtractor->GetPointVectorOptimalRadius( points,
          m_DynamicScaleUsed, radiusMin, radiusMax, radiusStep,
          radiusTolerance ) )
          {
          m_DynamicScaleUsed = ( tmpPoint.GetRadius()
            + m_DynamicScaleUsed ) / 2;
          }
        if( m_StatusCallBack )
          {
          char s[80];
          std::sprintf( s, "Extract: Ridge: DS = %1.1f",
            m_DynamicScaleUsed );
          m_StatusCallBack( s, NULL, 0 );
          }
        else if( this->GetDebug() || verbose )
          {
          std::cout << "Dynamic Scale = " << m_DynamicScaleUsed
            << std::endl;
          }
        SetScale( m_DynamicScaleUsed );
        m_RadiusExtractor->SetRadiusStart( m_DynamicScaleUsed );
        }
      }
    }

  if( dir == -1 )
    {
    std::vector< TubePointType > * curPoints = &( m_Tube->GetPoints() );
    std::vector< TubePointType > newPoints;
    newPoints.clear();
    for( int i=pnts.size()-1; i>=0; --i )
      {
      newPoints.push_back( pnts[i] );
      }
    for( unsigned int i=0; i<curPoints->size(); ++i )
      {
      newPoints.push_back( ( *curPoints )[i] );
      }
    m_Tube->SetPoints( newPoints );
    }
  else
    {
    for( unsigned int i=0; i<pnts.size(); i++ )
      {
      m_Tube->GetPoints().push_back( pnts[i] );
      }
    }

  if( verbose || this->GetDebug() )
    {
    std::cout << "*** Ridge terminated: Cannot recover" << std::endl;
    std::cout << "    Added " << pnts.size() << " points." << std::endl;
    }
  SetScale( iScale0 );

  if( !pnts.empty() )
    {
    return true;
    }
  else
    {
    return false;
    }
}

template< class TInputImage >
unsigned int
RidgeExtractor<TInputImage>
::GetNumberOfFailureCodes( void ) const
{
  return 10;
}

template< class TInputImage >
const std::string
RidgeExtractor<TInputImage>
::GetFailureCodeName( FailureCodeEnum code ) const
{
  switch( code )
    {
    case SUCCESS:
      {
      return "SUCCESS";
      }
    case EXITED_IMAGE:
      {
      return "EXITED_IMAGE";
      }
    case REVISITED_VOXEL:
      {
      return "REVISITED_VOXEL";
      }
    case RIDGE_FAIL:
      {
      return "RIDGE_FAIL";
      }
    case ROUND_FAIL:
      {
      return "ROUND_FAIL";
      }
    case CURVE_FAIL:
      {
      return "CURVE_FAIL";
      }
    case LEVEL_FAIL:
      {
      return "LEVEL_FAIL";
      }
    case TANGENT_FAIL:
      {
      return "TANGENT_FAIL";
      }
    case DISTANCE_FAIL:
      {
      return "DISTANCE_FAIL";
      }
    default:
    case OTHER_FAIL:
      {
      return "OTHER_FAIL";
      }
    }
}

template< class TInputImage >
unsigned int
RidgeExtractor<TInputImage>
::GetFailureCodeCount( FailureCodeEnum code ) const
{
  return m_FailureCodeCount[ code ];
}

template< class TInputImage >
void
RidgeExtractor<TInputImage>
::ResetFailureCodeCounts( void )
{
  m_FailureCodeCount.fill( 0 );
}


/**
 * Compute the local ridge
 */
template< class TInputImage >
typename RidgeExtractor<TInputImage>::FailureCodeEnum
RidgeExtractor<TInputImage>
::LocalRidge( ContinuousIndexType & newX, bool verbose )
{
  //if( this->GetDebug() )
    {
    std::cout << "Ridge::LocalRidge" << std::endl;
    std::cout << "  x = " << newX << std::endl;
    }

  typename ImageType::IndexType indx;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indx[i] = ( int )( newX[i] + 0.5 ); // rounding
    if( newX[i] < ( double )m_ExtractBoundMin[i]
      || newX[i] + 0.5 > ( double )m_ExtractBoundMax[i] )
      {
      if( m_StatusCallBack )
        {
        m_StatusCallBack( NULL, "Exited Image", 0 );
        }
      if( verbose || this->GetDebug() )
        {
        std::cout << "RidgeExtractor::LocalRidge() : Exited Image 2"
          << std::endl;
        }
      return EXITED_IMAGE;
      }
    }

  double intensity;
  double roundness;
  double curvature;
  double levelness;
  double ridgeness = Ridgeness( newX, intensity, roundness, curvature,
    levelness );
  if( ridgeness == 0 )
    {
    std::cout << "Ridgeness = 0, aborting" << std::endl;
    return RIDGE_FAIL;
    }

  double     val;
  MatrixType lN( ImageDimension, ImageDimension-1 );
  VectorType pX( ImageDimension );
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    pX[ i ] = newX[ i ];
    }

  // Switch to searching within the normal plane
  for( unsigned int loop=0; loop<ImageDimension; loop++ )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      for( unsigned int j=0; j<ImageDimension-1; j++ )
        {
        lN( i, j ) = m_XHEVect( i, j );
        }
      }

    // Local 1D Ridge
    if( this->GetDebug() )
      {
      std::cout << "LocalRidge: Start px = " << pX << std::endl;
      std::cout << "  lN = " << lN << std::endl;
      std::cout << "  val = " << m_DataSpline->Value( pX ) << std::endl;
      }
    m_DataSpline->Extreme( pX, &val, ImageDimension-1, lN );
    if( this->GetDebug() )
      {
      std::cout << "End px = " << pX << std::endl;
      std::cout << "  val = " << val << std::endl;
      }
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      newX[i] = pX[i];
      }

    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      indx[i]=( int )( newX[i] + 0.5 ); // rounding
      if( newX[i] < ( double )m_ExtractBoundMin[i]
        || newX[i] + 0.5 > ( double )m_ExtractBoundMax[i] )
        {
        if( m_StatusCallBack )
          {
          m_StatusCallBack( NULL, "Exited Image", 0 );
          }
        if( verbose || this->GetDebug() )
          {
          std::cout << "RidgeExtractor::LocalRidge() : Exited Image 5"
            << std::endl;
          }
        return EXITED_IMAGE;
        }
      }

    if( m_TubeMaskImage->GetPixel( indx ) != 0 )
      {
      if( m_StatusCallBack )
        {
        m_StatusCallBack( NULL, "Revisited voxel", 0 );
        }
      if( verbose || this->GetDebug() )
        {
        std::cout << "RidgeExtractor::LocalRidge() : Revisited voxel 3"
          << m_TubeMaskImage->GetPixel( indx ) << std::endl;
        }
      return REVISITED_VOXEL;
      }
    if( this->GetDebug() )
      {
      std::cout << "newX = " << newX << std::endl;
      }
    ridgeness = Ridgeness( newX, intensity, roundness, curvature,
      levelness );
    if( ridgeness >= m_MinRidgenessStart &&
      roundness >= m_MinRoundnessStart &&
      curvature >= m_MinCurvatureStart &&
      levelness >= m_MinLevelnessStart )
      {
      if( this->GetDebug() )
        {
        std::cout << " Success: Local norm max: " << newX << std::endl;
        std::cout << "  Ridgeness: " << ridgeness << " >= "
          << m_MinRidgenessStart << std::endl;
        std::cout << "  Roundness: " << roundness << " >= "
          << m_MinRoundnessStart << std::endl;
        std::cout << "  Curvature: " << curvature << " >= "
          << m_MinCurvatureStart << std::endl;
        std::cout << "  Levelness: " << levelness << " >= "
          << m_MinLevelnessStart << std::endl;
        }
      return SUCCESS;
      }

    if( this->GetDebug() )
      {
      std::cout << " Local norm max: " << newX << std::endl;
      std::cout << "  Ridgeness: " << ridgeness << " >= "
        << m_MinRidgenessStart << std::endl;
      std::cout << "  Roundness: " << roundness << " >= "
        << m_MinRoundnessStart << std::endl;
      std::cout << "  Curvature: " << curvature << " >= "
        << m_MinCurvatureStart << std::endl;
      std::cout << "  Levelness: " << levelness << " >= "
        << m_MinLevelnessStart << std::endl;
      }
    }

  if( this->GetDebug() )
    {
    std::cout << " FAIL: Local norm max: " << newX << std::endl;
    std::cout << "  Ridgeness: " << ridgeness << " >= "
      << m_MinRidgenessStart << std::endl;
    std::cout << "  Roundness: " << roundness << " >= "
      << m_MinRoundnessStart << std::endl;
    std::cout << "  Curvature: " << curvature << " >= "
      << m_MinCurvatureStart << std::endl;
    std::cout << "  Levelness: " << levelness << " >= "
      << m_MinLevelnessStart << std::endl;
    }

  if( ridgeness < m_MinRidgenessStart )
    {
    if( m_StatusCallBack )
      {
      m_StatusCallBack( NULL, "Ridgeness failure", 0 );
      }
    if( this->GetDebug() )
      {
      std::cout << "LocalRidge : Ridgeness failure" << std::endl;
      }
    return RIDGE_FAIL;
    }

  if( roundness < m_MinRoundnessStart )
    {
    if( m_StatusCallBack )
      {
      m_StatusCallBack( NULL, "Roundness failure", 0 );
      }
    if( this->GetDebug() )
      {
      std::cout << "LocalRidge : Roundness failure" << std::endl;
      }
    return ROUND_FAIL;
    }

  if( curvature < m_MinCurvatureStart )
    {
    if( m_StatusCallBack )
      {
      m_StatusCallBack( NULL, "Curvature failure", 0 );
      }
    if( this->GetDebug() )
      {
      std::cout << "LocalRidge : Curvature failure" << std::endl;
      }
    return CURVE_FAIL;
    }

  if( levelness < m_MinLevelnessStart )
    {
    if( m_StatusCallBack )
      {
      m_StatusCallBack( NULL, "Levelness failure", 0 );
      }
    if( this->GetDebug() )
      {
      std::cout << "LocalRidge : Levelness failure" << std::endl;
      }
    return LEVEL_FAIL;
    }

  return OTHER_FAIL;
}

/**
 * Extract a tube
 */
template< class TInputImage >
typename RidgeExtractor<TInputImage>::TubeType::Pointer
RidgeExtractor<TInputImage>
::ExtractRidge( const ContinuousIndexType & newX, int tubeId,
  bool verbose )
{
  ContinuousIndexType lX;
  lX = newX;

  double scaleOriginal = GetScale();
  double scale0 = scaleOriginal;
  double radiusOriginal = scaleOriginal;
  if( m_RadiusExtractor )
    {
    radiusOriginal = m_RadiusExtractor->GetRadiusStart();
    }

  // Try to find a ridge voxel close by
  m_CurrentFailureCode = LocalRidge( lX, verbose );
  if( m_CurrentFailureCode != SUCCESS )
    {
    ++m_FailureCodeCount[ m_CurrentFailureCode ];
    if( verbose || this->GetDebug() )
      {
      std::cout << "LocalRidge fails at " << lX << std::endl;
      }
    return nullptr;
    }
  if( verbose || this->GetDebug() )
    {
    std::cout << "*** Ridge found at " << lX << std::endl;
    }

  typename TubeMaskImageType::IndexType indx;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    indx[i] = ( int )( lX[i] + 0.5 );
    }
  typename TubeMaskImageType::PixelType value =
    m_TubeMaskImage->GetPixel( indx );
  if( value != 0 && ( int )value != tubeId )
    {
    m_CurrentFailureCode = REVISITED_VOXEL;
    ++m_FailureCodeCount[ m_CurrentFailureCode ];
    return nullptr;
    }

  MatrixType lN( ImageDimension, ImageDimension-1 );
  VectorType lT( ImageDimension );

  if( m_DynamicScale && m_RadiusExtractor )
    {
    TubePointType tmpPoint;
    typename TubePointType::PointType tubeX;
    typename TubePointType::VectorType tubeV;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      tubeX[i] = lX[i];
      lT[i] = m_XHEVect( i, ImageDimension-1 );
      tubeV[i] = lT[i];
      }
    tmpPoint.SetPosition( tubeX );
    tmpPoint.SetTangent( tubeV );
    tmpPoint.SetRadius( scale0 );
    double radiusMin = m_RadiusExtractor->GetRadiusMin();
    double radiusMax = m_RadiusExtractor->GetRadiusMax();
    double radiusStep = m_RadiusExtractor->GetRadiusStep();
    double radiusTolerance = m_RadiusExtractor->GetRadiusTolerance();
    std::vector< TubePointType > points;
    points.push_back( tmpPoint );
    if( !m_RadiusExtractor->GetPointVectorOptimalRadius( points,
      scale0, radiusMin, radiusMax, radiusStep, radiusTolerance ) )
      {
      if( this->GetDebug() && m_StatusCallBack )
        {
        m_StatusCallBack( "Extract: Ridge: AS = ?",
          "Error: Medial Max Not Found", 0 );
        }
      m_DynamicScaleUsed = scale0;
      }
    else
      {
      m_DynamicScaleUsed = ( tmpPoint.GetRadius() + GetScale() ) / 2;
      }

    SetScale( m_DynamicScaleUsed );
    m_RadiusExtractor->SetRadiusStart( m_DynamicScaleUsed );
    if( verbose || this->GetDebug() )
      {
      std::cout << "DynamicScale = " << GetScale() << std::endl;
      std::cout << "  x =  " << lX << std::endl;
      std::cout << "  newX =  " << newX << std::endl;
      }
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      lX[i] = ( lX[i] + newX[i] ) / 2;
      }
    m_CurrentFailureCode = LocalRidge( lX, verbose );
    if( m_CurrentFailureCode != SUCCESS )
      {
      ++m_FailureCodeCount[ m_CurrentFailureCode ];
      if( m_StatusCallBack )
        {
        m_StatusCallBack( "AS Failure", NULL, 0 );
        }
      if( verbose || this->GetDebug() )
        {
        std::cout << "RidgeExtractor:Extract(): AS Failure" << std::endl;
        }
      m_DynamicScaleUsed = scaleOriginal;
      SetScale( scaleOriginal );
      m_RadiusExtractor->SetRadiusStart( radiusOriginal );
      return nullptr;
      }
    scale0 = m_DynamicScaleUsed;
    SetScale( scale0 );
    m_RadiusExtractor->SetRadiusStart( scale0 );
    }

  m_Tube = TubeType::New();
  m_Tube->SetId( tubeId );
  m_Tube->GetPoints().clear();

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    lT[i] = m_XHEVect( i, ImageDimension-1 );
    for( unsigned int j=0; j<ImageDimension-1; j++ )
      {
      lN( i, j ) = m_XHEVect( i, j );
      }
    }
  if( verbose || this->GetDebug() )
    {
    std::cout << "Traversing one way" << std::endl;
    }
  TraverseOneWay( lX, lT, lN, 1, verbose );
  if( verbose || this->GetDebug() )
    {
    std::cout << "End traversing one way" << std::endl;
    }

  SetScale( scale0 );
  if( m_RadiusExtractor )
    {
    m_RadiusExtractor->SetRadiusStart( scale0 );
    }

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    lT[i] = -1 * lT[i];
    }

  if( verbose || this->GetDebug() )
    {
    std::cout << "Traversing the other way" << std::endl;
    }
  TraverseOneWay( lX, lT, lN, -1, verbose );
  if( verbose || this->GetDebug() )
    {
    std::cout << "End traversing the other way" << std::endl;
    }

  // return to user defaults
  SetScale( scaleOriginal );
  if( m_RadiusExtractor )
    {
    m_RadiusExtractor->SetRadiusStart( radiusOriginal );
    }

  if( m_Tube->GetPoints().size() < 2.0/m_StepX )
    {
    std::cout << "Ridge too short, deleting." << std::endl;
    if( m_StatusCallBack )
      {
      m_StatusCallBack( "Extract: Ridge", "Too short", 0 );
      }
    DeleteTube( m_Tube );
    m_Tube = NULL;
    std::cout << "Ridge returning null." << std::endl;
    return m_Tube;
    }

  if( verbose || this->GetDebug() )
    {
    std::cout << "*** Extracted ridge of " << m_Tube->GetPoints().size()
              << " points." << std::endl;
    }

  //
  // Calculate tangents
  //
  if( m_Tube && m_Tube->GetPoints().size() > 0 )
    {
    if( this->GetDebug() )
      {
      std::cout << "Calculating tangents." << std::endl;
      }
    ::tube::ComputeTubeTangentsAndNormals< TubeType >( m_Tube );
    }

  if( m_StatusCallBack )
    {
    char s[80];
    std::sprintf( s, "%d points", ( int )( m_Tube->GetPoints().size() ) );
    m_StatusCallBack( "Extract: Ridge", s, 0 );
    }

  return m_Tube;
}

template< class TInputImage >
template< class TDrawMask >
bool
RidgeExtractor<TInputImage>
::DeleteTube( const TubeType * tube, TDrawMask * drawMask )
{
  typedef typename TDrawMask::PixelType      DrawPixelType;
  typedef NeighborhoodIterator< TDrawMask >  NeighborhoodIteratorType;

  if( tube->GetPoints().size() == 0 )
    {
    return true;
    }

  if( drawMask == NULL )
    {
    drawMask = m_TubeMaskImage;
    }

  DrawPixelType zero = 0;
  typename std::vector< TubePointType >::const_iterator pnt;
  VectorType x( ImageDimension );
  double r;
  for( pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end();
    ++pnt )
    {
    if( this->GetDebug() )
      {
      std::cout << "Del pnt = " << pnt->GetPosition() << std::endl;
      }
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      x[i] = ( int )( ( *pnt ).GetPosition()[i]+0.5 );
      }

    bool inside = true;
    typename ImageType::IndexType indx;
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      indx[i] = ( int )( x[i] + 0.5 );
      if( x[i] < ( double )m_ExtractBoundMin[i]
        || x[i] + 0.5 > ( double )m_ExtractBoundMax[i] )
        {
        inside = false;
        break;
        }
      }

    if( inside )
      {
      drawMask->SetPixel( indx, zero );
      r = ( *pnt ).GetRadius() + 0.5;
      if( r > 1 )
        {
        for( unsigned int i=0; i<ImageDimension; ++i )
          {
          if( indx[i]-r-1 < m_ExtractBoundMin[i]
            || indx[i]-r-1 > m_ExtractBoundMax[i]
            || indx[i]+r+1 < m_ExtractBoundMin[i]
            || indx[i]+r+1 > m_ExtractBoundMax[i] )
            {
            inside = false;
            break;
            }
          }
        typename NeighborhoodIteratorType::RadiusType rad;
        rad.Fill( r );
        NeighborhoodIteratorType it( rad, drawMask,
          drawMask->GetLargestPossibleRegion() );
        double rr = r * r;
        it.SetLocation( indx );
        if( inside )
          {
          for( unsigned int i=0; i<it.Size(); ++i )
            {
            double dist = 0;
            for( unsigned int j=0; j<ImageDimension; j++ )
              {
              double tf = it.GetOffset( i )[j];
              dist += tf * tf;
              }
            if( dist <= rr )
              {
              it.SetPixel( i, zero );
              }
            }
          }
        else
          {
          for( unsigned int i=0; i<it.Size(); ++i )
            {
            double dist = 0;
            for( unsigned int j=0; j<ImageDimension; j++ )
              {
              double tf = it.GetOffset( i )[j];
              dist += tf * tf;
              }
            if( dist <= rr )
              {
              it.SetPixel( i, zero, inside );
              }
            }
          }
        }
      }
    }
  return true;
}

/**
 * Delete a tube */
template< class TInputImage >
bool
RidgeExtractor<TInputImage>
::DeleteTube( const TubeType * tube )
{
  return this->DeleteTube< TubeMaskImageType >( tube, m_TubeMaskImage );
}


/**
 * Add a tube */
template< class TInputImage >
template< class TDrawMask >
bool
RidgeExtractor<TInputImage>
::AddTube( const TubeType * tube, TDrawMask * drawMask )
{
  if( this->GetDebug() )
    {
    std::cout << "*** START: AddTube" << std::endl;
    }

  typedef NeighborhoodIterator< TDrawMask > NeighborhoodIteratorType;

  int tubeId = tube->GetId();
  int tubePointCount = 0;

  if( drawMask == NULL )
    {
    drawMask = m_TubeMaskImage;
    }

  VectorType x( ImageDimension );
  double r;

  typename std::vector< TubePointType >::const_iterator pnt;

  for( pnt = tube->GetPoints().begin(); pnt != tube->GetPoints().end();
    pnt++ )
    {
    if( this->GetDebug() )
      {
      std::cout << "Add pnt = " << pnt->GetPosition() << std::endl;
      }
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      x[i] = ( int )( ( *pnt ).GetPosition()[i]+0.5 );
      }

    bool inside = true;
    typename ImageType::IndexType indx;
    for( unsigned int i=0; i<ImageDimension; ++i )
      {
      indx[i] = ( int )( x[i] + 0.5 );
      if( x[i] < ( double )m_ExtractBoundMin[i]
        || x[i] + 0.5 > ( double )m_ExtractBoundMax[i] )
        {
        inside = false;
        break;
        }
      }

    if( inside )
      {
      drawMask->SetPixel( indx, ( PixelType )( tubeId +
          ( tubePointCount/10000.0 ) ) );
      r = ( *pnt ).GetRadius() + 0.5;
      if( r > 1 )
        {
        inside = true;
        for( unsigned int i=0; i<ImageDimension; ++i )
          {
          if( indx[i]-r-1 < m_ExtractBoundMin[i]
            || indx[i]-r-1 > m_ExtractBoundMax[i]
            || indx[i]+r+1 < m_ExtractBoundMin[i]
            || indx[i]+r+1 > m_ExtractBoundMax[i] )
            {
            inside = false;
            break;
            }
          }
        typename NeighborhoodIteratorType::RadiusType rad;
        rad.Fill( r );
        NeighborhoodIteratorType it( rad, drawMask,
          drawMask->GetLargestPossibleRegion() );
        double rr = r * r;
        it.SetLocation( indx );
        if( inside )
          {
          for( unsigned int i=0; i<it.Size(); ++i )
            {
            double dist = 0;
            for( unsigned int j=0; j<ImageDimension; j++ )
              {
              double tf = it.GetOffset( i )[j];
              dist += tf * tf;
              }
            if( dist <= rr )
              {
              it.SetPixel( i, ( PixelType )( tubeId +
                  ( tubePointCount/10000.0 ) ) );
              }
            }
          }
        else
          {
          for( unsigned int i=0; i<it.Size(); ++i )
            {
            double dist = 0;
            for( unsigned int j=0; j<ImageDimension; j++ )
              {
              double tf = it.GetOffset( i )[j];
              dist += tf * tf;
              }
            if( dist <= rr )
              {
              it.SetPixel( i, ( PixelType )( tubeId +
                  ( tubePointCount/10000.0 ) ), inside );
              }
            }
          }
        }
      }
    tubePointCount++;
    }

  if( this->GetDebug() )
    {
    std::cout << "*** END: AddTube" << std::endl;
    }

  return true;
}

/**
 * Add a tube */
template< class TInputImage >
bool
RidgeExtractor<TInputImage>
::AddTube( const TubeType * tube )
{
  return this->AddTube< TubeMaskImageType >( tube, m_TubeMaskImage );
}

/** Set the idle call back */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::IdleCallBack( bool ( *idleCallBack )() )
{
  m_IdleCallBack = idleCallBack;
}

/** Set the status callback  */
template< class TInputImage >
void
RidgeExtractor<TInputImage>
::StatusCallBack( void ( *statusCallBack )( const char *, const char *,
    int ) )
{
  m_StatusCallBack = statusCallBack;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeRidgeExtractor_hxx )
