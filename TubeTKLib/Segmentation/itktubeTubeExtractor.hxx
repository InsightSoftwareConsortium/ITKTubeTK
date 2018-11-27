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

#ifndef __itktubeTubeExtractor_hxx
#define __itktubeTubeExtractor_hxx

#include "itktubeTubeExtractor.h"

namespace itk
{

namespace tube
{

/**
 * Constructor */
template< class TInputImage >
TubeExtractor<TInputImage>
::TubeExtractor( void )
{
  m_RidgeOp = NULL;
  m_RadiusOp = NULL;

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
  m_NewTubeCallBack = NULL;
  m_AbortProcess = NULL;

  m_InputImage = NULL;
  m_RadiusInputImage = NULL;

  m_TubeColor.set_size( 4 );
  m_TubeColor[0] = 1.0f;
  m_TubeColor[1] = 0.0f;
  m_TubeColor[2] = 0.0f;
  m_TubeColor[3] = 1.0f;

  m_TubeGroup = TubeGroupType::New();
}

/**
 * Destructor */
template< class TInputImage >
TubeExtractor<TInputImage>
::~TubeExtractor( void )
{
}

/**
 * Set the input image */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetInputImage( ImageType * inputImage )
{
  this->m_InputImage = inputImage;
  this->m_RidgeOp = RidgeExtractor<ImageType>::New();
  this->m_RidgeOp->SetInputImage( this->m_InputImage );

  this->m_RadiusInputImage = inputImage;
  this->m_RadiusOp = RadiusExtractor2<ImageType>::New();
  this->m_RadiusOp->SetInputImage( this->m_InputImage );

  this->m_RidgeOp->SetRadiusExtractor( this->m_RadiusOp );
}

/**
 * Optionally set a different image for radius estimation */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetRadiusInputImage( ImageType * inputImage )
{
  if( m_RadiusOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }
  this->m_RadiusInputImage = inputImage;
  this->m_RadiusOp->SetInputImage( this->m_RadiusInputImage );

  this->m_RidgeOp->SetRadiusExtractor( this->m_RadiusOp );
}

/**
 * Set the tube mask image */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetTubeMaskImage( typename
  TubeExtractor<TInputImage>::TubeMaskImageType * mask )
{
  m_RidgeOp->SetTubeMaskImage( mask );
}

/**
 * Get the tube mask image */
template< class TInputImage >
typename TubeExtractor<TInputImage>::TubeMaskImageType *
TubeExtractor<TInputImage>
::GetTubeMaskImage( void )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return m_RidgeOp->GetTubeMaskImage();
}

/**
 * Set Data Min value */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetDataMin( double dataMin )
{
  if( this->m_RidgeOp.IsNull() || this->m_RadiusOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetDataMin( dataMin );
  this->m_RadiusOp->SetDataMin( dataMin );
}

/**
 * Get Data Min value */
template< class TInputImage >
double
TubeExtractor<TInputImage>
::GetDataMin( void )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->GetDataMin();
}

/**
 * Set Data Max value */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetDataMax( double dataMax )
{
  if( this->m_RidgeOp.IsNull() || this->m_RadiusOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetDataMax( dataMax );
  this->m_RadiusOp->SetDataMax( dataMax );
}


/**
 * Get Data Max value */
template< class TInputImage >
double
TubeExtractor<TInputImage>
::GetDataMax( void )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->GetDataMax();
}

/**
 * Set Data Min value */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetExtractBoundMin( const typename TInputImage::IndexType & dataMin )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetExtractBoundMin( dataMin );
}

/**
 * Get Data Min value */
template< class TInputImage >
typename TInputImage::IndexType
TubeExtractor<TInputImage>
::GetExtractBoundMin( void ) const
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->GetExtractBoundMin();
}

/**
 * Set Data Max value */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetExtractBoundMax( const typename TInputImage::IndexType & dataMax )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetExtractBoundMax( dataMax );
}


/**
 * Get Data Max value */
template< class TInputImage >
typename TInputImage::IndexType
TubeExtractor<TInputImage>
::GetExtractBoundMax( void ) const
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->GetExtractBoundMax();
}

/**
 * Set Radius */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetRadius( double radius )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetScale( radius );
  this->m_RadiusOp->SetRadiusStart( radius );
}

/**
 * Get Radius */
template< class TInputImage >
double
TubeExtractor<TInputImage>
::GetRadius( void )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->GetScale();
}

/**
 * Get the ridge extractor */
template< class TInputImage >
typename RidgeExtractor<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRidgeOp( void )
{
  return this->m_RidgeOp;
}

/**
 * Get the radius extractor */
template< class TInputImage >
typename RadiusExtractor2<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRadiusOp( void )
{
  return this->m_RadiusOp;
}

/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template< class TInputImage >
bool
TubeExtractor<TInputImage>
::LocalTube( ContinuousIndexType & x )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->LocalRidge( x );
}

/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template< class TInputImage >
typename TubeExtractor< TInputImage >::TubeType *
TubeExtractor<TInputImage>
::ExtractTube( const ContinuousIndexType & x, unsigned int tubeID,
  bool verbose )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  IndexType xi;
  for( unsigned int i=0; i<ImageDimension; ++i )
    {
    xi[i] = x[i];
    }
  if( this->m_RidgeOp->GetTubeMaskImage()->GetPixel( xi ) != 0 )
    {
    if( this->GetDebug() )
      {
      std::cout << "Initial pixel on prior tube." << std::endl;
      std::cout << "  x = " << x << std::endl;
      }
    return NULL;
    }

  typename TubeType::Pointer tube = this->m_RidgeOp->ExtractRidge( x,
    tubeID, verbose );

  if( tube.IsNull() )
    {
    if( this->GetDebug() )
      {
      std::cout << "m_RidgeOp->Extract() fails!" << std::endl;
      std::cout << "  x = " << x << std::endl;
      }
    return tube;
    }

  if( this->m_AbortProcess != NULL )
    {
    if( this->m_AbortProcess() )
      {
      if( this->m_StatusCallBack )
        {
        this->m_StatusCallBack( "Extract: Ridge", "Aborted", 0 );
        }
      return NULL;
      }
    }

  if( !this->m_RadiusOp->ExtractRadii( tube ) )
    {
    return NULL;
    }

  if( this->m_NewTubeCallBack != NULL )
    {
    this->m_NewTubeCallBack( tube );
    }

  if( this->m_StatusCallBack )
    {
    char s[80];
    std::sprintf( s, "%zd points", tube->GetPoints().size() );
    this->m_StatusCallBack( "Extract: Ridge", s, 0 );
    }

  // Set the Spacing of the tube as the same spacing of the image
  typename ImageType::SpacingType spacing;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    spacing[i] = this->m_InputImage->GetSpacing()[i];
    }
  tube->GetIndexToObjectTransform()->SetScaleComponent( spacing );

  return tube;

}

/**
 * Get list of extracted tubes */
template< class TInputImage >
typename TubeExtractor< TInputImage >::TubeGroupType *
TubeExtractor<TInputImage>
::GetTubeGroup( void )
{
  return m_TubeGroup;
}

/**
 * Set list of extracted tubes */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetTubeGroup( TubeGroupType * tubes )
{
  m_TubeGroup = tubes;
  typename TubeGroupType::ChildrenListType * cList =
    tubes->GetChildren( 9999 );
  typename TubeGroupType::ChildrenListType::iterator iter = cList->begin();
  while( iter != cList->end() )
    {
    this->AddTube( *iter );
    ++iter;
    }
}

/**
 * Smooth a tube */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SmoothTube( TubeType * tube, int h )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  tube = ::tube::SmoothTube< TubeType >( tube, h );
}

/**
 * Add a tube */
template< class TInputImage >
bool
TubeExtractor<TInputImage>
::AddTube( TubeType * tube )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  bool result = this->m_RidgeOp->AddTube( tube );
  if( result )
    {
    m_TubeGroup->AddSpatialObject( tube );
    }

  return result;
}

/**
 * Delete a tube */
template< class TInputImage >
bool
TubeExtractor<TInputImage>
::DeleteTube( TubeType * tube )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  bool result = this->m_RidgeOp->DeleteTube( tube );
  if( result )
    {
    m_TubeGroup->RemoveSpatialObject( tube );
    }

  return result;
}

/**
 * Set the tube color */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::SetTubeColor( const vnl_vector<double> & color )
{
  int nc = color.size();
  if( nc > 4 )
    {
    nc = 4;
    }
  else if( nc < 4 )
    {
    this->m_TubeColor[3] = 1.0;
    }
  for( int i=0; i<nc; i++ )
    {
    this->m_TubeColor[i] = color[i];
    }
}

template< class TInputImage >
vnl_vector<double> &
TubeExtractor<TInputImage>
::GetTubeColor( void )
{
  return m_TubeColor;
}

/**
 * Set the idle call back */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::IdleCallBack( bool ( *idleCallBack )() )
{
  this->m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template< class TInputImage >
void
TubeExtractor<TInputImage>
::StatusCallBack( void ( *statusCallBack )( const char *, const char *,
    int ) )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_StatusCallBack = statusCallBack;
  this->m_RidgeOp->StatusCallBack( statusCallBack );
  this->m_RadiusOp->StatusCallBack( statusCallBack );
}

/**
 * Set the status callback  */
template< class TInputImage >
void
  TubeExtractor<TInputImage>
::NewTubeCallBack( void ( *newTubeCallBack )( TubeType * ) )
{
  this->m_NewTubeCallBack = newTubeCallBack;
}

/**
 * Abort the process  */
template< class TInputImage >
void
  TubeExtractor<TInputImage>
::AbortProcess( bool ( *abortProcess )() )
{
  this->m_AbortProcess = abortProcess;
}

/**
 * PrintSelf */
template< class TInputImage >
void TubeExtractor<TInputImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( this->m_InputImage.IsNotNull() )
    {
    os << indent << "Input Image = " << this->m_InputImage << std::endl;
    }
  else
    {
    os << indent << "Input Image = NULL" << std::endl;
    }

  if( this->m_RadiusInputImage.IsNotNull() )
    {
    os << indent << "Radius Input Image = " << this->m_RadiusInputImage
      << std::endl;
    }
  else
    {
    os << indent << "Radius Input Image = NULL" << std::endl;
    }

  os << indent << "TubeColor.r = " << this->m_TubeColor[0] << std::endl;
  os << indent << "TubeColor.g = " << this->m_TubeColor[1] << std::endl;
  os << indent << "TubeColor.b = " << this->m_TubeColor[2] << std::endl;
  os << indent << "TubeColor.a = " << this->m_TubeColor[3] << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeExtractor_hxx )
