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
#ifndef __itkTubeTubeExtractor_txx
#define __itkTubeTubeExtractor_txx

#include "itkTubeTubeExtractor.h"

namespace itk
{

namespace tube
{

/**
 * Constructor */
template<class TInputImage>
TubeExtractor<TInputImage>
::TubeExtractor( void )
{
  m_RidgeOp = NULL;
  m_RadiusOp = NULL;

  m_IdleCallBack = NULL;
  m_StatusCallBack = NULL;
  m_NewTubeCallBack = NULL;
  m_AbortProcess = NULL;

  m_Color[0] = 0.0f;
  m_Color[1] = 0.0f;
  m_Color[2] = 0.0f;
  m_Color[3] = 0.0f;
}

/**
 * Destructor */
template<class TInputImage>
TubeExtractor<TInputImage>
::~TubeExtractor( void )
{
}

/**
 * Set the input image */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetInputImage( typename ImageType::Pointer inputImage )
{
  this->m_InputImage = inputImage;

  this->m_RidgeOp = RidgeExtractor<ImageType>::New();
  this->m_RidgeOp->SetInputImage( this->m_InputImage );

  this->m_RadiusOp = RadiusExtractor<ImageType>::New();
  this->m_RadiusOp->SetInputImage( this->m_InputImage );
  this->m_RidgeOp->SetRadiusExtractor( this->m_RadiusOp );
}

/**
 * Set the input image */
template<class TInputImage>
typename TubeExtractor<TInputImage>::TubeMaskImageType::Pointer
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
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetDataMin( double dataMin )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetDataMin( dataMin );
  this->m_RadiusOp->SetDataMin( dataMin );
}

/**
 * Get Data Min value */
template<class TInputImage>
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
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetDataMax( double dataMax )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetDataMax( dataMax );
  this->m_RadiusOp->SetDataMax( dataMax );
}


/**
 * Get Data Max value */
template<class TInputImage>
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
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetExtractBoundMin( typename TInputImage::IndexType dataMin )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetExtractBoundMin( dataMin );
}

/**
 * Get Data Min value */
template<class TInputImage>
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
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetExtractBoundMax( typename TInputImage::IndexType dataMax )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetExtractBoundMax( dataMax );
}


/**
 * Get Data Max value */
template<class TInputImage>
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
template<class TInputImage>
void
TubeExtractor<TInputImage>
::SetRadius( double radius )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SetScale( radius );
  this->m_RadiusOp->SetRadius0( radius );
}

/**
 * Get Radius */
template<class TInputImage>
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
 * Extract the ridge */
template<class TInputImage>
void
TubeExtractor<TInputImage>
::ExtractBrightTube( bool extractRidge )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RadiusOp->SetExtractBrightTube( extractRidge );
}

/**
 * Get Extract ridge */
template<class TInputImage>
bool
TubeExtractor<TInputImage>
::ExtractBrightTube( void )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RadiusOp->GetExtractBrightTube();
}

/**
 * Get the ridge extractor */
template<class TInputImage>
typename RidgeExtractor<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRidgeOp( void )
{
  return this->m_RidgeOp;
}

/**
 * Get the radius extractor */
template<class TInputImage>
typename RadiusExtractor<TInputImage>::Pointer
TubeExtractor<TInputImage>
::GetRadiusOp( void )
{
  return this->m_RadiusOp;
}

/**
 * Extract the tube given the position of the first point
 * and the tube ID */
template<class TInputImage>
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
template<class TInputImage>
typename TubeExtractor< TInputImage >::TubeType::Pointer
TubeExtractor<TInputImage>
::ExtractTube( ContinuousIndexType & x, unsigned int tubeID )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  typename TubeType::Pointer tube = this->m_RidgeOp->ExtractRidge( x,
    tubeID );

  if( tube.IsNull() )
    {
    if( this->GetDebug() )
      {
      std::cout << "m_RidgeOp->Extract() fails!" << std::endl;
      std::cout << "  x = " << x << std::endl;
      }
    return tube;
    }

  // Set the Spacing of the tube as the same spacing of the image
  typename ImageType::SpacingType spacing;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    spacing[i] = this->m_InputImage->GetSpacing()[i];
    }
  tube->GetIndexToObjectTransform()->SetScaleComponent( spacing );

  this->m_RidgeOp->SmoothTube( tube, 15 );

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

  double tR = this->m_RadiusOp->GetRadius0();

  if( this->m_RidgeOp->GetDynamicScale() )
    {
    this->m_RadiusOp->SetRadius0( this->m_RidgeOp->GetDynamicScaleUsed() );
    }

  if( !this->m_RadiusOp->ExtractRadii( tube ) )
    {
    this->m_RadiusOp->SetRadius0( tR );
    return NULL;
    }
  this->m_RadiusOp->SetRadius0( tR );

  if( this->m_NewTubeCallBack != NULL )
    {
    this->m_NewTubeCallBack( tube );
    }

  if( this->m_StatusCallBack )
    {
    char s[80];
    sprintf( s, "%ld points", tube->GetPoints().size() );
    this->m_StatusCallBack( "Extract: Ridge", s, 0 );
    }

  return tube;

}

/**
 * Smooth a tube */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::SmoothTube( TubeType * tube, int h )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  this->m_RidgeOp->SmoothTube( tube, h );
}

/**
 * Add a tube */
template<class TInputImage>
bool
  TubeExtractor<TInputImage>
::AddTube( TubeType * tube )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->AddTube( tube );
}

/**
 * Delete a tube */
template<class TInputImage>
bool
  TubeExtractor<TInputImage>
::DeleteTube( TubeType * tube )
{
  if( this->m_RidgeOp.IsNull() )
    {
    throw( "Input data must be set first in TubeExtractor" );
    }

  return this->m_RidgeOp->DeleteTube( tube );
}

/**
 * Set the tube color */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::SetColor( float color[4] )
{
  for( unsigned int i=0; i<4; i++ )
    {
    this->m_Color[i] = color[i];
    }
}

/**
 * Set the idle call back */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::IdleCallBack( bool ( *idleCallBack )() )
{
  this->m_IdleCallBack = idleCallBack;
}

/**
 * Set the status callback  */
template<class TInputImage>
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
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::NewTubeCallBack( void ( *newTubeCallBack )( TubeType * ) )
{
  this->m_NewTubeCallBack = newTubeCallBack;
}

/**
 * Abort the process  */
template<class TInputImage>
void
  TubeExtractor<TInputImage>
::AbortProcess( bool ( *abortProcess )() )
{
  this->m_AbortProcess = abortProcess;
}

/**
 * PrintSelf */
template<class TInputImage>
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

  os << indent << "Color.r = " << this->m_Color[0] << std::endl;
  os << indent << "Color.g = " << this->m_Color[1] << std::endl;
  os << indent << "Color.b = " << this->m_Color[2] << std::endl;
  os << indent << "Color.a = " << this->m_Color[3] << std::endl;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkTubeTubeExtractor_txx)
