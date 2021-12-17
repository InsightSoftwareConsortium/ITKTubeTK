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

#ifndef __itktubeBlurImageFunction_hxx
#define __itktubeBlurImageFunction_hxx


#include <itkContinuousIndex.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>

#include <cmath>
#include <algorithm>

namespace itk
{

namespace tube
{

/**
 * Set the input Image */
template< class TInputImage >
BlurImageFunction<TInputImage>
::BlurImageFunction( void )
{
  this->m_Image = NULL;

  m_Spacing.Fill( 0 );
  m_OriginalSpacing.Fill( 0 );
  m_UseRelativeSpacing = true;

  m_Scale = 1;
  m_Extent = 3.1;

  m_KernelTotal = 0;
  m_KernelMin.Fill( 0 );
  m_KernelMax.Fill( 0 );
  m_KernelSize.Fill( 0 );
  m_KernelX.clear();
  m_KernelWeights.clear();

  m_ImageIndexMin.Fill( 0 );
  m_ImageIndexMax.Fill( 0 );
}

/**
 * Set the input Image */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::SetInputImage( const InputImageType * ptr )
{
  this->Superclass::SetInputImage( ptr );
  m_Spacing  = this->GetInputImage()->GetSpacing();
  m_OriginalSpacing  = this->GetInputImage()->GetSpacing();
  if( m_UseRelativeSpacing )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      m_Spacing[i] = m_OriginalSpacing[i] / m_OriginalSpacing[0];
      }
    }

  m_ImageIndexMin =
     this->GetInputImage()->GetLargestPossibleRegion().GetIndex();
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    m_ImageIndexMax[i] = m_ImageIndexMin[i]
      + this->GetInputImage()->GetLargestPossibleRegion().GetSize()[i]
      - 1;
    }

  /* Values by default */
  this->RecomputeKernel();

}


/**
 * Print */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::SetUseRelativeSpacing( bool useRelativeSpacing )
{
  m_UseRelativeSpacing = useRelativeSpacing;
  if( m_UseRelativeSpacing )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      m_Spacing[i] = m_OriginalSpacing[i] / m_OriginalSpacing[0];
      }
    }
  else
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      m_Spacing[i] = m_OriginalSpacing[i];
      }
    }

}

/**
 * Print */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "calculate Blurring value at point:" << std::endl;
  os << indent << "UseRelativeSpacing = " << m_UseRelativeSpacing
    << std::endl;
  os << indent << "Spacing = " << m_Spacing << std::endl;
  os << indent << "OriginalSpacing = " << m_OriginalSpacing << std::endl;
  os << indent << "Scale = " << m_Scale << std::endl;
  os << indent << "Extent = " << m_Extent << std::endl;

  os << indent << "KernelWeights.size = " << m_KernelWeights.size()
    << std::endl;
  os << indent << "KernelX.size = " << m_KernelX.size() << std::endl;
  os << indent << "KernelMin = " << m_KernelMin << std::endl;
  os << indent << "KernelMax = " << m_KernelMax << std::endl;
  os << indent << "KernelSize = " << m_KernelSize << std::endl;
  os << indent << "KernelTotal = " << m_KernelTotal << std::endl;

  os << indent << "ImageIndexMin = " << m_ImageIndexMin << std::endl;
  os << indent << "ImageIndexMax = " << m_ImageIndexMax << std::endl;
}


/**
 * SetScale
 * Pre-compute kernel weights */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::RecomputeKernel( void )
{
  if( this->GetDebug() )
    {
    std::cout << "RecomputeKernel" << std::endl;
    }
  double gfact = -0.5/( m_Scale*m_Scale );

  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    m_KernelMax[i] = ( int )( ( m_Scale*m_Extent )/m_Spacing[i] );
    if( m_KernelMax[i]<1 )
      {
      m_KernelMax[i] = 1;
      }
    m_KernelMin[i] = -m_KernelMax[i];
    m_KernelSize[i] = m_KernelMax[i] - m_KernelMin[i] + 1;
    }
  if( this->GetDebug() )
    {
    std::cout << "  Scale = " << m_Scale << std::endl;
    std::cout << "  Extent = " << m_Extent << std::endl;
    std::cout << "  KernelMin = " << m_KernelMin << std::endl;
    std::cout << "  KernelMax = " << m_KernelMax << std::endl;
    std::cout << "  KernelSize = " << m_KernelSize << std::endl;
    }

  m_KernelWeights.clear();
  m_KernelX.clear();

  IndexType index;
  m_KernelTotal = 0;
  if( ImageDimension == 3 )
    {
    for( index[2] = m_KernelMin[2]; index[2]<=m_KernelMax[2]; index[2]++ )
      {
      double distZ = index[2] * m_Spacing[2];
      distZ = distZ * distZ;
      for( index[1] = m_KernelMin[1]; index[1]<=m_KernelMax[1]; index[1]++ )
        {
        double distY = index[1] * m_Spacing[1];
        distY = distY * distY + distZ;
        for( index[0] = m_KernelMin[0]; index[0]<=m_KernelMax[0]; index[0]++ )
          {
          double dist = index[0] * m_Spacing[0];
          dist = dist * dist + distY;
          double w = std::exp( gfact*( dist ) );
          m_KernelWeights.push_back( w );
          m_KernelX.push_back( index );
          m_KernelTotal += w;
          }
        }
      }
    }
  else if( ImageDimension == 2 )
    {
    for( index[1] = m_KernelMin[1]; index[1]<=m_KernelMax[1]; index[1]++ )
      {
      double distY = index[1] * m_Spacing[1];
      distY = distY * distY;
      for( index[0] = m_KernelMin[0]; index[0]<=m_KernelMax[0]; index[0]++ )
        {
        double dist = index[0] * m_Spacing[0];
        dist = dist * dist + distY;
        double w = std::exp( gfact*( dist ) );
        m_KernelWeights.push_back( w );
        m_KernelX.push_back( index );
        m_KernelTotal += w;
        }
      }
    }
}

/**
 * SetScale
 * Pre-compute kernel weights */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::SetScale( double scale )
{
  if( m_Scale != scale )
    {
    m_Scale = scale;
    this->RecomputeKernel();
    }
}

/**
 * SetExtent
 * Pre-compute kernel weights */
template< class TInputImage >
void
BlurImageFunction<TInputImage>
::SetExtent( double extent )
{
  if( m_Extent != extent )
    {
    m_Extent = extent;
    this->RecomputeKernel();
    }
}

/**
 * Evaluate the function at the specified point */
template< class TInputImage >
double
BlurImageFunction<TInputImage>
::Evaluate( const PointType& point ) const
{
  if( this->GetDebug() )
    {
    std::cout << "BlurImageFunction::Evaluate" << std::endl;
    }

  ContinuousIndexType index;
  if( !this->m_Image )
    {
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      index[i] = point[i];
      }
    }
  else
    {
    this->m_Image->TransformPhysicalPointToContinuousIndex( point, index );
    }

  if( this->GetDebug() )
    {
    std::cout << "  Calling EvaluateAtContinuousIndex " << std::endl;
    }
  return this->EvaluateAtContinuousIndex( index );
}

template< class TInputImage >
double
BlurImageFunction<TInputImage>
::EvaluateAtIndex( const IndexType & point ) const
{
  if( this->GetDebug() )
    {
    std::cout << "BlurImageFunction::EvaluateAtIndex" << std::endl;
    std::cout << "  Point = " << point << std::endl;
    }

  if( !this->m_Image )
    {
    return 0.0;
    }

  double res = 0;
  double wTotal = 0;

  IndexType kernelX;

  bool boundary = false;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    if( point[i]+m_KernelMin[i]<m_ImageIndexMin[i]
       || point[i]+m_KernelMax[i]>m_ImageIndexMax[i] )
      {
      boundary = true;
      break;
      }
    }

  if( !boundary )
    {
    itk::ImageRegionConstIterator< InputImageType > imIt( this->m_Image,
      this->m_Image->GetLargestPossibleRegion() );
    imIt.GoToBegin();
    KernelWeightsListType::const_iterator it = m_KernelWeights.begin();
    KernelWeightsListType::const_iterator itEnd = m_KernelWeights.end();
    typename KernelXListType::const_iterator  itX = m_KernelX.begin();
    int skipX = ( *itX )[0];
    while( it != itEnd )
      {
      if( ( *itX )[0] == skipX )
        {
        for( unsigned int i=0; i<ImageDimension; i++ )
          {
          kernelX[i] = point[i] + ( *itX )[i];
          }
        imIt.SetIndex( kernelX );
        }
      res += ( imIt.Get() ) * ( *it );
      ++imIt;
      ++it;
      ++itX;
      ++kernelX[0];
      }
    wTotal = m_KernelTotal;
    }
  else
    {
    if( this->GetDebug() )
      {
      std::cout << "  Boundary point" << std::endl;
      }
    KernelWeightsListType::const_iterator it;
    KernelWeightsListType::const_iterator itEnd;
    typename KernelXListType::const_iterator  itX;
    it = m_KernelWeights.begin();
    itEnd = m_KernelWeights.end();
    itX = m_KernelX.begin();
    wTotal = 0;
    double w;
    while( it != itEnd )
      {
      bool valid = true;
      for( unsigned int i=0; i<ImageDimension; i++ )
        {
        kernelX[i] = point[i] + ( *itX )[i];
        if( kernelX[i] < m_ImageIndexMin[i] ||
            kernelX[i] > m_ImageIndexMax[i] )
          {
          valid = false;
          break;
          }
        }
      if( valid )
        {
        w = *it;
        res += this->m_Image->GetPixel( kernelX ) * w;
        wTotal += w;
        }
      ++it;
      ++itX;
      }
    }

  if( wTotal < *( m_KernelWeights.begin() ) || wTotal == 0 ) 
    {
    return 0;
    }

  if( this->GetDebug() )
    {
    std::cout << "  result = " << res/wTotal << std::endl;
    }
  return res/wTotal;
}

template< class TInputImage >
double
BlurImageFunction<TInputImage>
::EvaluateAtContinuousIndex( const ContinuousIndexType & point ) const
{
  if( this->GetDebug() )
    {
    std::cout << "BlurImageFunction::EvaluateAtContinuousIndex"
      << std::endl;
    std::cout << "  Point = " << point << std::endl;
    }

  if( !this->m_Image )
    {
    return 0.0;
    }

  double w;
  double res = 0;
  double wTotal = 0;
  double gfact = -0.5/( m_Scale*m_Scale );
  double kernrad = m_Scale*m_Extent*m_Scale*m_Extent;

  IndexType kernelX;

  bool boundary = false;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    if( ( int )( point[i] )+m_KernelMin[i]<m_ImageIndexMin[i]
      || ( int )( point[i] )+m_KernelMax[i]>m_ImageIndexMax[i] )
      {
      boundary = true;
      break;
      }
    }

  if( !boundary )
    {
    if( ImageDimension == 3 )
      {
      for( int z = m_KernelMin[2]; z<=m_KernelMax[2]; z++ )
        {
        kernelX[2] = ( int )( point[2] )+z;
        double distZ = ( kernelX[2]-point[2] )*m_Spacing[2];
        distZ = distZ * distZ;
        for( int y = m_KernelMin[1]; y<=m_KernelMax[1]; y++ )
          {
          kernelX[1] = ( int )( point[1] )+y;
          double distY = ( kernelX[1]-point[1] )*m_Spacing[1];
          distY = distY * distY + distZ;
          for( int x = m_KernelMin[0]; x<=m_KernelMax[0]; x++ )
            {
            kernelX[0] = ( int )( point[0] )+x;
            double distX = ( kernelX[0]-point[0] )*m_Spacing[0];
            double dist = distX * distX + distY;
            if( dist <= kernrad )
              {
              w = std::exp( gfact*dist );
              wTotal += w;
              res += this->m_Image->GetPixel( kernelX ) * w;
              }
            }
          }
        }
      }
    else if( ImageDimension == 2 )
      {
      for( int y = m_KernelMin[1]; y<=m_KernelMax[1]; y++ )
        {
        kernelX[1] = ( int )( point[1] )+y;
        double distY = ( kernelX[1]-point[1] )*m_Spacing[1];
        distY = distY * distY;
        for( int x = m_KernelMin[0]; x<=m_KernelMax[0]; x++ )
          {
          kernelX[0] = ( int )( point[0] )+x;
          double distX = ( kernelX[0]-point[0] )*m_Spacing[0];
          double dist = distX * distX + distY;
          if( dist <= kernrad )
            {
            w = std::exp( gfact*dist );
            wTotal += w;
            res += this->m_Image->GetPixel( kernelX ) * w;
            }
          }
        }
      }
    }
  else
    {
    if( this->GetDebug() )
      {
      std::cout << "  Boundary point" << std::endl;
      }
    IndexType minX;
    IndexType maxX;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      minX[i] = std::max( ( int )( ( int )( point[i] )+m_KernelMin[i] ),
        ( int )( m_ImageIndexMin[i] ) );
      maxX[i] = std::min( ( int )( ( int )( point[i] )+m_KernelMax[i] ),
        ( int )( m_ImageIndexMax[i] ) );
      }

    if( ImageDimension == 3 )
      {
      for( kernelX[2] = minX[2]; kernelX[2]<=maxX[2]; kernelX[2]++ )
        {
        double distZ = ( kernelX[2]-point[2] )*m_Spacing[2];
        distZ = distZ * distZ;
        for( kernelX[1] = minX[1]; kernelX[1]<=maxX[1]; kernelX[1]++ )
          {
          double distY = ( kernelX[1]-point[1] )*m_Spacing[1];
          distY = distY * distY + distZ;
          for( kernelX[0] = minX[0]; kernelX[0]<=maxX[0]; kernelX[0]++ )
            {
            double distX = ( kernelX[0]-point[0] )*m_Spacing[0];
            double dist = distX * distX + distY;
            if( dist <= kernrad )
              {
              w = std::exp( gfact*( dist ) );
              wTotal += w;
              res += this->m_Image->GetPixel( kernelX ) * w;
              }
            }
          }
        }
      }
    else if( ImageDimension == 2 )
      {
      for( kernelX[1] = minX[1]; kernelX[1]<=maxX[1]; kernelX[1]++ )
        {
        double distY = ( kernelX[1]-point[1] )*m_Spacing[1];
        distY = distY * distY;
        for( kernelX[0] = minX[0]; kernelX[0]<=maxX[0]; kernelX[0]++ )
          {
          double distX = ( kernelX[0]-point[0] )*m_Spacing[0];
          double dist = distX * distX + distY;
          if( dist <= kernrad )
            {
            w = std::exp( gfact*( dist ) );
            wTotal += w;
            res += this->m_Image->GetPixel( kernelX ) * w;
            }
          }
        }
      }
    }

  if( wTotal < *( m_KernelWeights.begin() ) || wTotal == 0 )
    {
    return 0;
    }
  if( this->GetDebug() )
    {
    std::cout << "  result = " << res/wTotal << std::endl;
    }
  return res/wTotal;
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeBlurImageFunction_hxx )
