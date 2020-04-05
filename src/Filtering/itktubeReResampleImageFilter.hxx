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

#ifndef __itktubeReResampleImageFilter_hxx
#define __itktubeReResampleImageFilter_hxx

#include "itktubeReResampleImageFilter.h"

namespace itk
{

namespace tube
{

/** Constructor */
template< class TPixel, unsigned int VDimension >
ReResampleImageFilter< TPixel, VDimension >::
ReResampleImageFilter( void )
{
  m_MatchImage = nullptr;
  m_Spacing.clear();
  m_Origin.clear();
  m_Index.clear();
  m_ResampleFactor.clear();
  m_MakeIsotropic = false;
  m_MakeHighResIso = false;
  m_Interpolator = "Linear";
  m_LoadTransform = false;
  m_Transform = nullptr;

  m_Input = nullptr;
  m_Output = nullptr;
}

/** GenerateData */
template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
Update( void )
{
  m_Filter = ResampleFilterType::New();

  m_Filter->SetInput( m_Input );
  
  typename ImageType::SpacingType     inSpacing = m_Input->GetSpacing();
  typename ImageType::PointType       inOrigin = m_Input->GetOrigin();
  typename ImageType::SizeType        inSize =
    m_Input->GetLargestPossibleRegion().GetSize();
  typename ImageType::IndexType       inIndex =
    m_Input->GetLargestPossibleRegion().GetIndex();
  typename ImageType::DirectionType   inDirection =
    m_Input->GetDirection();

  typename ImageType::SizeType       outSize;
  typename ImageType::SpacingType    outSpacing;
  typename ImageType::PointType      outOrigin;
  typename ImageType::IndexType      outIndex;
  typename ImageType::DirectionType  outDirection;

  for( unsigned int i = 0; i< VDimension; i++ )
    {
    outSpacing[i] = inSpacing[i];
    outOrigin[i] = inOrigin[i];
    outIndex[i] = inIndex[i];
    outSize[i] = inSize[i];
    }
  outDirection = inDirection;

  if( m_MatchImage )
    {
    outSpacing = m_MatchImage->GetSpacing();
    outOrigin = m_MatchImage->GetOrigin();
    outDirection = m_MatchImage->GetDirection();
    outSize = m_MatchImage->GetLargestPossibleRegion().GetSize();
    outIndex = m_MatchImage->GetLargestPossibleRegion().GetIndex();
    }

  if( m_Origin.size() > 0 )
    {
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outOrigin[i] = m_Origin[i];
      }
    }
  if( m_Index.size() > 0 )
    {
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outIndex[i] = m_Index[i];
      }
    }
  if( m_Spacing.size() > 0 )
    {
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outSpacing[i] = m_Spacing[i];
      }
    }
  else if( m_ResampleFactor.size() > 0 )
    {
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outSpacing[i] = outSpacing[i] / m_ResampleFactor[i];
      }
    }
  else if( m_MakeIsotropic )
    {
    double iso = outSpacing[0];
    for( unsigned int i = 1; i < VDimension - 1; i++ )
      {
      iso += outSpacing[i];
      }
    iso /= ( VDimension - 1 );
    iso += outSpacing[ VDimension - 1 ];
    iso /= 2;
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outSpacing[i] = iso;
      }
    }
  else if( m_MakeHighResIso )
    {
    double iso = outSpacing[0];
    for( unsigned int i = 1; i < VDimension; i++ )
      {
      if( outSpacing[i] < iso )
        {
        iso = outSpacing[i];
        }
      }
    for( unsigned int i = 0; i < VDimension; i++ )
      {
      outSpacing[i] = iso;
      }
    }

  for( unsigned int i = 0; i < VDimension; i++ )
    {
    if( outSpacing[i] <= 0 )
      {
      std::cerr << "ERROR: Illegal or missing output spacing specified."
                << std::endl;
      return;
      }
    }

  if( !m_MatchImage )
    {
    std::vector< double > outResampleFactor;
    outResampleFactor.resize( VDimension );
    for( unsigned int i = 0; i< VDimension; i++ )
      {
      outResampleFactor[i] = inSpacing[i] / outSpacing[i];
      outSize[i] = static_cast<unsigned long>( inSize[i]
                                              * outResampleFactor[i] );
      }
    }

  typedef typename itk::InterpolateImageFunction< ImageType, double >
    InterpType;
  typename InterpType::Pointer interp;
  if( m_Interpolator == "Sinc" )
    {
    typedef typename itk::WindowedSincInterpolateImageFunction<
      ImageType, VDimension >     SincInterpType;
    interp = SincInterpType::New();
    }
  else if( m_Interpolator == "BSpline" )
    {
    typedef typename itk::BSplineInterpolateImageFunction<
      ImageType, double >         BSplineInterpType;
    interp = BSplineInterpType::New();
    }
  else if( m_Interpolator == "NearestNeighbor" )
    {
    typedef typename itk::NearestNeighborInterpolateImageFunction<
      ImageType, double >    NearestNeighborInterpType;
    interp = NearestNeighborInterpType::New();
    }
  else // default = if( interpolator == "Linear" )
    {
    typedef typename itk::LinearInterpolateImageFunction<
      ImageType, double >    LinearInterpType;
    interp = LinearInterpType::New();
    }
  m_Filter->SetInterpolator( interp );

  if( m_LoadTransform )
    {
    m_Filter->SetTransform( m_Transform );
    }

  m_Filter->SetSize( outSize );
  m_Filter->SetOutputStartIndex( outIndex );
  m_Filter->SetOutputOrigin( outOrigin );
  m_Filter->SetOutputSpacing( outSpacing );
  m_Filter->SetOutputDirection( outDirection );
  m_Filter->SetDefaultPixelValue( 0 );

  m_Filter->Update();

  m_Output = m_Filter->GetOutput();
}

template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
SetTransform( TransformType* t )
{
  m_Transform = t;
}

template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
SetSpacing( std::vector<double> s )
{
  m_Spacing = s;
}

template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
SetOrigin( std::vector<double> o )
{
  m_Origin = o;
}

template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
SetIndex( std::vector<int> i )
{
  m_Index = i;
}

template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
SetResampleFactor( std::vector<double> rf )
{
  m_ResampleFactor = rf;
}

/** PrintSelf */
template< class TPixel, unsigned int VDimension >
void
ReResampleImageFilter< TPixel, VDimension >::
PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Match Image = " << m_MatchImage << std::endl;
  os << indent << "Spacing size = " << m_Spacing.size() << std::endl;
  if( m_Spacing.size() > 0 )
    {
    os << indent << "Spacing[0] = " << m_Spacing[0] << std::endl;
    }
  os << indent << "Origin size = " << m_Origin.size() << std::endl;
  if( m_Origin.size() > 0 )
    {
    os << indent << "Origin[0] = " << m_Origin[0] << std::endl;
    }
  os << indent << "Index size = " << m_Index.size() << std::endl;
  if( m_Index.size() > 0 )
    {
    os << indent << "Index[0] = " << m_Index[0] << std::endl;
    }
  os << indent << "ResampleFactor size = " << m_ResampleFactor.size()
    << std::endl;
  if( m_ResampleFactor.size() > 0 )
    {
    os << indent << "ResampleFactor[0] = " << m_ResampleFactor[0]
      << std::endl;
    }
  if( m_MakeIsotropic )
    {
    os << indent << "MakeIsotropic = True" << std::endl;
    }
  else
    {
    os << indent << "MakeIsotropic = False" << std::endl;
    }
  if( m_MakeHighResIso )
    {
    os << indent << "MakeHighResIso = True" << std::endl;
    }
  else
    {
    os << indent << "MakeHighResIso = False" << std::endl;
    }
  os << indent << "Interpolator = " << m_Interpolator << std::endl;
  if( m_LoadTransform )
    {
    os << indent << "LoadTransform = True" << std::endl;
    }
  else
    {
    os << indent << "LoadTransform = False" << std::endl;
    }
  os << indent << "Transform = " << m_Transform << std::endl;

}

} // End namespace tube

} // End namespace itk
#endif // End !defined( __itktubeReResampleImageFilter_hxx )
