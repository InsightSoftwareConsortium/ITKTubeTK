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

#ifndef __itktubeTubeSpatialObjectToDensityImageFilter_hxx
#define __itktubeTubeSpatialObjectToDensityImageFilter_hxx

#include "itktubeTubeSpatialObjectToDensityImageFilter.h"

/** Constructor */
template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
TubeSpatialObjectToDensityImageFilter< TDensityImageType, TRadiusImageType,
                                 TTangentImageType >
::TubeSpatialObjectToDensityImageFilter( void )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_Size[i] = 0;
    m_Spacing[i] = 1;
    }
  m_MaxDensityIntensity = 255;   //NumericTraits<DensityPixelType>::max();
  m_UseSquaredDistance = false;
}

/** Destructor */
template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
TubeSpatialObjectToDensityImageFilter< TDensityImageType, TRadiusImageType,
                                 TTangentImageType>
::~TubeSpatialObjectToDensityImageFilter( void )
{
}

template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
void
TubeSpatialObjectToDensityImageFilter< TDensityImageType, TRadiusImageType,
                                 TTangentImageType>
::SetSpacing( SpacingType s )
{
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_Spacing[i] = s[i];
    }
}


template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
void
TubeSpatialObjectToDensityImageFilter< TDensityImageType, TRadiusImageType,
                                 TTangentImageType >
::Update( void )
{
  if( m_Size[0] == 0 )
    {
    std::cerr << "Error, no size parameters given " << std::endl;
    return;
    }
  try
    {
    TubeGroupPointer tubes = this->GetInputTubeGroup();

    tubes->ComputeFamilyBoundingBox( 99999 );

    typedef TubeSpatialObjectToImageFilter<ImageDimension, DensityImageType,
      RadiusImageType, TangentImageType> FilterType;
    typename FilterType::Pointer tubefilter = FilterType::New();

    //Set record radius value for filter
    tubefilter->SetBuildRadiusImage( true );
    tubefilter->SetBuildTangentImage( true );
    tubefilter->SetUseRadius( true );
    tubefilter->SetInput( tubes );
    tubefilter->SetSize( m_Size );
    tubefilter->SetSpacing( m_Spacing );
    tubefilter->Update();

    typename DanielssonFilterType::Pointer danFilter =
      DanielssonFilterType::New();

    danFilter->SetInput( tubefilter->GetOutput() );

    danFilter->SetUseImageSpacing( true );
    danFilter->SetInputIsBinary( true );

    if( m_UseSquaredDistance )
      {
      danFilter->SetSquaredDistance( true );
      }

    danFilter->Update();

    VectorImagePointer  vectorImage = danFilter->GetVectorDistanceMap();
    m_RadiusMapImage   = tubefilter->GetRadiusImage();
    m_TangentMapImage  = tubefilter->GetTangentImage();
    m_DensityMapImage  = danFilter->GetDistanceMap();

    typedef ImageRegionIterator<VectorImageType>   VectorIteratorType;
    typedef ImageRegionIterator<RadiusImageType>   RadiusIteratorType;

    VectorIteratorType it_vector( vectorImage,
      vectorImage->GetLargestPossibleRegion() );
    RadiusIteratorType it_radius( m_RadiusMapImage,
      m_RadiusMapImage->GetLargestPossibleRegion() );

    it_vector.GoToBegin();
    it_radius.GoToBegin();

    //Use Vector Image to add the radius values
    while( !it_vector.IsAtEnd() )
      {
      VectorPixelType v = it_vector.Value();
      typename DensityImageType::IndexType index = it_vector.GetIndex();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        index[i] += v[i];
        }
      RadiusPixelType radius = m_RadiusMapImage->GetPixel( index );

      it_radius.Set( radius );

      ++it_vector;
      ++it_radius;
      }

    // Use Vector Image to find the closest vessel and add the tangent
    // direction
    typedef ImageRegionIterator<TangentImageType>   TangentIteratorType;
    TangentIteratorType it_tangent( m_TangentMapImage,
      m_TangentMapImage->GetLargestPossibleRegion() );

    it_vector.GoToBegin();
    it_tangent.GoToBegin();

    //Use Vector Image to add the radius values
    while( !it_vector.IsAtEnd() )
      {
      VectorPixelType v = it_vector.Value();
      typename DensityImageType::IndexType index = it_vector.GetIndex();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        index[i] += v[i];
        }
      TangentPixelType tangent = m_TangentMapImage->GetPixel( index );

      it_tangent.Set( tangent );

      ++it_vector;
      ++it_tangent;
      }

    //**** Invert the Distance Map image
    typename InverseIntensityImageFilter<DensityImageType>::Pointer
             InverseFilter;
    InverseFilter = InverseIntensityImageFilter<DensityImageType>::New();

    InverseFilter->SetInput( m_DensityMapImage );  //  inverse intensity
    InverseFilter->SetInverseMaximumIntensity( m_MaxDensityIntensity );
    InverseFilter->Update();

    m_DensityMapImage = InverseFilter->GetOutput();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr
      << "\n Error caught in TubeSpatialObjectToDensityImageFilter Class"
      << std::endl;
    std::cerr << e.GetDescription() <<std::endl;
    }
}

#endif // End !defined( __itktubeTubeSpatialObjectToDensityImageFilter_hxx )
