/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeTubeSpatialObjectToDensityImage_hxx
#define __itktubeTubeSpatialObjectToDensityImage_hxx

#include "itktubeTubeSpatialObjectToDensityImage.h"

/** Constructor */
template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
TubeSpatialObjectToDensityImage< TDensityImageType, TRadiusImageType,
                                 TTangentImageType >
::TubeSpatialObjectToDensityImage( void )
{
  for(unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_Size[i] = 0;
    m_Spacing[i] = 1;
    }
  m_Max = 255;                   //NumericTraits<DensityPixelType>::max();
  m_UseSquareDistance = false;
}

/** Destructor */
template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
TubeSpatialObjectToDensityImage< TDensityImageType, TRadiusImageType,
                                 TTangentImageType>
::~TubeSpatialObjectToDensityImage( void )
{
}

template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
void
TubeSpatialObjectToDensityImage< TDensityImageType, TRadiusImageType,
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
TubeSpatialObjectToDensityImage< TDensityImageType, TRadiusImageType,
                                 TTangentImageType >
::Update( void )
{
  if(m_Size[0] == 0)
    {
    std::cerr << "Error, no size parameters given " << std::endl;
    return;
    }
  try
    {
    TubeGroupPointer tubes = this->GetTubes();

    tubes->SetBoundingBoxChildrenDepth( tubes->GetMaximumDepth() );
    tubes->ComputeBoundingBox();

    typedef TubeSpatialObjectToImageFilter<ImageDimension,DensityImageType>
                 FilterType;
    typename FilterType::Pointer tubefilter = FilterType::New();

    //Set record radius value for filter
    tubefilter->SetBuildRadiusImage( true );
    tubefilter->SetBuildTangentImage( true );
    tubefilter->SetUseRadius( true );
    tubefilter->SetInput( tubes );
    tubefilter->SetSize( m_Size );
    tubefilter->SetSpacing( m_Spacing );
    tubefilter->Update();

    typename DanielssonFilterType::Pointer
                danFilter = DanielssonFilterType::New();

    danFilter->SetInput( tubefilter->GetOutput() );

    danFilter->Update();

    VectorImagePointer  vectorImage = danFilter->GetVectorDistanceMap();
    m_RadiusImage   = tubefilter->GetRadiusImage();
    m_TangentImage  = tubefilter->GetTangentImage();
    m_DensityImage  = danFilter->GetDistanceMap();

    // ** If Requested: Square the Dan.Dis. image values to get the
    //      Squared distance image ** //
    if( m_UseSquareDistance )
      {
      typedef ImageRegionIterator<DensityImageType>   DistanceIteratorType;
      DistanceIteratorType  it_dis( m_DensityImage,
                                    m_DensityImage->GetLargestPossibleRegion());
      it_dis.GoToBegin();
      while( !it_dis.IsAtEnd() )
        {
        DensityPixelType number = it_dis.Get();
        DensityPixelType square = number*number;
        it_dis.Set( square );
        ++it_dis;
        }
      }

    typedef ImageRegionIterator<VectorImageType>   VectorIteratorType;
    typedef ImageRegionIterator<RadiusImageType>   RadiusIteratorType;

    VectorIteratorType it_vector( vectorImage,
                                  vectorImage->GetLargestPossibleRegion() );
    RadiusIteratorType it_radius( m_RadiusImage,
                                  m_RadiusImage->GetLargestPossibleRegion() );

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
      RadiusPixelType radius = m_RadiusImage->GetPixel( index );

      it_radius.Set( radius );

      ++it_vector;
      ++it_radius;
      }

    // Use Vector Image to find the closest vessel and add the tangent direction
    typedef ImageRegionIterator<TangentImageType>   TangentIteratorType;
    TangentIteratorType it_tangent( m_TangentImage,
                                    m_TangentImage->GetLargestPossibleRegion());

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
      TangentPixelType tangent = m_TangentImage->GetPixel( index );

      it_tangent.Set( tangent );

      ++it_vector;
      ++it_tangent;
      }

    //**** Invert the Distance Map image
    typename InverseIntensityImageFilter<DensityImageType>::Pointer
             InverseFilter;
    InverseFilter = InverseIntensityImageFilter<DensityImageType>::New();

    InverseFilter->SetInput( m_DensityImage );  //  inverse intensity
    InverseFilter->SetInverseMaximumIntensity( m_Max );
    InverseFilter->Update();

    m_DensityImage = InverseFilter->GetOutput();
    }
  catch( itk::ExceptionObject &e )
    {
    std::cerr << "\n Error caught in TubeSpatialObjectToDensityImage Class"
              << std::endl;
    std::cerr << e.GetDescription() <<std::endl;
    }
}

#endif // End !defined(__itktubeTubeSpatialObjectToDensityImage_hxx)
