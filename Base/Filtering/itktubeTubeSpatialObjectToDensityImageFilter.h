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

#ifndef __itktubeTubeSpatialObjectToDensityImageFilter_h
#define __itktubeTubeSpatialObjectToDensityImageFilter_h

#include "itktubeInverseIntensityImageFilter.h"
#include "itktubeTubeSpatialObjectToImageFilter.h"

#include <itkDanielssonDistanceMapImageFilter.h>
#include <itkGroupSpatialObject.h>
#include <itkImage.h>
#include <itkImageFileWriter.h>
#include <itkVesselTubeSpatialObject.h>

namespace itk
{

namespace tube
{

template< class TDensityImageType, class TRadiusImageType = Image< float, 3 >,
          class TTangentImageType = Image< Vector< float, 3 >, 3 > >
class TubeSpatialObjectToDensityImageFilter : public Object
{
public:

  typedef TubeSpatialObjectToDensityImageFilter   Self;
  typedef Object                                  Superclass;
  typedef SmartPointer< Self >                    Pointer;

  itkNewMacro( Self );
  itkTypeMacro( TubeSpatialObjectToDensityImageFilter, Object );

  /** Typdefs */
  typedef TDensityImageType                              DensityImageType;
  typedef typename DensityImageType::PixelType           DensityPixelType;
  typedef typename DensityImageType::Pointer             DensityImagePointer;

  typedef TRadiusImageType                               RadiusImageType;
  typedef typename RadiusImageType::PixelType            RadiusPixelType;
  typedef typename RadiusImageType::Pointer              RadiusImagePointer;

  typedef TTangentImageType                              TangentImageType;
  typedef typename TangentImageType::PixelType           TangentPixelType;
  typedef typename TangentImageType::Pointer             TangentImagePointer;

  /** Define the Dimension variable */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       DensityImageType::ImageDimension );

  typedef GroupSpatialObject<
    itkGetStaticConstMacro( ImageDimension ) >  TubeGroupType;
  typedef typename TubeGroupType::Pointer       TubeGroupPointer;

  typedef VesselTubeSpatialObject<
    itkGetStaticConstMacro( ImageDimension ) >  TubeType;

  typedef typename DensityImageType::OffsetType VectorPixelType;
  typedef Image<
    VectorPixelType,
    itkGetStaticConstMacro( ImageDimension ) >  VectorImageType;
  typedef typename VectorImageType::Pointer     VectorImagePointer;

  typedef typename DensityImageType::SizeType     SizeType;
  typedef typename DensityImageType::SpacingType  SpacingType;

  typedef TubeSpatialObjectToImageFilter<
    itkGetStaticConstMacro( ImageDimension ),
    DensityImageType > TubetoImageFilterType;

  typedef DanielssonDistanceMapImageFilter<
    DensityImageType, DensityImageType > DanielssonFilterType;

  /** Retrieve Density map created by inverted Danielsson Distance Map */
  itkSetMacro( DensityMapImage, DensityImagePointer );
  itkGetMacro( DensityMapImage, DensityImagePointer );
  itkSetMacro( RadiusMapImage, RadiusImagePointer );
  itkGetMacro( RadiusMapImage, RadiusImagePointer );
  itkSetMacro( TangentMapImage, TangentImagePointer );
  itkGetMacro( TangentMapImage, TangentImagePointer );

  /** Use square distance instead of linear distance */
  itkSetMacro( UseSquareDistance, bool );
  itkGetMacro( UseSquareDistance, bool );
  itkSetMacro( MaxDensityIntensity, DensityPixelType);
  itkGetMacro( MaxDensityIntensity, DensityPixelType);
  itkSetMacro( Size, SizeType );
  itkGetMacro( Size, SizeType );
  /** Sets the input tubes */
  itkSetMacro( InputTubeGroup, TubeGroupPointer );
  itkGetMacro( InputTubeGroup, TubeGroupPointer );

  /** Sets the element spacing */
  void SetSpacing( SpacingType );
  void Update( void );

protected:

  TubeSpatialObjectToDensityImageFilter( void );
  ~TubeSpatialObjectToDensityImageFilter( void );

private:

  TubeGroupPointer                  m_InputTubeGroup;
  DensityImagePointer               m_DensityMapImage;
  RadiusImagePointer                m_RadiusMapImage;
  TangentImagePointer               m_TangentMapImage;
  SizeType                          m_Size;
  SpacingType                       m_Spacing;

  /** Max value allowed for inverse intensity filter */
  DensityPixelType                  m_MaxDensityIntensity;
  bool                              m_UseSquareDistance;

}; // End class TubeSpatialObjectToDensityImageFilter

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTubeSpatialObjectToDensityImageFilter.hxx"
#endif

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeTubeSpatialObjectToDensityImageFilter_h)
