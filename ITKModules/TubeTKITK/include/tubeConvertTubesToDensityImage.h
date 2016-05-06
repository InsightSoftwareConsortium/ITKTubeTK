/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeConvertTubesToDensityImage_h
#define __tubeConvertTubesToDensityImage_h

// ITK includes
#include <itkGroupSpatialObject.h>
#include <itkMacro.h>
#include <itkProcessObject.h>


// TubeTK includes
#include "itktubeTubeSpatialObjectToDensityImageFilter.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ConvertTubesToDensityImage
 *
 *  \ingroup TubeTKITK
 */

template< class TDensityImageType, class TRadiusImageType,
          class TTangentImageType >
class ConvertTubesToDensityImage:
  public itk::ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef ConvertTubesToDensityImage                 Self;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  /** Typdefs */
  typedef TDensityImageType                              DensityImageType;
  typedef typename DensityImageType::PixelType           DensityPixelType;
  typedef typename DensityImageType::Pointer             DensityImagePointer;
  typedef typename DensityImageType::SizeType            SizeType;
  typedef typename DensityImageType::SpacingType         SpacingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertTubesToDensityImage, Object);

  /** Define the Dimension variable */
  itkStaticConstMacro( ImageDimension, unsigned int,
                       DensityImageType::ImageDimension );

  typedef TRadiusImageType                               RadiusImageType;
  typedef typename RadiusImageType::Pointer              RadiusImagePointer;

  typedef TTangentImageType                              TangentImageType;
  typedef typename TangentImageType::Pointer             TangentImagePointer;

  typedef itk::GroupSpatialObject< ImageDimension >      TubeGroupType;
  typedef typename TubeGroupType::Pointer                TubeGroupPointer;

  typedef itk::tube::TubeSpatialObjectToDensityImageFilter<
  DensityImageType, RadiusImageType, TangentImageType >  FilterType;

 // typedef typename FilterType::TubeGroupType              TubeGroupType;
 // typedef typename FilterType::TangentImageType           TangentImageType;
 // typedef typename FilterType::

  /** Set maximum density intensity value. Its a constant */
  tubeWrapSetMacro( MaxDensityIntensity, DensityPixelType, Filter );
  tubeWrapGetMacro( MaxDensityIntensity, DensityPixelType, Filter );

  /** Set the size of the output image volumes  */
  tubeWrapSetMacro( Size, SizeType, Filter );
  tubeWrapGetMacro( Size, SizeType, Filter );

  /** Set whether to use squared or actual distances in calculations. */
  tubeWrapSetMacro( UseSquareDistance, bool, Filter );
  tubeWrapGetMacro( UseSquareDistance, bool, Filter );

  /** Set the input tubes */
  tubeWrapSetMacro( InputTubeGroup, TubeGroupPointer, Filter );
  tubeWrapGetMacro( InputTubeGroup, TubeGroupPointer, Filter );

  /* Runs tubes to density image conversion */
  tubeWrapUpdateMacro(Filter);

  /* Get the generated density image */
  tubeWrapGetMacro( DensityMapImage, DensityImagePointer, Filter );

  /* Get the generated radius image */
  tubeWrapGetMacro( RadiusMapImage, RadiusImagePointer, Filter );

  /* Get the generated tangent map image */
  tubeWrapGetMacro( TangentMapImage, TangentImagePointer, Filter );

  /** Sets the element spacing */
  tubeWrapSetMacro( Spacing, SpacingType, Filter );
  tubeWrapGetMacro( Spacing, SpacingType, Filter );

protected:
  ConvertTubesToDensityImage( void );
  ~ConvertTubesToDensityImage() {}
  void PrintSelf( std::ostream & os, itk::Indent indent ) const;

private:
  /** itkConvertTubesToImageFilter parameters **/
  ConvertTubesToDensityImage(const Self &);
  void operator=(const Self &);

  typename FilterType::Pointer m_Filter;
};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertTubesToDensityImage.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToDensityImage_h )
