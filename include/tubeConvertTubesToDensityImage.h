/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef tubeConvertTubesToDensityImage_h
#define tubeConvertTubesToDensityImage_h

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
 *  \ingroup TubeTK
 */

template <class TOutputPixel, unsigned int Dimension>
class ConvertTubesToDensityImage : public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = ConvertTubesToDensityImage;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;

  /** Typdefs */
  using DensityImageType = itk::Image<TOutputPixel, Dimension>;
  using DensityPixelType = typename DensityImageType::PixelType;
  using DensityImagePointer = typename DensityImageType::Pointer;
  using SizeType = typename DensityImageType::SizeType;
  using SpacingType = typename DensityImageType::SpacingType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkOverrideGetNameOfClassMacro(ConvertTubesToDensityImage);

  using RadiusImageType = itk::Image<TOutputPixel, Dimension>;
  using RadiusImagePointer = typename RadiusImageType::Pointer;
  using TangentPixelType = itk::Vector<TOutputPixel, Dimension>;
  using TangentImageType = itk::Image<TangentPixelType, Dimension>;
  using TangentImagePointer = typename TangentImageType::Pointer;

  using TubeGroupType = itk::GroupSpatialObject<Dimension>;
  using TubeGroupPointer = typename TubeGroupType::Pointer;

  using FilterType =
    itk::tube::TubeSpatialObjectToDensityImageFilter<DensityImageType, RadiusImageType, TangentImageType>;

  /** Set maximum density intensity value. Its a constant */
  tubeWrapSetMacro(MaxDensityIntensity, DensityPixelType, Filter);
  tubeWrapGetMacro(MaxDensityIntensity, DensityPixelType, Filter);

  /** Set the size of the output image volumes  */
  tubeWrapSetMacro(Size, SizeType, Filter);
  tubeWrapGetMacro(Size, SizeType, Filter);

  /** Set whether to use squared or actual distances in calculations. */
  tubeWrapSetMacro(UseSquaredDistance, bool, Filter);
  tubeWrapGetMacro(UseSquaredDistance, bool, Filter);

  /** Set the input tubes */
  tubeWrapSetMacro(InputTubeGroup, TubeGroupPointer, Filter);
  tubeWrapGetMacro(InputTubeGroup, TubeGroupPointer, Filter);

  /* Runs tubes to density image conversion */
  tubeWrapUpdateMacro(Filter);

  /* Get the generated density image */
  tubeWrapGetMacro(DensityMapImage, DensityImagePointer, Filter);

  /* Get the generated radius image */
  tubeWrapGetMacro(RadiusMapImage, RadiusImagePointer, Filter);

  /* Get the generated tangent map image */
  tubeWrapGetMacro(TangentMapImage, TangentImagePointer, Filter);

  /** Sets the element spacing */
  tubeWrapSetMacro(Spacing, SpacingType, Filter);
  tubeWrapGetMacro(Spacing, SpacingType, Filter);

protected:
  ConvertTubesToDensityImage(void);
  ~ConvertTubesToDensityImage() {}
  void
  PrintSelf(std::ostream & os, itk::Indent indent) const override;

private:
  /** itkConvertTubesToImageFilter parameters **/
  ConvertTubesToDensityImage(const Self &);
  void
  operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void
  SetInput(const DataObjectIdentifierType &, itk::DataObject *) override {};

  typename FilterType::Pointer m_Filter;
};

} // End namespace tube


#ifndef ITK_MANUAL_INSTANTIATION
#  include "tubeConvertTubesToDensityImage.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToDensityImage_h )
