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
#ifndef tubeCropTubes_h
#define tubeCropTubes_h

#include "itktubeCropTubesFilter.h"
#include "itkProcessObject.h"
#include "tubeWrappingMacros.h"

namespace tube
{
/** \class ComputeBinaryImageSimilarityMetrics
 *
 *  \ingroup TubeTK
 */

template <unsigned int VDimension>
class CropTubes : public itk::ProcessObject
{
public:
  /** Standard class type alias. */
  using Self = CropTubes;
  using Superclass = itk::ProcessObject;
  using Pointer = itk::SmartPointer<Self>;
  using ConstPointer = itk::SmartPointer<const Self>;
  using FilterType = itk::tube::CropTubesFilter<VDimension>;

  using TubeGroupType = typename FilterType::TubeGroupType;
  using TubeGroupPointer = typename TubeGroupType::Pointer;
  using ImageType = typename FilterType::ImageType;
  using PointType = typename FilterType::PointType;
  using VectorType = typename FilterType::VectorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ). */
  itkTypeMacro(CropTubes, ProcessObject);

  /* Set input tubes */
  tubeWrapSetMacro(Input, TubeGroupPointer, Filter);
  tubeWrapGetConstObjectMacro(Input, TubeGroupType, Filter);

  /** Set/Get bounding box corner */
  tubeWrapSetMacro(BoxPositionInWorldSpace, PointType, Filter);
  tubeWrapGetMacro(BoxPositionInWorldSpace, PointType, Filter);

  /** Set/Get bounding box size */
  tubeWrapSetMacro(BoxSizeInWorldSpace, VectorType, Filter);
  tubeWrapGetMacro(BoxSizeInWorldSpace, VectorType, Filter);

  /** Set/Get mask image */
  tubeWrapSetObjectMacro(MaskImage, ImageType, Filter);
  tubeWrapGetObjectMacro(MaskImage, ImageType, Filter);

  /** Set/Get Use mask image to crop tubes */
  tubeWrapSetMacro(UseMaskImage, bool, Filter);
  tubeWrapGetMacro(UseMaskImage, bool, Filter);

  /** Set/Get bool to crop tubes or not */
  tubeWrapSetMacro(CropTubes, bool, Filter);
  tubeWrapGetMacro(CropTubes, bool, Filter);

  /** Run Crop Tubes application */
  tubeWrapUpdateMacro(Filter);

  /* Get the crop tubes output */
  tubeWrapGetMacro(Output, TubeGroupPointer, Filter);

protected:
  CropTubes(void);
  ~CropTubes() {}

  void
  PrintSelf(std::ostream & os, itk::Indent indent) const override;

private:
  /** itktubeCropTubesFilter parameters **/
  CropTubes(const Self &);

  void
  operator=(const Self &);

  // To remove warning "was hidden [-Woverloaded-virtual]"
  void
  SetInput(const DataObjectIdentifierType &, itk::DataObject *) override {};

  typename FilterType::Pointer m_Filter;
};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#  include "tubeCropTubes.hxx"
#endif

#endif // End !defined( __tubeCropTubes_h )
