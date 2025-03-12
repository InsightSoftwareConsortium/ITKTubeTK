/*=========================================================================

   Library:   TubeTK

   Copyright Kitware Inc.

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

#ifndef itktubeResampleTubesFilter_h
#define itktubeResampleTubesFilter_h

#include "itkSpatialObject.h"
#include "itkImage.h"
#include <itkDisplacementFieldTransform.h>
#include "itktubeSpatialObjectFilter.h"
#include "itktubePointBasedSpatialObjectTransformFilter.h"
#include "itkTransformBase.h"
namespace itk
{
namespace tube
{

/** \class ResampleTubesFilter
 * \brief resamples a given tube spatial object.
 *
 */
template <unsigned int ObjectDimension>
class ResampleTubesFilter : public SpatialObjectFilter<ObjectDimension>
{
public:
  /** Standard class type alias. */
  using SpatialObjectType = itk::SpatialObject<ObjectDimension>;

  using Self = ResampleTubesFilter<ObjectDimension>;
  using Superclass = SpatialObjectFilter<ObjectDimension>;
  using Pointer = SmartPointer<Self>;
  using ConstPointer = SmartPointer<const Self>;

  using PixelType = char;
  using ImageType = itk::Image<PixelType, ObjectDimension>;

  /** Typedefs for Displacement field tranform.    */
  using DisplacementFieldTransformType = itk::DisplacementFieldTransform<double, ObjectDimension>;
  typedef typename DisplacementFieldTransformType::DisplacementFieldType DisplacementFieldType;

  /** Typedefs for transform read from a file    */
  using BaseTransformType = itk::TransformBaseTemplate<double>;
  using BaseTransformListType = std::list<BaseTransformType::Pointer>;

  /** Run-time type information ( and related methods ).   */
  itkTypeMacro(ResampleTubesFilter, SpatialObjectFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set/Get match image */
  itkSetObjectMacro(MatchImage, ImageType);
  itkGetModifiableObjectMacro(MatchImage, ImageType);

  /** Set/Get sampling factor */
  itkSetMacro(SamplingFactor, int);
  itkGetMacro(SamplingFactor, int);

  /** Set/Get  use Inverse Transform */
  itkSetMacro(UseInverseTransform, bool);
  itkGetMacro(UseInverseTransform, bool);

  void
  SetDisplacementField(DisplacementFieldType * field);
  void
  SetReadTransformList(const BaseTransformListType * tList);

protected:
  ResampleTubesFilter(void);
  virtual ~ResampleTubesFilter(void);

  virtual void
  GenerateData(void) override;

  void
  PrintSelf(std::ostream & os, Indent indent) const override;

private:
  // purposely not implemented
  ResampleTubesFilter(const Self &);
  // purposely not implemented
  void
  operator=(const Self &);

  typename ImageType::Pointer             m_MatchImage;
  int                                     m_SamplingFactor;
  bool                                    m_UseInverseTransform;
  const BaseTransformListType *           m_ReadTransformList;
  typename DisplacementFieldType::Pointer m_DisplacementField;

  void
  ReadImageTransform(typename SpatialObjectType::TransformType::Pointer & outputTransform);
  typename SpatialObjectType::Pointer
  ApplyDisplacementFieldTransform(typename SpatialObjectType::TransformType::ConstPointer & outputTransform);
  typename SpatialObjectType::Pointer
  ApplyInputTransform(typename SpatialObjectType::TransformType::ConstPointer & outputTransform);
  typename SpatialObjectType::Pointer
  ApplyIdentityAffineTransform(typename SpatialObjectType::TransformType::ConstPointer & outputTransform);
}; // End class ResampleTubesFilter

} // End namespace tube
} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#  include "itktubeResampleTubesFilter.hxx"
#endif

#endif // End !defined( __itktubeResampleTubesFilter_h )
