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

#ifndef __itktubeResampleTubesFilter_hxx
#define __itktubeResampleTubesFilter_hxx

#include <iterator>


#include "itkMath.h"

#include "tubeMacro.h"

#include "itkTubeSpatialObject.h"
#include "itktubeSubSampleSpatialObjectFilter.h"

namespace itk
{
namespace tube
{

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
ResampleTubesFilter<ObjectDimension>::ResampleTubesFilter(void)
{
  m_UseInverseTransform = false;
  m_MatchImage = nullptr;
  m_DisplacementField = nullptr;
  m_ReadTransformList = nullptr;
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
ResampleTubesFilter<ObjectDimension>::~ResampleTubesFilter(void)
{}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
void
ResampleTubesFilter<ObjectDimension>::SetDisplacementField(DisplacementFieldType * field)
{
  m_DisplacementField = field;
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
void
ResampleTubesFilter<ObjectDimension>::SetReadTransformList(const BaseTransformListType * tList)
{
  m_ReadTransformList = tList;
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
void
ResampleTubesFilter<ObjectDimension>::ReadImageTransform(
  typename SpatialObjectType::TransformType::Pointer & outputTransform)
{
  typename ImageType::PointType     origin = m_MatchImage->GetOrigin();
  typename ImageType::DirectionType directions = m_MatchImage->GetDirection();

  outputTransform = SpatialObjectType::TransformType::New();
  outputTransform->SetIdentity();
  itk::Vector<double, ObjectDimension> offset;
  for (unsigned int i = 0; i < ObjectDimension; ++i)
  {
    offset[i] = origin[i];
  }
  outputTransform->SetMatrix(directions);
  outputTransform->SetOffset(offset);
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
typename itk::SpatialObject<ObjectDimension>::Pointer
ResampleTubesFilter<ObjectDimension>::ApplyDisplacementFieldTransform(
  typename SpatialObjectType::TransformType::ConstPointer & outputTransform)
{
  typename SpatialObjectType::ConstPointer inSO = this->GetInput();

  /** Typedefs for Displacement field tranform filter.    */
  typedef itk::tube::PointBasedSpatialObjectTransformFilter<DisplacementFieldTransformType, ObjectDimension>
    DisplacementFieldTransformFilterType;

  // Create new transform
  typename DisplacementFieldTransformType::Pointer transform = DisplacementFieldTransformType::New();
  transform->SetDisplacementField(m_DisplacementField);

  // Create the filter and apply
  typename DisplacementFieldTransformFilterType::Pointer filter = DisplacementFieldTransformFilterType::New();
  filter->SetInput(inSO);
  filter->SetTransform(transform);
  filter->SetOutputObjectToParentTransform(outputTransform.GetPointer());
  filter->Update();

  this->GraftOutput(filter->GetOutput());

  return this->GetOutput();
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
typename itk::SpatialObject<ObjectDimension>::Pointer
ResampleTubesFilter<ObjectDimension>::ApplyInputTransform(
  typename SpatialObjectType::TransformType::ConstPointer & outputTransform)
{
  typename SpatialObjectType::ConstPointer inSO = this->GetInput();

  /** Typedefs for transform read from a file    */
  typedef itk::MatrixOffsetTransformBase<double, ObjectDimension, ObjectDimension> MatrixOffsetTransformType;
  typedef itk::tube::PointBasedSpatialObjectTransformFilter<MatrixOffsetTransformType, ObjectDimension>
    MatrixOffsetTransformFilterType;

  BaseTransformListType::const_iterator tListIt;
  tListIt = m_ReadTransformList->begin();
  while (tListIt != m_ReadTransformList->end())
  {
    typename MatrixOffsetTransformType::Pointer transform =
      dynamic_cast<MatrixOffsetTransformType *>((*tListIt).GetPointer());
    typename MatrixOffsetTransformFilterType::Pointer filter = MatrixOffsetTransformFilterType::New();
    if (m_UseInverseTransform)
    {
      typename MatrixOffsetTransformType::InverseTransformBaseType::Pointer ivT = transform->GetInverseTransform();
      transform = (MatrixOffsetTransformType *)ivT.GetPointer();
    }

    filter->SetInput(inSO);
    filter->SetTransform(transform);
    filter->SetOutputObjectToParentTransform(outputTransform);
    filter->Update();

    inSO = filter->GetOutput();
    this->GraftOutput(filter->GetOutput());
    ++tListIt;
  }

  return this->GetOutput();
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
typename itk::SpatialObject<ObjectDimension>::Pointer
ResampleTubesFilter<ObjectDimension>::ApplyIdentityAffineTransform(
  typename SpatialObjectType::TransformType::ConstPointer & outputTransform)
{
  typename SpatialObjectType::ConstPointer inSO = this->GetInput();

  /** Typedefs for Affine Transform */
  typedef itk::AffineTransform<double, ObjectDimension> AffineTransformType;
  typedef itk::tube::PointBasedSpatialObjectTransformFilter<AffineTransformType, ObjectDimension>
    AffineTransformFilterType;

  typename AffineTransformType::Pointer identityAffineTransform = AffineTransformType::New();
  identityAffineTransform->SetIdentity();

  typename AffineTransformFilterType::Pointer filter = AffineTransformFilterType::New();
  filter->SetInput(inSO);
  filter->SetTransform(identityAffineTransform);
  filter->SetOutputObjectToParentTransform(outputTransform);
  filter->Update();

  this->GraftOutput(filter->GetOutput());

  return this->GetOutput();
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
void
ResampleTubesFilter<ObjectDimension>::GenerateData(void)
{
  typename SpatialObjectType::ConstPointer inSO = this->GetInput();

  typename SpatialObjectType::Pointer outSO = this->GetOutput();

  typename SpatialObjectType::TransformType::ConstPointer outputTransformConst;
  if (m_MatchImage)
  {
    typename SpatialObjectType::TransformType::Pointer outputTransform;
    this->ReadImageTransform(outputTransform);
    outputTransformConst = outputTransform.GetPointer();
  }
  else
  {
    outputTransformConst = inSO->GetObjectToWorldTransform();
  }

  bool outSOUpdated = false;
  if (m_DisplacementField)
  {
    outSO = this->ApplyDisplacementFieldTransform(outputTransformConst);
    outSOUpdated = true;
  }
  else if (m_ReadTransformList)
  {
    outSO = this->ApplyInputTransform(outputTransformConst);
    outSOUpdated = true;
  }
  else if (m_MatchImage)
  {
    outSO = this->ApplyIdentityAffineTransform(outputTransformConst);
    outSOUpdated = true;
  }

  if (m_SamplingFactor != 1)
  {
    /** Typedefs for Sub samppling filter     */
    typedef itk::tube::SubSampleSpatialObjectFilter<ObjectDimension> SubSampleFilterType;

    typename SubSampleFilterType::Pointer subSampleFilter = SubSampleFilterType::New();
    if (outSOUpdated)
    {
      subSampleFilter->SetInput(outSO);
    }
    else
    {
      subSampleFilter->SetInput(inSO);
      subSampleFilter->GraftOutput(outSO);
    }
    subSampleFilter->SetSampling(m_SamplingFactor);

    try
    {
      subSampleFilter->Update();
    }
    catch (const std::exception & e)
    {
      std::cout << e.what();
      return;
    }

    this->GraftOutput(subSampleFilter->GetOutput());
    outSO = subSampleFilter->GetOutput();
  }
}

//--------------------------------------------------------------------------
template <unsigned int ObjectDimension>
void
ResampleTubesFilter<ObjectDimension>::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);
}

} // End namespace tube
} // End namespace itk

#endif // End !defined( __itktubeResampleTubesFilter_hxx )
