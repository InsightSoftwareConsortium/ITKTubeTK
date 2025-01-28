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


#ifndef __itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod_txx
#define __itktubeScaleSkewVersor3DSpatialObjectToImageRegistrationMethod_txx


namespace itk
{

namespace tube
{

template <class TImage>
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<
  TImage>::ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod(void)
{
  this->SetTransform(ScaleSkewVersor3DTransformType::New());
  this->GetTypedTransform()->SetIdentity();

  this->SetInitialTransformParameters(this->GetTypedTransform()->GetParameters());
  this->SetInitialTransformFixedParameters(this->GetTypedTransform()->GetFixedParameters());

  typename Superclass::TransformParametersScalesType scales;
  scales.set_size(this->GetTypedTransform()->GetNumberOfParameters());
  if (scales.size() != 12)
  {
    std::cerr << "ERROR: number of parameters not standard for transform" << std::endl;
    std::cout << "   # = " << scales.size() << ", expecting 12" << std::endl;
  }
  unsigned int scaleNum = 0;
  // Versor
  for (unsigned int d1 = 0; d1 < ImageDimension; d1++)
  {
    scales[scaleNum] = 1000;
    ++scaleNum;
  }
  // Offset
  for (unsigned int d1 = 0; d1 < ImageDimension; d1++)
  {
    scales[scaleNum] = 1;
    ++scaleNum;
  }
  // Scale
  for (unsigned int d1 = 0; d1 < ImageDimension; d1++)
  {
    scales[scaleNum] = 100;
    ++scaleNum;
  }
  // Skew
  for (unsigned int d1 = 0; d1 < ImageDimension; d1++)
  {
    scales[scaleNum] = 1000;
    ++scaleNum;
  }
  this->SetTransformParametersScales(scales);
  this->SetTransformMethodEnum(Superclass::AFFINE_TRANSFORM);
}

template <class TImage>
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<
  TImage>::~ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod(void)
{}

template <class TImage>
void
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::GenerateData(void)
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters(this->GetInitialTransformFixedParameters());

  Superclass::GenerateData();
}

template <class TImage>
typename ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::TransformType *
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::GetTypedTransform(void)
{
  return static_cast<TransformType *>(Superclass::GetTransform());
}

template <class TImage>
const typename ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::TransformType *
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::GetTypedTransform(void) const
{
  return static_cast<const TransformType *>(Superclass::GetTransform());
}

template <class TImage>
typename ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::AffineTransformPointer
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::GetAffineTransform(void) const
{
  AffineTransformPointer trans = AffineTransformType::New();

  const TransformType * typedTransform = this->GetTypedTransform();

  trans->SetIdentity();
  trans->SetCenter(typedTransform->GetCenter());
  trans->SetMatrix(typedTransform->GetMatrix());
  trans->SetOffset(typedTransform->GetOffset());

  return trans;
}

template <class TImage>
void
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::SetInitialTransformParametersFromAffineTransform(
  const AffineTransformType * tfm)
{
  this->SetInitialTransformFixedParameters(tfm->GetFixedParameters());
  this->SetInitialTransformParameters(tfm->GetParameters());
}

template <class TImage>
void
ScaleSkewVersor3DSpatialObjectToImageRegistrationMethod<TImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
}

} // namespace tube

} // namespace itk

#endif
