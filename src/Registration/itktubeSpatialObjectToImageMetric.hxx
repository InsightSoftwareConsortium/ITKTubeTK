/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef itkSpatialObjectToImageMetric_hxx
#define itkSpatialObjectToImageMetric_hxx

#include "itkSpatialObjectToImageMetric.h"

namespace itk
{
/** Constructor */
template <unsigned int ObjectDimension, typename TFixedImage>
SpatialObjectToImageMetric<TMovingSpatialObject, TFixedImage>::SpatialObjectToImageMetric()

{
  m_FixedImage = nullptr;          // has to be provided by the user.
  m_MovingSpatialObject = nullptr; // has to be provided by the user.
  m_Transform = nullptr;           // has to be provided by the user.

  m_SamplingRatio = 0.1;

  this->m_UseFixedImageRegionOfInterest = false;
  this->m_FixedImageRegionOfInterestPoint1.Fill(0);
  this->m_FixedImageRegionOfInterestPoint2.Fill(0);
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetFixedImage( const ImageType * fixedImage )
{
  if( this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;

    this->ProcessObject::SetNthInput(0, const_cast<ImageType *>( fixedImage ) );

    this->Modified();
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetMovingSpatialObject( const SpatialObjectType * movingSpatialObject )
{
  this->m_MovingSpatialObject = movingSpatialObject;

  this->ProcessObject::SetNthInput(1, const_cast<SpatialObjectType *>(
      m_MovingSpatialObject ) );

  this->Modified();
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetFixedImageRegionOfInterest( const PointType & point1,
  const PointType & point2 )
{
  m_FixedImageRegionOfInterestPoint1 = point1;
  m_FixedImageRegionOfInterestPoint2 = point2;
  m_UseFixedImageRegionOfInterest = true;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetFixedImageMaskObject( const MaskObjectType * maskObject )
{
  if( this->m_FixedImageMaskObject.GetPointer() != maskObject )
    {
    this->m_FixedImageMaskObject = maskObject;

    this->Modified();

    if( maskObject )
      {
      m_UseFixedImageMaskObject = true;
      }
    else
      {
      m_UseFixedImageMaskObject = false;
      }
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetMovingSpatialObjectMaskObject( const MaskObjectType * maskObject )
{
  if( this->m_MovingSpatialObjectMaskObject.GetPointer() != maskObject )
    {
    this->m_MovingSpatialObjectMaskObject = maskObject;

    this->Modified();

    if( maskObject )
      {
      m_UseMovingSpatialObjectMaskObject = true;
      }
    else
      {
      m_UseMovingSpatialObjectMaskObject = false;
      }
    }
}

/** Return the number of parameters required by the Transform */
template <typename TMovingSpatialObject, typename TFixedImage>
unsigned int
SpatialObjectToImageMetric<TMovingSpatialObject, TFixedImage>::GetNumberOfParameters() const
{
  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform is not present");
  }
  return m_Transform->GetNumberOfParameters();
}

/**
 * Initialize
 */

template <typename TMovingSpatialObject, typename TFixedImage>
void
SpatialObjectToImageMetric<TMovingSpatialObject, TFixedImage>::Initialize()
{
  if (!m_Transform)
  {
    itkExceptionMacro(<< "Transform is not present");
  }

  if (!m_MovingSpatialObject)
  {
    itkExceptionMacro(<< "MovingSpatialObject is not present");
  }

  if (!m_FixedImage)
  {
    itkExceptionMacro(<< "FixedImage is not present");
  }

  // If the image is provided by a source, update the source.
  if (m_FixedImage->GetSource())
  {
    m_FixedImage->GetSource()->Update();
  }

  // If there are any observers on the metric, call them to give the
  // user code a chance to set parameters on the metric
  this->InvokeEvent(InitializeEvent());
}

/** PrintSelf */
template <typename TMovingSpatialObject, typename TFixedImage>
void
SpatialObjectToImageMetric<TMovingSpatialObject, TFixedImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving Spatial Object: " << m_MovingSpatialObject.GetPointer() << std::endl;
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer() << std::endl;
  os << indent << "Last Transform parameters = " << m_LastTransformParameters << std::endl;

  os << indent << "Sampling Ratio:    " << m_SamplingRatio << std::endl;
}
} // end namespace itk

#endif
