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

#ifndef itktubeSpatialObjectToImageMetric_hxx
#define itktubeSpatialObjectToImageMetric_hxx


namespace itk
{

namespace tube
{

/** Constructor */
template <unsigned int ObjectDimension, typename TFixedImage>
SpatialObjectToImageMetric<ObjectDimension, TFixedImage>::SpatialObjectToImageMetric()

{
  m_FixedImage = nullptr;          // has to be provided by the user.

  m_MovingSpatialObject = nullptr; // has to be provided by the user.

  m_Transform = nullptr;           // has to be provided by the user.

  m_UseFixedImageMaskObject = false;
  m_FixedImageMaskObject = nullptr;
  
  m_UseMovingSpatialObjectMaskObject = false;
  m_MovingSpatialObjectMaskObject = nullptr;

  m_LastTransformParameters.fill(0);
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetFixedImage( const FixedImageType * fixedImage )
{
  if( this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;

    //this->ProcessObject::SetNthInput(0, const_cast<ImageType *>( fixedImage ) );

    this->Modified();
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetMovingSpatialObject( const MovingSpatialObjectType * movingSpatialObject )
{
  this->m_MovingSpatialObject = movingSpatialObject;

  //this->ProcessObject::SetNthInput(1, const_cast<SpatialObjectType *>(
      //m_MovingSpatialObject ) );

  this->Modified();
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageMetric<ObjectDimension, TImage>
::SetFixedImageMaskObject( const ImageMaskObjectType * maskObject )
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
::SetMovingSpatialObjectMaskObject(
  const SpatialObjectMaskObjectType * maskObject )
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
template <unsigned int ObjectDimension, typename TFixedImage>
unsigned int
SpatialObjectToImageMetric<ObjectDimension, TFixedImage>::GetNumberOfParameters() const
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

template <unsigned int ObjectDimension, typename TFixedImage>
void
SpatialObjectToImageMetric<ObjectDimension, TFixedImage>::Initialize()
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
template <unsigned int ObjectDimension, typename TFixedImage>
void
SpatialObjectToImageMetric<ObjectDimension, TFixedImage>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Moving Spatial Object: "
    << m_MovingSpatialObject.GetPointer() << std::endl;
  os << indent << "Fixed  Image: " << m_FixedImage.GetPointer() << std::endl;
  os << indent << "Transform:    " << m_Transform.GetPointer() << std::endl;
  os << indent << "Last Transform parameters = " << m_LastTransformParameters
    << std::endl;
  if(m_UseMovingSpatialObjectMaskObject)
    {
    os << indent << "Use Moving Spatial Object Mask Object: True" << std::endl;
    }
  else
    {
    os << indent << "Use Moving Spatial Object Mask Object: False" << std::endl;
    }
  os << indent << "Moving Spatial Object Mask Object: "
    << m_MovingSpatialObjectMaskObject.GetPointer() << std::endl;
  if(m_UseFixedImageMaskObject)
    {
    os << indent << "Use Fixed Image Mask Object: True" << std::endl;
    }
  else
    {
    os << indent << "Use Fixed Image Mask Object: False" << std::endl;
    }
  os << indent << "Fixed Image Mask Object: "
    << m_FixedImageMaskObject.GetPointer() << std::endl;
}

} // end namespace tube

} // end namespace itk

#endif
