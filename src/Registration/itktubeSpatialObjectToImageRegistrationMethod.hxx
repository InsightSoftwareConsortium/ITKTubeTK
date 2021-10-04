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

#ifndef __itktubeSpatialObjectToImageRegistrationMethod_txx
#define __itktubeSpatialObjectToImageRegistrationMethod_txx

#include "itktubeSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::SpatialObjectToImageRegistrationMethod( void )
{
  this->SetNumberOfRequiredOutputs( 1 ); // the transform

  this->m_Transform = 0;
  typename TransformOutputType::Pointer transformDecorator =
    static_cast<TransformOutputType *>
    ( this->MakeOutput(static_cast<DataObjectPointerArraySizeType>(0)).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );

  this->m_RegistrationNumberOfWorkUnits = this->GetNumberOfWorkUnits();
  this->SetNumberOfWorkUnits( this->m_RegistrationNumberOfWorkUnits );

  this->m_FixedImage = 0;
  this->m_MovingSpatialObject = 0;

  this->m_UseFixedImageMaskObject = false;
  this->m_FixedImageMaskObject = 0;

  this->m_UseMovingSpatialObjectMaskObject = false;
  this->m_MovingSpatialObjectMaskObject = 0;

  this->m_Observer = 0;
  this->m_ReportProgress = false;

}

template <unsigned int ObjectDimension, class TImage>
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::~SpatialObjectToImageRegistrationMethod( void )
{
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
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
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::SetMovingSpatialObject( const SpatialObjectType * movingSpatialObject )
{
  if( this->m_MovingSpatialObject != movingSpatialObject )
    {
    this->m_MovingSpatialObject = movingSpatialObject;

    //this->ProcessObject::SetNthInput(1, m_MovingSpatialObject );

    this->Modified();
    }
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
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
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::SetMovingSpatialObjectMaskObject( const SpatialObjectMaskObjectType *
  maskObject )
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

template <unsigned int ObjectDimension, class TImage>
const typename SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::TransformOutputType
* SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetOutput() const
  {
  return static_cast<const TransformOutputType *>(
    this->ProcessObject::GetOutput( 0 ) );
  }

template <unsigned int ObjectDimension, class TImage>
DataObject::Pointer
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::MakeOutput( DataObjectPointerArraySizeType idx )
{
  switch( idx )
    {
    case 0:
      return static_cast<DataObject *>(
               TransformOutputType::New().GetPointer() );
      break;
    default:
      itkExceptionMacro(
        "MakeOutput request for an output number larger than the expected number of outputs" );
      return nullptr;
    }
}

template <unsigned int ObjectDimension, class TImage>
ModifiedTimeType
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetMTime( void ) const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long m;

  if( m_Transform.IsNotNull() )
    {
    m = m_Transform->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_FixedImage.IsNotNull() )
    {
    m = m_FixedImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_FixedImageMaskObject.IsNotNull() )
    {
    m = m_FixedImageMaskObject->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_MovingSpatialObject.IsNotNull() )
    {
    m = m_MovingSpatialObject->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_MovingSpatialObjectMaskObject.IsNotNull() )
    {
    m = m_MovingSpatialObjectMaskObject->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  return mtime;
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::Initialize( void )
{
  this->GetMultiThreader()->SetNumberOfWorkUnits( m_RegistrationNumberOfWorkUnits );

  if( m_Transform.IsNull() )
    {
    itkExceptionMacro( << "Transform is not set" );
    }

  if( m_FixedImage.IsNull() )
    {
    itkExceptionMacro( << "Fixed image is not set" );
    }

  if( m_MovingSpatialObject.IsNull() )
    {
    itkExceptionMacro( << "Moving image is not set" );
    }

  TransformOutputType * transformOutput =
    static_cast<TransformOutputType *>( this->ProcessObject::GetOutput( 0 ) );

  transformOutput->Set( m_Transform.GetPointer() );

}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GenerateData( void )
{
  this->Update();
}

template <unsigned int ObjectDimension, class TImage>
void
SpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Number of threads = " << this->m_RegistrationNumberOfWorkUnits
     << std::endl;
  if( this->m_Transform.IsNotNull() )
    {
    os << indent << "Transform = " << this->m_Transform << std::endl;
    }
  else
    {
    os << indent << "Transform = 0" << std::endl;
    }

  if( this->m_Observer.IsNotNull() )
    {
    os << indent << "Observer = " << this->m_Observer << std::endl;
    }
  else
    {
    os << indent << "Observer = 0" << std::endl;
    }

  if( this->m_FixedImage.IsNotNull() )
    {
    os << indent << "Fixed Image = " << this->m_FixedImage << std::endl;
    }
  else
    {
    os << indent << "Fixed Image = 0" << std::endl;
    }

  if( this->m_MovingSpatialObject.IsNotNull() )
    {
    os << indent << "Moving Image = " << this->m_FixedImage << std::endl;
    }
  else
    {
    os << indent << "Moving Image = 0" << std::endl;
    }

  if( this->m_FixedImageMaskObject.IsNotNull() )
    {
    os << indent << "Fixed Image Mask Object = " << this->m_FixedImageMaskObject
       << std::endl;
    }
  else
    {
    os << indent << "Fixed image mask = 0" << std::endl;
    }

  if( this->m_MovingSpatialObjectMaskObject.IsNotNull() )
    {
    os << indent << "Moving Image Mask Object = " << this->m_MovingSpatialObjectMaskObject
       << std::endl;
    }
  else
    {
    os << indent << "Moving image mask = 0" << std::endl;
    }

  os << indent << "Report progress = " << this->m_ReportProgress << std::endl;

}

}; // tube

}; // itk

#endif
