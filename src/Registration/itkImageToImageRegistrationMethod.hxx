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

#ifndef __ImageToImageRegistrationMethod_txx
#define __ImageToImageRegistrationMethod_txx

#include "itkImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
ImageToImageRegistrationMethod<TImage>
::ImageToImageRegistrationMethod( void )
{
  this->SetNumberOfRequiredOutputs( 1 ); // the transform

  this->m_Transform = 0;
  typename TransformOutputType::Pointer transformDecorator =
    static_cast<TransformOutputType *>
    ( this->MakeOutput(static_cast<DataObjectPointerArraySizeType>(0)).GetPointer() );

  this->ProcessObject::SetNthOutput( 0, transformDecorator.GetPointer() );

  this->m_RegistrationNumberOfWorkUnits = this->GetNumberOfWorkUnits();
  this->GetMultiThreader()->SetNumberOfWorkUnits( this->m_RegistrationNumberOfWorkUnits );

  this->m_FixedImage = 0;
  this->m_MovingImage = 0;
  this->m_UseFixedImageMaskObject = false;
  this->m_FixedImageMaskObject = 0;
  this->m_UseMovingImageMaskObject = false;
  this->m_MovingImageMaskObject = 0;
  this->m_Observer = 0;
  this->m_ReportProgress = false;

  this->m_UseRegionOfInterest = false;
  this->m_RegionOfInterestPoint1.Fill(0);
  this->m_RegionOfInterestPoint2.Fill(0);

}

template <class TImage>
ImageToImageRegistrationMethod<TImage>
::~ImageToImageRegistrationMethod( void )
{
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
::SetFixedImage( const ImageType * fixedImage )
{
  if( this->m_FixedImage.GetPointer() != fixedImage )
    {
    this->m_FixedImage = fixedImage;

    this->ProcessObject::SetNthInput(0, const_cast<ImageType *>( fixedImage ) );
    this->Modified();
    }
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
::SetMovingImage( const ImageType * movingImage )
{
  if( this->m_MovingImage.GetPointer() != movingImage )
    {
    this->m_MovingImage = movingImage;

    this->ProcessObject::SetNthInput(1, const_cast<ImageType *>( movingImage ) );
    this->Modified();
    }
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
::SetRegionOfInterest( const PointType & point1, const PointType & point2 )
{
  m_RegionOfInterestPoint1 = point1;
  m_RegionOfInterestPoint2 = point2;
  m_UseRegionOfInterest = true;
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
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

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
::SetMovingImageMaskObject( const MaskObjectType * maskObject )
{
  if( this->m_MovingImageMaskObject.GetPointer() != maskObject )
    {
    this->m_MovingImageMaskObject = maskObject;

    this->Modified();

    if( maskObject )
      {
      m_UseMovingImageMaskObject = true;
      }
    else
      {
      m_UseMovingImageMaskObject = false;
      }
    }
}

template <class TImage>
const typename ImageToImageRegistrationMethod<TImage>::TransformOutputType
* ImageToImageRegistrationMethod<TImage>
::GetOutput() const
  {
  return static_cast<const TransformOutputType *>( this->ProcessObject::GetOutput( 0 ) );
  }

template <class TImage>
DataObject::Pointer
ImageToImageRegistrationMethod<TImage>
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

template <class TImage>
ModifiedTimeType
ImageToImageRegistrationMethod<TImage>
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

  if( m_MovingImage.IsNotNull() )
    {
    m = m_MovingImage->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  if( m_MovingImageMaskObject.IsNotNull() )
    {
    m = m_MovingImageMaskObject->GetMTime();
    mtime = (m > mtime ? m : mtime);
    }

  return mtime;
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
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

  if( m_MovingImage.IsNull() )
    {
    itkExceptionMacro( << "Moving image is not set" );
    }

  TransformOutputType * transformOutput =
    static_cast<TransformOutputType *>( this->ProcessObject::GetOutput( 0 ) );

  transformOutput->Set( m_Transform.GetPointer() );

}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
::GenerateData( void )
{
  this->Update();
}

template <class TImage>
void
ImageToImageRegistrationMethod<TImage>
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

  if( this->m_MovingImage.IsNotNull() )
    {
    os << indent << "Moving Image = " << this->m_FixedImage << std::endl;
    }
  else
    {
    os << indent << "Moving Image = 0" << std::endl;
    }

  os << indent << "Use region of interest = " << m_UseRegionOfInterest
     << std::endl;
  os << indent << "Region of interest point1 = " << m_RegionOfInterestPoint1
     << std::endl;
  os << indent << "Region of interest point2 = " << m_RegionOfInterestPoint2
     << std::endl;

  if( this->m_FixedImageMaskObject.IsNotNull() )
    {
    os << indent << "Fixed Image Mask Object = " << this->m_FixedImageMaskObject
       << std::endl;
    }
  else
    {
    os << indent << "Fixed image mask = 0" << std::endl;
    }

  if( this->m_MovingImageMaskObject.IsNotNull() )
    {
    os << indent << "Moving Image Mask Object = " << this->m_MovingImageMaskObject
       << std::endl;
    }
  else
    {
    os << indent << "Moving image mask = 0" << std::endl;
    }

  os << indent << "Report progress = " << this->m_ReportProgress << std::endl;

}

};

#endif
