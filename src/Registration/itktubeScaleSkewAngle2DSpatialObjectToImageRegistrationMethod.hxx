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

#ifndef __itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod_txx
#define __itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod_txx

#include "itktubeScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod( void )
{
  this->SetTransform( ScaleSkewAngle2DTransformType::New() );
  this->GetTypedTransform()->SetUseSingleScale( true );
  this->GetTypedTransform()->SetIdentity();

  this->SetInitialTransformParameters( this->GetTypedTransform()
                                       ->GetParameters() );
  this->SetInitialTransformFixedParameters( this->GetTypedTransform()
                                            ->GetFixedParameters() );

  typename Superclass::TransformParametersScalesType scales;
  scales.set_size( this->GetTypedTransform()->GetNumberOfParameters() );
  if( scales.size() != 6 )
    {
    std::cerr << "ERROR: number of parameters not standard for transform"
              << std::endl;
    }
  unsigned int scaleNum = 0;
  // Angle
  scales[scaleNum++] = 1000;
  // Offset
  scales[scaleNum++] = 1;
  scales[scaleNum++] = 1;
  // Scale
  scales[scaleNum++] = 100;
  scales[scaleNum++] = 100;
  // Skew
  scales[scaleNum++] = 1000;
  this->SetTransformParametersScales( scales );

  this->SetTransformMethodEnum( Superclass::AFFINE_TRANSFORM );
}

template <unsigned int ObjectDimension, class TImage>
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::~ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod( void )
{
}

template <unsigned int ObjectDimension, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters(
    this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <unsigned int ObjectDimension, class TImage>
typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::TransformType
* ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <unsigned int ObjectDimension, class TImage>
const typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::TransformType
* ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <unsigned int ObjectDimension, class TImage>
typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::AffineTransformPointer
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetAffineTransform( void ) const
{
  AffineTransformPointer trans = AffineTransformType::New();

  const TransformType * typedTransform = this->GetTypedTransform();

  trans->SetIdentity();
  trans->SetCenter( typedTransform->GetCenter() );
  trans->SetMatrix( typedTransform->GetMatrix() );
  trans->SetOffset( typedTransform->GetOffset() );

  return trans;
}

template <unsigned int ObjectDimension, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::SetInitialTransformParametersFromAffineTransform(
  const AffineTransformType * tfm )
{
  this->SetInitialTransformFixedParameters( tfm->GetFixedParameters() );
  this->SetInitialTransformParameters( tfm->GetParameters() );
}

template <unsigned int ObjectDimension, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

} // tube

} // itk

#endif
