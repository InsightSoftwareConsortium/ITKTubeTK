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

#ifndef __itkScaleSkewAngle2DSpatialObjectToImageRegistrationMethod_txx
#define __itkScaleSkewAngle2DSpatialObjectToImageRegistrationMethod_txx

#include "itkScaleSkewAngle2DSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

template <class TSpatialObject, class TImage>
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
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

  this->SetMaxIterations( 150 );
  this->SetNumberOfSamples( 150000 );
}

template <class TSpatialObject, class TImage>
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::~ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod( void )
{
}

template <class TSpatialObject, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters(
    this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <class TSpatialObject, class TImage>
typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::TransformType
* ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <class TSpatialObject, class TImage>
const typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::TransformType
* ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <class TSpatialObject, class TImage>
typename ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::AffineTransformPointer
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
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

template <class TSpatialObject, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::SetInitialTransformParametersFromAffineTransform(
  const AffineTransformType * tfm )
{
  this->SetInitialTransformFixedParameters( tfm->GetFixedParameters() );
  this->SetInitialTransformParameters( tfm->GetParameters() );
}

template <class TSpatialObject, class TImage>
void
ScaleSkewAngle2DSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif
