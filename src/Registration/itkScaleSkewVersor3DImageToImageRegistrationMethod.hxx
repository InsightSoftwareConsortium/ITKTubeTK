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

#ifndef __itkScaleSkewVersor3DImageToImageRegistrationMethod_txx
#define __itkScaleSkewVersor3DImageToImageRegistrationMethod_txx

#include "itkScaleSkewVersor3DImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::ScaleSkewVersor3DImageToImageRegistrationMethod( void )
{
  this->SetTransform( ScaleSkewVersor3DTransformType::New() );
  this->GetTypedTransform()->SetIdentity();

  this->SetInitialTransformParameters( this->GetTypedTransform()
                                       ->GetParameters() );
  this->SetInitialTransformFixedParameters( this->GetTypedTransform()
                                            ->GetFixedParameters() );

  typename Superclass::TransformParametersScalesType scales;
  scales.set_size( this->GetTypedTransform()->GetNumberOfParameters() );
  if( scales.size() != 12 )
    {
    std::cerr << "ERROR: number of parameters not standard for transform"
              << std::endl;
    std::cout << "   # = " << scales.size() << ", expecting 12" << std::endl;
    }
  unsigned int scaleNum = 0;
  // Versor
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    scales[scaleNum] = 1000;
    ++scaleNum;
    }
  // Offset
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    scales[scaleNum] = 1;
    ++scaleNum;
    }
  // Scale
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    scales[scaleNum] = 100;
    ++scaleNum;
    }
  // Skew
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    scales[scaleNum] = 1000;
    ++scaleNum;
    }
  this->SetTransformParametersScales( scales );
  this->SetTransformMethodEnum( Superclass::AFFINE_TRANSFORM );

  this->SetMaxIterations( 150 );
  this->SetNumberOfSamples( 150000 );
}

template <class TImage>
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::~ScaleSkewVersor3DImageToImageRegistrationMethod( void )
{
}

template <class TImage>
void
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters( this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <class TImage>
typename ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>::TransformType
* ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
const typename ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>::TransformType
* ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
typename ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>::AffineTransformPointer
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
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

template <class TImage>
void
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::SetInitialTransformParametersFromAffineTransform(
  const AffineTransformType * tfm )
{
  this->SetInitialTransformFixedParameters( tfm->GetFixedParameters() );
  this->SetInitialTransformParameters( tfm->GetParameters() );
}

template <class TImage>
void
ScaleSkewVersor3DImageToImageRegistrationMethod<TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif
