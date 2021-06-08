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

#ifndef __itkAffineSpatialObjectToImageRegistrationMethod_txx
#define __itkAffineSpatialObjectToImageRegistrationMethod_txx

#include "itkAffineSpatialObjectToImageRegistrationMethod.h"

namespace itk
{

template <class TSpatialObject, class TImage>
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::AffineSpatialObjectToImageRegistrationMethod( void )
{
  this->SetTransform( AffineTransformType::New() );
  this->GetTypedTransform()->SetIdentity();

  this->SetInitialTransformParameters( this->GetTypedTransform()
                                       ->GetParameters() );
  this->SetInitialTransformFixedParameters( this->GetTypedTransform()
                                            ->GetFixedParameters() );

  typename Superclass::TransformParametersScalesType scales;
  scales.set_size( this->GetTypedTransform()->GetNumberOfParameters() );
  if( scales.size() != ImageDimension * (ImageDimension + 1) )
    {
    std::cerr << "ERROR: number of parameters not standard for affine transform"
              << std::endl;
    }
  unsigned int scaleNum = 0;
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    for( unsigned int d2 = 0; d2 < ImageDimension; d2++ )
      {
      if( d1 == d2 )
        {
        scales[scaleNum] = 100;
        }
      else
        {
        scales[scaleNum] = 1000;
        }
      ++scaleNum;
      }
    }
  for( unsigned int d1 = 0; d1 < ImageDimension; d1++ )
    {
    scales[scaleNum] = 1;
    ++scaleNum;
    }
  this->SetTransformParametersScales( scales );

  this->SetTransformMethodEnum( Superclass::AFFINE_TRANSFORM );

  this->SetMaxIterations( 150 );
  this->SetNumberOfSamples( 150000 );
}

template <class TSpatialObject, class TImage>
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::~AffineSpatialObjectToImageRegistrationMethod( void )
{
}

template <class TSpatialObject, class TImage>
void
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters( this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <class TSpatialObject, class TImage>
typename AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::TransformType
* AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <class TSpatialObject, class TImage>
const typename AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::TransformType
* AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <class TSpatialObject, class TImage>
typename AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>::AffineTransformPointer
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
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
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::SetInitialTransformParametersFromAffineTransform( const AffineTransformType * affine )
{
  this->SetInitialTransformFixedParameters( affine->GetFixedParameters() );
  this->SetInitialTransformParameters( affine->GetParameters() );
}

template <class TSpatialObject, class TImage>
void
AffineSpatialObjectToImageRegistrationMethod<TSpatialObject, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif
