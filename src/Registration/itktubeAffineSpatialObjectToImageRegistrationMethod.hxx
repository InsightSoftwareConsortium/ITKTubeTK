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

#ifndef __itktubeAffineSpatialObjectToImageRegistrationMethod_hxx
#define __itktubeAffineSpatialObjectToImageRegistrationMethod_hxx


namespace itk
{

namespace tube
{

template <unsigned int ObjectDimension, class TImage>
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
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

template <unsigned int ObjectDimension, class TImage>
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::~AffineSpatialObjectToImageRegistrationMethod( void )
{
}

template <unsigned int ObjectDimension, class TImage>
void
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters( this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <unsigned int ObjectDimension, class TImage>
typename AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::TransformType
* AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <unsigned int ObjectDimension, class TImage>
const typename AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::TransformType
* AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <unsigned int ObjectDimension, class TImage>
typename AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>::AffineTransformPointer
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
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
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::SetInitialTransformParametersFromAffineTransform( const AffineTransformType * affine )
{
  this->SetInitialTransformFixedParameters( affine->GetFixedParameters() );
  this->SetInitialTransformParameters( affine->GetParameters() );
}

template <unsigned int ObjectDimension, class TImage>
void
AffineSpatialObjectToImageRegistrationMethod<ObjectDimension, TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

} // tube

} // itk

#endif
