/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: ITKHeader.h,v $
  Language:  C++
  Date:      $Date: 2007-07-10 11:35:36 -0400 (Tue, 10 Jul 2007) $
  Version:   $Revision: 0 $

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

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
