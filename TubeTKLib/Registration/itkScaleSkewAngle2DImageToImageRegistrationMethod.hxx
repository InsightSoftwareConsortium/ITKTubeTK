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

#ifndef __itkScaleSkewAngle2DImageToImageRegistrationMethod_txx
#define __itkScaleSkewAngle2DImageToImageRegistrationMethod_txx

#include "itkScaleSkewAngle2DImageToImageRegistrationMethod.h"

namespace itk
{

template <class TImage>
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::ScaleSkewAngle2DImageToImageRegistrationMethod( void )
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
  if( scales.size() != 7 )
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
  scales[scaleNum++] = 1000;
  this->SetTransformParametersScales( scales );

  this->SetTransformMethodEnum( Superclass::AFFINE_TRANSFORM );

  this->SetMaxIterations( 150 );
  this->SetNumberOfSamples( 150000 );
}

template <class TImage>
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::~ScaleSkewAngle2DImageToImageRegistrationMethod( void )
{
}

template <class TImage>
void
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::GenerateData( void )
{
  // Set the center of rotation
  this->GetTransform()->SetFixedParameters(
    this->GetInitialTransformFixedParameters() );

  Superclass::GenerateData();
}

template <class TImage>
typename ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>::TransformType
* ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void )
{
  return static_cast<TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
const typename ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>::TransformType
* ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::GetTypedTransform( void ) const
{
  return static_cast<const TransformType  *>( Superclass::GetTransform() );
}

template <class TImage>
typename ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>::AffineTransformPointer
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
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
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::SetInitialTransformParametersFromAffineTransform(
  const AffineTransformType * tfm )
{
  this->SetInitialTransformFixedParameters( tfm->GetFixedParameters() );
  this->SetInitialTransformParameters( tfm->GetParameters() );
}

template <class TImage>
void
ScaleSkewAngle2DImageToImageRegistrationMethod<TImage>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);
}

}

#endif
