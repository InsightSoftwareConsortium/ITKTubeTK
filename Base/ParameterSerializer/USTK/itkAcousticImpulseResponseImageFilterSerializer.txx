/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __itkAcousticImpulseResponseImageFilterSerializer_txx
#define __itkAcousticImpulseResponseImageFilterSerializer_txx

#include "itkAcousticImpulseResponseImageFilterSerializer.h"

#include "itkGradientMagnitudeImageFilter.h"
#include "itkGradientMagnitudeImageFilterSerializer.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilterSerializer.h"
#include "itkParameterSerializerValue.h"

namespace itk
{

namespace tube
{

template< class TAcousticImpulseResponseImageFilter >
AcousticImpulseResponseImageFilterSerializer< TAcousticImpulseResponseImageFilter >
::AcousticImpulseResponseImageFilterSerializer( void )
{
  this->m_AngleDependence = new DoubleValue;
  this->m_AngleDependence->SetValue( 1.0 );
  this->m_Parameters["AngleDependence"] = this->m_AngleDependence;

  typedef GradientMagnitudeImageFilter
    < typename AcousticImpulseResponseImageFilterType::OperatorImageType,
      typename AcousticImpulseResponseImageFilterType::OperatorImageType >
        DefaultGradientMagnitudeFilterType;
  typedef GradientMagnitudeImageFilterSerializer< DefaultGradientMagnitudeFilterType >
    DefaultGradientMagnitudeFilterSerializerType;
  this->m_GradientMagnitudeFilterSerializer = DefaultGradientMagnitudeFilterSerializerType::New();
  this->m_GradientMagnitudeFilter = new ParameterSerializerValue;
  this->m_GradientMagnitudeFilter->SetValue( m_GradientMagnitudeFilterSerializer.GetPointer() );
  this->m_Parameters["GradientMagnitudeFilter"] = this->m_GradientMagnitudeFilter;
}


template< class TAcousticImpulseResponseImageFilter >
AcousticImpulseResponseImageFilterSerializer< TAcousticImpulseResponseImageFilter >
::~AcousticImpulseResponseImageFilterSerializer( void )
{
  delete m_AngleDependence;
  delete m_GradientMagnitudeFilter;
}


template< class TAcousticImpulseResponseImageFilter >
void
AcousticImpulseResponseImageFilterSerializer< TAcousticImpulseResponseImageFilter >
::AssignGradientMagnitudeFilter( AcousticImpulseResponseImageFilterType * filter )
{
  typedef typename AcousticImpulseResponseImageFilterType::GradientMagnitudeFilterType
    GradientMagnitudeFilterBaseType;
  GradientMagnitudeFilterBaseType * gradientMagnitudeFilterBase =
    filter->GetGradientMagnitudeFilter();

  // Is it a GradientMagnitudeImageFilter?
  typedef GradientMagnitudeImageFilter
    < typename AcousticImpulseResponseImageFilterType::OperatorImageType,
      typename AcousticImpulseResponseImageFilterType::OperatorImageType >
        GradientMagnitudeImageFilterType;
  const GradientMagnitudeImageFilterType * gradientMagnitudeImageFilter =
    dynamic_cast< const GradientMagnitudeImageFilterType * >( gradientMagnitudeFilterBase );
  if( gradientMagnitudeImageFilter != NULL )
    {
    // Replace the serializer with one of the appropriate type.
    typedef GradientMagnitudeImageFilterSerializer< GradientMagnitudeImageFilterType >
      GradientMagnitudeImageFilterSerializerType;
    GradientMagnitudeImageFilterSerializerType * serializer =
      dynamic_cast< GradientMagnitudeImageFilterSerializerType * >
        ( m_GradientMagnitudeFilterSerializer.GetPointer() );
    if( serializer == NULL )
      {
      m_GradientMagnitudeFilterSerializer =
        GradientMagnitudeImageFilterSerializerType::New();
      m_GradientMagnitudeFilter->
        SetValue( m_GradientMagnitudeFilterSerializer.GetPointer() );
      }
    m_GradientMagnitudeFilterSerializer->SetTargetObject( gradientMagnitudeFilterBase );
    }
  else
    {
    // Is it a GradientMagnitudeRecursiveGaussianImageFilter?
    typedef GradientMagnitudeRecursiveGaussianImageFilter
      < typename AcousticImpulseResponseImageFilterType::OperatorImageType,
        typename AcousticImpulseResponseImageFilterType::OperatorImageType >
          GradientMagnitudeRecursiveGaussianImageFilterType;
    const GradientMagnitudeRecursiveGaussianImageFilterType * gradientMagnitudeRecursiveGaussianImageFilter =
      dynamic_cast< const GradientMagnitudeRecursiveGaussianImageFilterType * >( gradientMagnitudeFilterBase );
    if( gradientMagnitudeRecursiveGaussianImageFilter != NULL )
      {
      // Replace the serializer with one of the appropriate type.
      typedef GradientMagnitudeRecursiveGaussianImageFilterSerializer
        < GradientMagnitudeRecursiveGaussianImageFilterType >
        GradientMagnitudeRecursiveGaussianImageFilterSerializerType;
      GradientMagnitudeRecursiveGaussianImageFilterSerializerType * serializer =
        dynamic_cast< GradientMagnitudeRecursiveGaussianImageFilterSerializerType * >
          ( m_GradientMagnitudeFilterSerializer.GetPointer() );
      if( serializer == NULL )
        {
        m_GradientMagnitudeFilterSerializer =
          GradientMagnitudeRecursiveGaussianImageFilterSerializerType::New().GetPointer();
        m_GradientMagnitudeFilter
          ->SetValue( m_GradientMagnitudeFilterSerializer.GetPointer() );
        }
      m_GradientMagnitudeFilterSerializer->SetTargetObject( gradientMagnitudeFilterBase );
      }
    else
      {
      m_GradientMagnitudeFilterSerializer->SetTargetObject( gradientMagnitudeFilterBase );
      }
    }
}


template< class TAcousticImpulseResponseImageFilter >
void
AcousticImpulseResponseImageFilterSerializer< TAcousticImpulseResponseImageFilter >
::Serialize( void )
{
  AcousticImpulseResponseImageFilterType * filter =
    dynamic_cast< AcousticImpulseResponseImageFilterType * >
      ( this->GetTargetObject() );
  if( filter == NULL )
    {
    itkWarningMacro("AcousticImpulseResponseImageFilterSerializer target object not set");
    }
  else
    {
    this->m_AngleDependence->SetValue( filter->GetAngleDependence() );
    this->m_Parameters["AngleDependence"] = this->m_AngleDependence;

    this->AssignGradientMagnitudeFilter( filter );
    m_GradientMagnitudeFilterSerializer->Serialize();
    }

  Superclass::Serialize();
}


template< class TAcousticImpulseResponseImageFilter >
void
AcousticImpulseResponseImageFilterSerializer< TAcousticImpulseResponseImageFilter >
::DeSerialize( void )
{
  AcousticImpulseResponseImageFilterType * filter =
    dynamic_cast< AcousticImpulseResponseImageFilterType * >
      ( this->GetTargetObject() );

  if( filter != NULL )
    {
    this->AssignGradientMagnitudeFilter( filter );
    }

  Superclass::DeSerialize();

  if( filter != NULL )
    {
    filter->SetAngleDependence( this->m_AngleDependence->GetValue() );
    }
}

} // end namespace tube

} // end namespace itk

#endif
