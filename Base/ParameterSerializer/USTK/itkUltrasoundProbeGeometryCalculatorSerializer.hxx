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

#ifndef __itkUltrasoundProbeGeometryCalculatorSerializer_hxx
#define __itkUltrasoundProbeGeometryCalculatorSerializer_hxx

#include "itkUltrasoundProbeGeometryCalculatorSerializer.h"

namespace itk
{

namespace tube
{

template< class TUltrasoundProbeGeometryCalculator >
UltrasoundProbeGeometryCalculatorSerializer< TUltrasoundProbeGeometryCalculator >
::UltrasoundProbeGeometryCalculatorSerializer( void )
{
  this->m_GeneralBeamDirection = new UnsignedIntegerValue;
  this->m_GeneralBeamDirection->SetValue( 1 );
  this->m_Parameters["GeneralBeamDirection"] = this->m_GeneralBeamDirection;
  this->m_BackgroundValue = new DoubleValue;
  this->m_BackgroundValue->SetValue( 0.0 );
  this->m_Parameters["BackgroundValue"] = this->m_BackgroundValue;
}


template< class TUltrasoundProbeGeometryCalculator >
UltrasoundProbeGeometryCalculatorSerializer< TUltrasoundProbeGeometryCalculator >
::~UltrasoundProbeGeometryCalculatorSerializer( void )
{
  delete m_GeneralBeamDirection;
  delete m_BackgroundValue;
}


template< class TUltrasoundProbeGeometryCalculator >
void
UltrasoundProbeGeometryCalculatorSerializer< TUltrasoundProbeGeometryCalculator >
::Serialize( void )
{
  UltrasoundProbeGeometryCalculatorType * calculator =
    dynamic_cast< UltrasoundProbeGeometryCalculatorType * >
      ( this->GetTargetObject() );
  if( calculator == NULL )
    {
    itkWarningMacro(
      << "UltrasoundProbeGeometryCalculatorSerializer target object not set" );
    }
  else
    {
    this->m_GeneralBeamDirection->SetValue( calculator->GetGeneralBeamDirection() );
    this->m_Parameters["GeneralBeamDirection"] = this->m_GeneralBeamDirection;

    this->m_BackgroundValue->SetValue(
      static_cast< double >( calculator->GetBackgroundValue() ) );
    this->m_Parameters["BackgroundValue"] = this->m_BackgroundValue;
    }

  Superclass::Serialize();
}


template< class TUltrasoundProbeGeometryCalculator >
void
UltrasoundProbeGeometryCalculatorSerializer< TUltrasoundProbeGeometryCalculator >
::DeSerialize( void )
{
  Superclass::DeSerialize();

  UltrasoundProbeGeometryCalculatorType * calculator =
    dynamic_cast< UltrasoundProbeGeometryCalculatorType * >
      ( this->GetTargetObject() );

  if( calculator != NULL )
    {
    calculator->SetGeneralBeamDirection( this->m_GeneralBeamDirection->GetValue() );
    calculator->SetBackgroundValue(
      static_cast< typename UltrasoundProbeGeometryCalculatorType::InputPixelType >(
        this->m_BackgroundValue->GetValue() ) );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itkUltrasoundProbeGeometryCalculatorSerializer_hxx)
