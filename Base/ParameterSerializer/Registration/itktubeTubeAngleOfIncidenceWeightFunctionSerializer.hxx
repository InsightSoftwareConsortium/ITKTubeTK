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

#ifndef __itktubeTubeAngleOfIncidenceWeightFunctionSerializer_hxx
#define __itktubeTubeAngleOfIncidenceWeightFunctionSerializer_hxx

#include "itktubeTubeAngleOfIncidenceWeightFunctionSerializer.h"

namespace itk
{

namespace tube
{

template< class TTubeAngleOfIncidenceWeightFunction >
TubeAngleOfIncidenceWeightFunctionSerializer< TTubeAngleOfIncidenceWeightFunction >
::TubeAngleOfIncidenceWeightFunctionSerializer( void )
{
  this->m_FractionalImportance = new DoubleValue;
  this->m_FractionalImportance->SetValue( 0.5 );
  this->m_Parameters["FractionalImportance"] = this->m_FractionalImportance;

  this->m_AngleDependence = new DoubleValue;
  this->m_AngleDependence->SetValue( 1.0 );
  this->m_Parameters["AngleDependence"] = this->m_AngleDependence;

  this->m_UltrasoundProbeOrigin = new DoubleArrayValue;
  const DoubleArrayValue::DoubleArrayType probeOrigin( 3, 0.0 );
  this->m_UltrasoundProbeOrigin->SetValue( probeOrigin );
  this->m_Parameters["UltrasoundProbeOrigin"] = this->m_UltrasoundProbeOrigin;
}


template< class TTubeAngleOfIncidenceWeightFunction >
TubeAngleOfIncidenceWeightFunctionSerializer< TTubeAngleOfIncidenceWeightFunction >
::~TubeAngleOfIncidenceWeightFunctionSerializer( void )
{
  delete m_FractionalImportance;
  delete m_AngleDependence;
  delete m_UltrasoundProbeOrigin;
}


template< class TTubeAngleOfIncidenceWeightFunction >
void
TubeAngleOfIncidenceWeightFunctionSerializer< TTubeAngleOfIncidenceWeightFunction >
::Serialize( void )
{
  TubeAngleOfIncidenceWeightFunctionType * object =
    dynamic_cast< TubeAngleOfIncidenceWeightFunctionType * >
      ( this->GetTargetObject() );
  if( object != NULL )
    {
    this->m_FractionalImportance->SetValue( object->GetFractionalImportance() );
    this->m_Parameters["FractionalImportance"] = this->m_FractionalImportance;

    this->m_AngleDependence->SetValue( object->GetAngleDependence() );
    this->m_Parameters["AngleDependence"] = this->m_AngleDependence;

    const typename TubeAngleOfIncidenceWeightFunctionType::PointType &
      targetOrigin = object->GetUltrasoundProbeOrigin();
    DoubleArrayValue::DoubleArrayType & serializerOrigin =
      this->m_UltrasoundProbeOrigin->GetModifiableValue();
    serializerOrigin.resize( targetOrigin.Size() );
    for( size_t ii = 0; ii < serializerOrigin.size(); ++ii )
      {
      serializerOrigin[ii] = targetOrigin[ii];
      }
    this->m_Parameters["UltrasoundProbeOrigin"] = this->m_UltrasoundProbeOrigin;
    }

  Superclass::Serialize();
}


template< class TTubeAngleOfIncidenceWeightFunction >
void
TubeAngleOfIncidenceWeightFunctionSerializer< TTubeAngleOfIncidenceWeightFunction >
::DeSerialize( void )
{
  Superclass::DeSerialize();

  TubeAngleOfIncidenceWeightFunctionType * object =
    dynamic_cast< TubeAngleOfIncidenceWeightFunctionType * >
      ( this->GetTargetObject() );

  if( object != NULL )
    {
    object->SetFractionalImportance( this->m_FractionalImportance->GetValue() );

    object->SetAngleDependence( this->m_AngleDependence->GetValue() );

    const DoubleArrayValue::DoubleArrayType & serializerOrigin =
      this->m_UltrasoundProbeOrigin->GetValue();
    typename TubeAngleOfIncidenceWeightFunctionType::PointType targetOrigin;
    for( size_t ii = 0; ii < targetOrigin.Size(); ++ ii )
      {
      targetOrigin[ii] = serializerOrigin[ii];
      }
    object->SetUltrasoundProbeOrigin( targetOrigin );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeTubeAngleOfIncidenceWeightFunctionSerializer_hxx)
