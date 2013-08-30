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

#ifndef __itktubeSubSampleTubeTreeSpatialObjectFilterSerializer_hxx
#define __itktubeSubSampleTubeTreeSpatialObjectFilterSerializer_hxx

#include "itktubeSubSampleTubeTreeSpatialObjectFilterSerializer.h"

namespace itk
{

namespace tube
{

template< class TSubSampleTubeTreeSpatialObjectFilter >
SubSampleTubeTreeSpatialObjectFilterSerializer< TSubSampleTubeTreeSpatialObjectFilter >
::SubSampleTubeTreeSpatialObjectFilterSerializer( void )
{
  this->m_Sampling = new UnsignedIntegerValue;
  this->m_Sampling->SetValue( 1 );
  this->m_Parameters["Sampling"] = this->m_Sampling;
}


template< class TSubSampleTubeTreeSpatialObjectFilter >
SubSampleTubeTreeSpatialObjectFilterSerializer< TSubSampleTubeTreeSpatialObjectFilter >
::~SubSampleTubeTreeSpatialObjectFilterSerializer( void )
{
  delete m_Sampling;
}


template< class TSubSampleTubeTreeSpatialObjectFilter >
void
SubSampleTubeTreeSpatialObjectFilterSerializer< TSubSampleTubeTreeSpatialObjectFilter >
::Serialize( void )
{
  SubSampleTubeTreeSpatialObjectFilterType * filter =
    dynamic_cast< SubSampleTubeTreeSpatialObjectFilterType * >
      ( this->GetTargetObject() );
  if( filter != NULL )
    {
    this->m_Sampling->SetValue( filter->GetSampling() );
    this->m_Parameters["Sampling"] = this->m_Sampling;
    }

  Superclass::Serialize();
}


template< class TSubSampleTubeTreeSpatialObjectFilter >
void
SubSampleTubeTreeSpatialObjectFilterSerializer< TSubSampleTubeTreeSpatialObjectFilter >
::DeSerialize( void )
{
  Superclass::DeSerialize();

  SubSampleTubeTreeSpatialObjectFilterType * filter =
    dynamic_cast< SubSampleTubeTreeSpatialObjectFilterType * >
      ( this->GetTargetObject() );

  if( filter != NULL )
    {
    filter->SetSampling( this->m_Sampling->GetValue() );
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeSubSampleTubeTreeSpatialObjectFilterSerializer_hxx)
