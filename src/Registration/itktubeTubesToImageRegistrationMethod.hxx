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

#ifndef __itktubeTubesToImageRegistrationMethod_hxx
#define __itktubeTubesToImageRegistrationMethod_hxx

#include "itktubeTubesToImageRegistrationMethod.h"

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"

#include <itkGradientDescentOptimizer.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>

namespace itk
{

namespace tube
{

template< class TMovingSpatialObject, class TFixedImage>
TubesToImageRegistrationMethod< TMovingSpatialObject, TFixedImage >
::TubesToImageRegistrationMethod( void )
{
  this->m_InitialTransformParameters = ParametersType( ParametersDimension );
  this->m_LastTransformParameters = ParametersType( ParametersDimension );

  this->m_InitialTransformParameters.Fill( 0.0f );
  this->m_LastTransformParameters.Fill( 0.0f );

  m_IsInitialized = false;
}

template< class TMovingSpatialObject, class TFixedImage>
void
TubesToImageRegistrationMethod< TMovingSpatialObject, TFixedImage >
::Initialize()
{
  typename TransformType::Pointer transform = TransformType::New();
  this->SetTransform( transform );

  try
    {
    Superclass::Initialize();
    DefaultMetricType * defaultMetric
      = dynamic_cast< DefaultMetricType * >( this->GetMetric() );
    if( defaultMetric != NULL )
      {
      defaultMetric->SetFeatureWeights( this->m_FeatureWeights );
      }
    }
  catch( const ExceptionObject& )
    {
    this->m_LastTransformParameters = ParametersType( 1 );
    this->m_LastTransformParameters.Fill( 0.0f );
    // pass exception to caller
    throw;
    }

  m_IsInitialized = true;
}


template< class TMovingSpatialObject, class TFixedImage>
void
TubesToImageRegistrationMethod< TMovingSpatialObject, TFixedImage >
::StartRegistration( void )
{
  if( !m_IsInitialized )
    {
    this->Initialize();
    }

  try
    {
    // do the optimization
    this->GetOptimizer()->StartOptimization();
    }
  catch( const ExceptionObject& err )
    {
    // An error has occurred in the optimization.
    // Update the parameters
    this->m_LastTransformParameters = this->GetOptimizer()
      ->GetCurrentPosition();
    // Pass exception to caller
    throw;
    }

  // give the result to the superclass
  this->m_LastTransformParameters = this->GetOptimizer()
    ->GetCurrentPosition();
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubesToImageRegistrationMethod_hxx )
