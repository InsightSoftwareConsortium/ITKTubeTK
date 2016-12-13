/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itktubeImageToTubeRigidRegistration_hxx
#define __itktubeImageToTubeRigidRegistration_hxx

#include "itktubeImageToTubeRigidRegistration.h"

#include "itktubeSubSampleTubeTreeSpatialObjectFilter.h"

#include <itkConjugateGradientOptimizer.h>
#include <itkGradientDescentOptimizer.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkSpatialObjectDuplicator.h>

namespace itk
{

namespace tube
{

template< class TFixedImage, class TMovingSpatialObject, class TMovingTube >
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::ImageToTubeRigidRegistration( void )
{
  this->m_InitialTransformParameters = ParametersType( ParametersDimension );
  this->m_LastTransformParameters = ParametersType( ParametersDimension );

  this->m_InitialTransformParameters.Fill( 0.0f );
  this->m_LastTransformParameters.Fill( 0.0f );

  m_IsInitialized = false;

  typename DefaultMetricType::Pointer metric = DefaultMetricType::New();
  this->SetMetric( metric );

  typedef LinearInterpolateImageFunction< FixedImageType > DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer interpolator = DefaultInterpolatorType::New();
  this->SetInterpolator( interpolator );

  typename Superclass::OptimizerType::ParametersType parameterScales( 6 );
  parameterScales[0] = 30.;
  parameterScales[1] = 30.;
  parameterScales[2] = 30.;
  parameterScales[3] = 1.;
  parameterScales[4] = 1.;
  parameterScales[5] = 1.;

  //typedef GradientDescentVariableStepOptimizer  DefaultOptimizerType;
  typedef GradientDescentOptimizer                DefaultOptimizerType;
  //typedef OnePlusOneEvolutionaryOptimizer       DefaultOptimizerType;
  typename DefaultOptimizerType::Pointer optimizer =
    DefaultOptimizerType::New();

  optimizer->MaximizeOn();
  optimizer->SetNumberOfIterations( 20 );
  optimizer->SetLearningRate( 0.1 );
  optimizer->SetScales( parameterScales );

  //optimizer->SetMaximumIteration( m_NumberOfIteration );

  //Statistics::NormalVariateGenerator::Pointer generator =
      //Statistics::NormalVariateGenerator::New();
  //generator->SetReferenceCount( 2 );
  //generator->Initialize( std::time( NULL ) );

  //optimizer->SetNormalVariateGenerator( generator );
  //optimizer->Initialize( 40 );

  this->SetOptimizer( optimizer );
}


template< class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::SetFeatureWeights( FeatureWeightsType & featureWeights )
{
  if( this->m_FeatureWeights.data_block() != featureWeights.data_block() ||
      this->m_FeatureWeights.GetSize() != featureWeights.GetSize() )
    {
    // m_FeatureWeights should treated as const.
    this->m_FeatureWeights.SetData( featureWeights.data_block(),
      featureWeights.GetSize(), false );
    this->Modified();
    }
}


template< class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::Initialize() throw ( ExceptionObject )
{
  typename TransformType::Pointer transform = TransformType::New();
  this->SetTransform( transform );

  try
    {
    // initialize the interconnects between components
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


template< class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
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

#endif // End !defined( __itktubeImageToTubeRigidRegistration_hxx )
