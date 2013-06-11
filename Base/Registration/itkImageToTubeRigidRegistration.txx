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

#ifndef __itkImageToTubeRigidRegistration_txx
#define __itkImageToTubeRigidRegistration_txx

#include "itkImageToTubeRigidRegistration.h"

#include "itkSubSampleTubeTreeSpatialObjectFilter.h"

#include <itkLinearInterpolateImageFunction.h>
#include <itkNormalVariateGenerator.h>
#include <itkSpatialObjectDuplicator.h>

namespace itk
{

template < class TFixedImage, class TMovingSpatialObject, class TMovingTube >
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::ImageToTubeRigidRegistration( void )
{
  this->m_InitialTransformParameters = ParametersType( ParametersDimension );
  this->m_LastTransformParameters = ParametersType( ParametersDimension );

  this->m_InitialTransformParameters.Fill( 0.0f );
  this->m_LastTransformParameters.Fill( 0.0f );

  m_NumberOfIteration = 100; //by default
  m_LearningRate = 0.1;
  m_IsInitialized = false;

  m_InitialPosition.set_size( 6 );
  m_InitialPosition.Fill( 0.0 );
  m_ParametersScale.set_size( 6 );
  m_ParametersScale[0] = 30.;
  m_ParametersScale[1] = 30.;
  m_ParametersScale[2] = 30.;
  m_ParametersScale[3] = 1.;
  m_ParametersScale[4] = 1.;
  m_ParametersScale[5] = 1.;

  typename DefaultMetricType::Pointer metric = DefaultMetricType::New();
  this->SetMetric( metric );

  typedef LinearInterpolateImageFunction< FixedImageType > DefaultInterpolatorType;
  typename DefaultInterpolatorType::Pointer interpolator = DefaultInterpolatorType::New();
  this->SetInterpolator( interpolator );
}


template < class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::SetInitialPosition( const double position[6] )
{
  m_InitialPosition.set_size( 6 );
  for( unsigned int ii = 0; ii < 6; ++ii )
    {
    m_InitialPosition[ii]=position[ii];
    }
}


template < class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::SetParametersScale( const double scales[6] )
{
  m_ParametersScale.set_size( 6 );
  for( unsigned int ii = 0; ii < 6; ++ii )
    {
    m_ParametersScale[ii]=scales[ii];
    }
}


template < class TFixedImage, class TMovingSpatialObject, class TMovingTube >
void
ImageToTubeRigidRegistration< TFixedImage, TMovingSpatialObject, TMovingTube >
::Initialize() throw ( ExceptionObject )
{
  typename OptimizerType::Pointer optimizer = OptimizerType::New();

  typename TransformType::Pointer transform = TransformType::New();
  this->SetTransform( transform );

  optimizer->MaximizeOn();
  optimizer->SetScales( m_ParametersScale );

  // Gradient descent stuff
  optimizer->SetLearningRate( m_LearningRate );
  optimizer->SetNumberOfIterations( m_NumberOfIteration );


  /*  optimizer->SetMaximumIteration( m_NumberOfIteration );

  Statistics::NormalVariateGenerator::Pointer generator =
      Statistics::NormalVariateGenerator::New();
  generator->SetReferenceCount( 2 );
  generator->Initialize( std::time( NULL ) );

  optimizer->SetNormalVariateGenerator( generator );
  optimizer->Initialize( 40 );
  */

  this->SetOptimizer( optimizer );

  this->SetInitialTransformParameters( m_InitialPosition );

  try
    {
    // initialize the interconnects between components
    Superclass::Initialize();
    }
  catch( const ExceptionObject& )
    {
    this->m_LastTransformParameters = ParametersType( 1 );
    this->m_LastTransformParameters.Fill( 0.0f );
    // pass exception to caller
    throw;
    }
  this->GetOptimizer()->SetCostFunction( this->GetMetric() );

  m_IsInitialized = true;
}


template < class TFixedImage, class TMovingSpatialObject, class TMovingTube >
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

} // End namespace itk

#endif // End !defined(__itkImageToTubeRigidRegistration_txx)
