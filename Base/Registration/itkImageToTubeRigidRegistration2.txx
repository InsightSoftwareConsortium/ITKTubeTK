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
#ifndef __itkImageToTubeRigidRegistration2_txx
#define __itkImageToTubeRigidRegistration2_txx

#include "itkImageToTubeRigidRegistration2.h"
#include <itkNormalVariateGenerator.h>

namespace itk
{

/** Constructor */
template <class TFixedImage, class TMovingTube>
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::ImageToTubeRigidRegistration2()
{
  this->m_InitialTransformParameters = ParametersType( ParametersDimension );
  this->m_LastTransformParameters = ParametersType( ParametersDimension );

  this->m_InitialTransformParameters.Fill( 0.0f );
  this->m_LastTransformParameters.Fill( 0.0f );

  m_NumberOfIteration = 100;
  m_LearningRate = 0.1;
  m_IsInitialized = false;

  m_InitialPosition.set_size( 6 );
  m_ParametersScale.set_size( 6 );
  m_InitialPosition.Fill( 0 );
  m_ParametersScale.Fill( 1 );

  m_MaskImage = 0;
  m_Extent = 3;
  m_Kappa = 1;
  m_Sampling = 30;
  m_Verbose = false;
}

/** Set the initial position */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::SetInitialPosition( double position[6] )
{
  m_InitialPosition.set_size( 6 );
  for( unsigned int i = 0; i < 6; ++i )
  {
    m_InitialPosition[i] = position[i];
  }
}

/** Set the parameters scales */
 template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::SetParametersScale( double scales[6] )
{
  m_ParametersScale.set_size( 6 );
  for( unsigned int i = 0; i < 6; ++i )
  {
    m_ParametersScale[i] = scales[i];
  }
}

/** Initialize by setting the interconnects
 *  between components. */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::Initialize() throw ( ExceptionObject )
{
  typename MetricType::Pointer metric = MetricType::New();

  metric->SetExtent( m_Extent );
  metric->SetKappa( m_Kappa );
  metric->SetSampling( m_Sampling );
  metric->SetVerbose( m_Verbose );

  this->SetMetric( metric );
  typename OptimizerType::Pointer optimizer = OptimizerType::New();
  typename TransformType::Pointer transform = TransformType::New();
  this->SetTransform( transform );
  typename InterpolatorType::Pointer interp = InterpolatorType::New();
  this->SetInterpolator( interp );

  optimizer->SetScales( m_ParametersScale );
  optimizer->MaximizeOn();
  optimizer->SetScales( m_ParametersScale );

  // Gradient descent parameters
  optimizer->SetLearningRate( m_LearningRate );
  optimizer->SetNumberOfIterations( m_NumberOfIteration );

  this->SetOptimizer( optimizer );
  this->SetInitialTransformParameters( m_InitialPosition );

  try
    {
    // initialize the interconnects between components
    Superclass::Initialize();
    }
  catch( ExceptionObject& err )
    {
    this->m_LastTransformParameters = ParametersType( 1 );
    this->m_LastTransformParameters.Fill( 0.0f );
    // pass exception to caller
    throw err;
    }
  this->GetOptimizer()->SetCostFunction( this->GetMetric() );

  m_IsInitialized = true;
}

/** Starts the Registration Process */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::StartRegistration( void )
{
  if( !m_IsInitialized )
    {
    this->Initialize();
    }

  try
    {
    this->SparseRegistration();
    this->GetOptimizer()->StartOptimization();
    }
  catch( ExceptionObject& err )
    {
    // An error has occurred in the optimization.
    // Update the parameters
    this->m_LastTransformParameters =
      this->GetOptimizer()->GetCurrentPosition();
    // Pass exception to caller
    throw err;
    }

  this->m_LastTransformParameters = this->GetOptimizer()->GetCurrentPosition();

  std::cout << "Solution: "
            << this->m_LastTransformParameters << std::endl;
}

/** Apply a sparse registration */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::SparseRegistration()
{
  if( !m_IsInitialized )
    {
    this->Initialize();
    }

  ParametersType parameters = this->GetLastTransformParameters();
  ParametersType params( MetricType::SpaceDimension );
  ParametersType optimalParameters( MetricType::SpaceDimension );

  double optimalValue = 0;
  for( float a = -0.1; a <= 0.1; a += 0.1 )
    {
    params[0] = parameters[0] + a;
    for( float b = -0.1; b <= 0.1; b += 0.1 )
      {
      params[1] = parameters[1] + b;
      for( float c = -0.1; c <= 0.1; c += 0.1 )
        {
        params[2] = parameters[2] + c;
        for( int x = -10; x <= 10; x += 10 )
          {
          params[3] = parameters[3] + x;
          for( int y = -10; y <= 10; y += 10 )
            {
            params[4] = parameters[4] + y;
            for( int z = -10; z <= 10;z += 10 )
              {
              params[5] = parameters[5] + z;
              double value = this->GetMetric()->GetValue( params );
              if( value > optimalValue )
                {
                optimalValue = value;
                optimalParameters[0] = params[0];
                optimalParameters[1] = params[1];
                optimalParameters[2] = params[2];
                optimalParameters[3] = params[3];
                optimalParameters[4] = params[4];
                optimalParameters[5] = params[5];
                }
              }
            }
          }
        }
      }
    }

  this->SetParametersScale( optimalParameters );
}

/** Get Center of Rotation */
template <class TFixedImage, class TMovingTube>
typename ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::MetricType::PointType
ImageToTubeRigidRegistration2<TFixedImage, TMovingTube>
::GetCenterOfRotation( void )
{
  typename MetricType::Pointer metric =
    dynamic_cast<MetricType*>( this->GetMetric() );

  return metric->GetCenterOfRotation();
}

} // end namespace itk


#endif
