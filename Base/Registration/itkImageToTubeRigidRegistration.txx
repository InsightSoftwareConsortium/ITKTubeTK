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
#ifndef __itkImageToTubeRigidRegistration_txx
#define __itkImageToTubeRigidRegistration_txx

#include "itkImageToTubeRigidRegistration.h"
#include <itkNormalVariateGenerator.h>


namespace itk
{

/** Constructor */
template <class TFixedImage, class TMovingTube>
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::ImageToTubeRigidRegistration()
{
  this->m_InitialTransformParameters = ParametersType( ParametersDimension );
  this->m_LastTransformParameters = ParametersType( ParametersDimension );

  this->m_InitialTransformParameters.Fill( 0.0f );
  this->m_LastTransformParameters.Fill( 0.0f );

  m_NumberOfIteration = 100; //by default
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
  m_Verbose = true;
}


/** Set the initial position */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::SetInitialPosition( double position[6] )
{
  m_InitialPosition.set_size( 6 );
  for( unsigned int i=0;i<6;i++ )
  {
    m_InitialPosition[i]=position[i];
  }
}

/** Set the parameters scales */
 template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::SetParametersScale( double scales[6] )
{
  m_ParametersScale.set_size( 6 );
  for( unsigned int i=0;i<6;i++ )
  {
    m_ParametersScale[i]=scales[i];
  }

}


/** Initialize by setting the interconnects
 *  between components. */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::Initialize() throw ( ExceptionObject )
{
  typename MetricType::Pointer metric = MetricType::New();

  metric->SetExtent( m_Extent );
  metric->SetKappa( m_Kappa );
  metric->SetSampling( m_Sampling );
  metric->SetVerbose( m_Verbose );
  metric->SetMaskImage( m_MaskImage );

  this->SetMetric( metric );
  typename OptimizerType::Pointer optimizer = OptimizerType::New();

  typename TransformType::Pointer transform = TransformType::New();
  this->SetTransform( transform );
  typename InterpolatorType::Pointer interp = InterpolatorType::New();
  this->SetInterpolator( interp );

  /*
  ParametersType  parametersScale( ParametersDimension );
  parametersScale[0] = 30; //20
  parametersScale[1] = 30;
  parametersScale[2] = 30;
  parametersScale[3] = 1;
  parametersScale[4] = 1;
  parametersScale[5] = 1;
  */

  optimizer->SetScales( m_ParametersScale );

  optimizer->MaximizeOn();
  optimizer->SetScales( m_ParametersScale );

  // Gradient descent stuff
  optimizer->SetLearningRate( m_LearningRate );
  optimizer->SetNumberOfIterations( m_NumberOfIteration );


  /*  optimizer->SetMaximumIteration( m_NumberOfIteration );

  Statistics::NormalVariateGenerator::Pointer generator =
      Statistics::NormalVariateGenerator::New();
  generator->SetReferenceCount( 2 );
  generator->Initialize( time( NULL ) );

  optimizer->SetNormalVariateGenerator( generator );
  optimizer->Initialize( 40 );
  */

  this->SetOptimizer( optimizer );

  /* ParametersType parameters = ParametersType( ParametersDimension );
  unsigned int k = 0;

  // Initialize the 3 rotation angles
  for ( unsigned int i=0; i<TFixedImage::ImageDimension; i++ )
    {
    parameters[ k++ ] = 0;
    }

  m_Parameters[ k++ ]=-0.0658653;
  m_Parameters[ k++ ]=0.0841746;
  m_Parameters[ k++ ]=-0.0794313;

  // Initialize the 3 translation offsets
  for ( unsigned int i=0; i<TFixedImage::ImageDimension; i++ )
    {
    parameters[ k++ ] = 0;
    }
  parameters[ k++ ]=30;//49.0362; //30
  parameters[ k++ ]=-20;//-24.9046; //-20
  parameters[ k++ ]=20;//26.0874; //20

  parameters[ k++ ]= -71;
  parameters[ k++ ]= -23;
  parameters[ k++ ]= -25.8;
  */

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


/**  */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::CreateMatlabMetric( const char* filename )
{
  FILE* fic = fopen( filename, "wb" );

  ParametersType params( ParametersDimension );
  params.Fill( 0 );

  unsigned long c0=clock();
  // Test for XY translation
  for( float x = -10; x<=10; x+=1.0 )
    {
    for( float y = -10; y<=10; y+=1.0 )
      {
        params[3] = x;
        params[4] = y;
        typename MetricType::DerivativeType  derivatives( 6 );
        derivatives.Fill( 0 );
        this->GetMetric()->GetDerivative( params, derivatives );
        fprintf( fic, "%f %f %f %f %f %f %f %f %f %f %f %f %f\n",
          params[0], params[1], params[2], params[3], params[4], params[5],
          derivatives[0], derivatives[1], derivatives[2], derivatives[3],
          derivatives[4], derivatives[5], 0 );
        //this->GetMetric()->GetValue( params ) );
      }
    }

  std::cout << "Total time to run MonteCarlo simulation: "
    << ( ( float )clock()-( float )c0 )/( float )CLOCKS_PER_SEC << " secs."
    << std::endl;

  fclose( fic );

}


/** Starts the Registration Process */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::StartRegistration( void )
{
  double c0 = clock();

  if( !m_IsInitialized )
    {
    this->Initialize();
    }

  /*
  if( m_IterationCommand )
    {
    m_IterationCommand->SetMaximumValue( m_NumberOfIteration );
    m_IterationCommand->SetOptimizer( ( OptimizerType* )
    this->GetOptimizer() );
    }
  */

  try
    {
    // do the optimization
    this->GetOptimizer()->StartOptimization();
    }
  catch( ExceptionObject& err )
    {
    // An error has occurred in the optimization.
    // Update the parameters
    this->m_LastTransformParameters = this->GetOptimizer()
      ->GetCurrentPosition();
    // Pass exception to caller
    throw err;
    }

  std::cout << "The Solution is : ";
  this->m_LastTransformParameters = this->GetOptimizer()
    ->GetCurrentPosition();
  // give the result to the superclass
  std::cout << this->m_LastTransformParameters << std::endl;

  std::cout << "Total Time = " << ( clock()-c0 )/( double )CLOCKS_PER_SEC
    << std::endl;

}

/** Apply the sparse registration */
template <class TFixedImage, class TMovingTube>
void
ImageToTubeRigidRegistration<TFixedImage, TMovingTube>
::SparseRegistration( ParametersType & parameters )
{
  if( !m_IsInitialized )
    {
    this->Initialize();
    }

  //this->GetMetric()->SetTransform( m_Transformation );
  //this->GetMetric()->SetExtent( 3 );
  //this->GetMetric()->SubSampleMovingTube( 40 );
  static_cast<itk::ImageToTubeRigidMetric<TFixedImage, TMovingTube>*>(
    this->GetMetric() )->SparseRegistration( parameters );
}

} // end namespace itk


#endif
