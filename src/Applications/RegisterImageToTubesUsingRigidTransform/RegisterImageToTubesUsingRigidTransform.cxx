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

#include "itktubeImageToTubeRigidRegistration.h"
#include "itktubeTubeToTubeTransformFilter.h"
#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkTransformFileWriter.h>
#include <itkTimeProbesCollectorBase.h>


#include "RegisterImageToTubesUsingRigidTransformCLP.h"
#include "tubePreProcessRegistrationInputs.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Does not currently use TPixel
#define PARSE_ARGS_FLOAT_ONLY 1

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter(
    "RegisterImageToTubesUsingRigidTransform",
    CLPProcessInformation );
  progressReporter.Start();

  const unsigned int Dimension = 3;
  typedef double     FloatType;

  typedef itk::TubeSpatialObject< Dimension >            TubeType;
  typedef itk::GroupSpatialObject< Dimension >           TubeNetType;
  typedef itk::Image< FloatType, Dimension >             ImageType;

  typedef itk::tube::ImageToTubeRigidRegistration< ImageType, TubeNetType,
   TubeType >           RegistrationMethodType;

  typedef typename RegistrationMethodType::TransformType TransformType;

  typedef RegistrationMethodType::FeatureWeightsType     PointWeightsType;

  typename ImageType::Pointer currentImage;
  typename TubeNetType::Pointer tubeNet;
  PointWeightsType pointWeights;

  if( tube::PreProcessRegistrationInputs< Dimension, FloatType, TubeType,
    TubeNetType, ImageType, RegistrationMethodType >( argc, argv,
    timeCollector, currentImage, tubeNet, pointWeights )
    == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }


  timeCollector.Start( "Register image to tube" );

  typename RegistrationMethodType::Pointer registrationMethod =
    RegistrationMethodType::New();

  registrationMethod->SetFixedImage( currentImage );
  registrationMethod->SetMovingSpatialObject( tubeNet );
  registrationMethod->SetFeatureWeights( pointWeights );

  // Set Optimizer parameters.
  typename RegistrationMethodType::OptimizerType::Pointer optimizer =
    registrationMethod->GetOptimizer();
  itk::GradientDescentOptimizer * gradientDescentOptimizer = dynamic_cast<
    itk::GradientDescentOptimizer * >( optimizer.GetPointer() );
  if( gradientDescentOptimizer )
    {
    gradientDescentOptimizer->SetLearningRate( 0.1 );
    gradientDescentOptimizer->SetNumberOfIterations( 1000 );
    }

  // TODO: This is hard-coded now, which is sufficient since
  // ImageToTubeRigidMetric only uses a Euler3DTransform.  Will need to
  // adjust
  // to the transform parameters in the future at compile time.
  const unsigned int NumberOfParameters = 6;

  try
    {
    registrationMethod->Initialize();
    registrationMethod->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Performing registration: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  double progress = 0.9;
  progressReporter.Report( progress );

  TransformType* registrationTransform =
    dynamic_cast<TransformType *>( registrationMethod->GetTransform() );
  registrationTransform->SetParameters(
    registrationMethod->GetLastTransformParameters() );
  std::ostringstream parametersMessage;
  parametersMessage << "Transform Parameters: " <<
    registrationMethod->GetLastTransformParameters();
  tube::InformationMessage( parametersMessage.str() );
  parametersMessage.str( "" );
  parametersMessage << "Transform Center Of Rotation: " <<
    registrationTransform->GetFixedParameters();
  tube::InformationMessage( parametersMessage.str() );
  timeCollector.Stop( "Register image to tube" );

  timeCollector.Start( "Save data" );

  itk::TransformFileWriter::Pointer writer =
    itk::TransformFileWriter::New();
  writer->SetFileName( outputTransform.c_str() );
  writer->SetInput( registrationTransform );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing transform: Exception caught: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Save data" );
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
