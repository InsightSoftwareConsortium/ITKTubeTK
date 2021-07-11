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

#include "itktubeSpatialObjectToImageRegistration.h"

#include <itkImageFileReader.h>
#include <itkMemoryProbesCollectorBase.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate( void ) {}
  ~CommandIterationUpdate( void )
    {
    if( measuresFileStream.is_open() )
      {
      measuresFileStream.close();
      }
    }

  std::ofstream measuresFileStream;

public:
  typedef itk::GradientDescentOptimizer OptimizerType;
  typedef const OptimizerType *         OptimizerPointer;

  void Execute( itk::Object *caller, const itk::EventObject & event ) override
    {
    Execute( ( const itk::Object * )caller, event );
    }

  void Execute( const itk::Object * object, const itk::EventObject & event )
    override
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( !( itk::IterationEvent().CheckEvent( &event ) ) )
      {
      return;
      }

    if( measuresFileStream.is_open() )
      {
      measuresFileStream << optimizer->GetCurrentIteration() << ";";
      measuresFileStream << optimizer->GetValue() << std::endl;
      }
    }

  void SetFileName( char* fileName )
    {
    measuresFileStream.open( fileName, std::ios::trunc );
    }

}; // End class CommandIterationUpdate

int itktubeSpatialObjectToImageRegistrationPerformanceTest( int argc, char * argv[] )
{
  if( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image " << "Input_Vessel " << "Output_Results"
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::TubeSpatialObject< 3 >                    TubeType;
  typedef itk::GroupSpatialObject<3>                     TubeNetType;
  typedef itk::SpatialObjectReader<3>                    TubeNetReaderType;
  typedef itk::Image<double, 3>                          ImageType;
  typedef itk::ImageFileReader<ImageType>                ImageReaderType;
  typedef itk::tube::ImageToTubeRigidRegistration<ImageType, TubeNetType, TubeType>
                                                         RegistrationMethodType;

  // read image
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName( argv[1] );
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Gaussian blur the original input image to increase the likelihood of vessel
  // spatial object overlapping with the vessel image at their initial alignment.
  // this enlarges the convergence zone.
  typedef itk::RecursiveGaussianImageFilter
    <ImageType, ImageType> GaussianBlurFilterType;

  GaussianBlurFilterType::Pointer blurFilters[3];
  for( int i = 0; i < 3; i++ )
    {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma( 3.0 );
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection( i );
    }
  blurFilters[0]->SetInput( imageReader->GetOutput() );
  blurFilters[1]->SetInput( blurFilters[0]->GetOutput() );
  blurFilters[2]->SetInput( blurFilters[1]->GetOutput() );
  try
    {
    blurFilters[0]->Update();
    blurFilters[1]->Update();
    blurFilters[2]->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // read vessel
  TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName( argv[2] );
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // subsample points in vessel
  typedef itk::tube::SubSampleTubeTreeSpatialObjectFilter< TubeNetType, TubeType >
    SubSampleTubeNetFilterType;
  SubSampleTubeNetFilterType::Pointer subSampleTubeNetFilter =
    SubSampleTubeNetFilterType::New();
  subSampleTubeNetFilter->SetInput( vesselReader->GetGroup() );
  subSampleTubeNetFilter->SetSampling( 30 );
  try
    {
    subSampleTubeNetFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // register the vessel and the image
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize( 137593424 );


  RegistrationMethodType::Pointer registrationMethod =
    RegistrationMethodType::New();

  registrationMethod->SetFixedImage( blurFilters[2]->GetOutput() );
  registrationMethod->SetMovingSpatialObject( subSampleTubeNetFilter->GetOutput() );

  // Set Optimizer parameters.
  RegistrationMethodType::OptimizerType::Pointer optimizer =
    registrationMethod->GetOptimizer();
  RegistrationMethodType::OptimizerType::ScalesType parameterScales( 6 );
  parameterScales[0] = 30.0;
  parameterScales[1] = 30.0;
  parameterScales[2] = 30.0;
  parameterScales[3] = 1.0;
  parameterScales[4] = 1.0;
  parameterScales[5] = 1.0;
  optimizer->SetScales( parameterScales );
  itk::GradientDescentOptimizer * gradientDescentOptimizer =
    dynamic_cast< itk::GradientDescentOptimizer * >( optimizer.GetPointer() );
  if( gradientDescentOptimizer )
    {
    gradientDescentOptimizer->SetLearningRate( 0.1 );
    gradientDescentOptimizer->SetNumberOfIterations( 1000 );
    }


  // Add a time probe
  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  char timeFile[250] = "";
  char valuesFile[250] = "";
  std::strcat ( timeFile, argv[3] );
  std::strcat ( valuesFile, argv[3] );
  std::strcat ( timeFile, "TimeMemory.txt" );
  std::strcat ( valuesFile, "Values.txt" );

  // Create stream to record the measure
  std::ofstream measuresFile;
  measuresFile.open( timeFile );
  if( !measuresFile.is_open() )
    {
    std::cerr << "Unable to open: " << timeFile << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    memorymeter.Start( "Registration" );
    chronometer.Start( "Registration" );

    registrationMethod->Initialize();

    // Create the Command observer and register it with the optimizer.
    //
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    observer->SetFileName( valuesFile );

    registrationMethod->GetOptimizer()->
      AddObserver( itk::IterationEvent(), observer );

    // Launch Registration
    registrationMethod->Update();

    chronometer.Stop( "Registration" );
    memorymeter.Stop( "Registration" );

    // Report the time and memory taken by the registration
    chronometer.Report( measuresFile );
    memorymeter.Report( measuresFile );

    std::cout << "Optimizer stop condition = "
              << registrationMethod->GetOptimizer()->GetStopConditionDescription()
              << std::endl;
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  measuresFile.close();
  return EXIT_SUCCESS;
}
