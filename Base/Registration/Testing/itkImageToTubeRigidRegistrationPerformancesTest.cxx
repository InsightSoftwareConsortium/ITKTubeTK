/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkImageToTubeRigidRegistration.h"
#include "itkSpatialObjectReader.h"
#include "itkImageFileReader.h"
#include "itkTubeToTubeTransformFilter.h"
#include "itkEuler3DTransform.h"
#include "itkMath.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"

// STD includes
#include <iostream>
#include <fstream>

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
#include "itkCommand.h"
class CommandIterationUpdate : public itk::Command
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>   Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate(){}
  ~CommandIterationUpdate()
    {
    if ( measuresFileStream.is_open() )
      {
      measuresFileStream.close();
      }
    }

  std::ofstream measuresFileStream;

public:
  typedef itk::GradientDescentOptimizer OptimizerType;
  typedef const OptimizerType *         OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

  void Execute(const itk::Object * object, const itk::EventObject & event)
    {
    OptimizerPointer optimizer =
      dynamic_cast< OptimizerPointer >( object );
    if( !(itk::IterationEvent().CheckEvent( &event )) )
      {
      return;
      }

    if ( measuresFileStream.is_open() )
      {
      measuresFileStream << optimizer->GetCurrentIteration() << ";";
      measuresFileStream << optimizer->GetValue() << std::endl;
      }
    }

  void SetFileName( char* fileName )
    {
    measuresFileStream.open( fileName, std::ios::trunc );
    }
};

int itkImageToTubeRigidRegistrationPerformancesTest(int argc, char* argv [] )
{
  if ( argc < 4 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image " << "Input_Vessel " << "Output_Results"
              << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::GroupSpatialObject<3>                      TubeNetType;
  typedef itk::SpatialObjectReader<3>                     TubeNetReaderType;
  typedef itk::Image<double, 3>                           ImageType;
  typedef itk::ImageFileReader<ImageType>                 ImageReaderType;
  typedef itk::ImageToTubeRigidRegistration
    <ImageType, TubeNetType> RegistrationFilterType;

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

  // guassian blur the original input image to increase the likelihood of vessel
  // spatial object overlapping with the vessel image at their initial alignment.
  // this enlarges the convergence zone.
  typedef itk::RecursiveGaussianImageFilter
    <ImageType, ImageType> GaussianBlurFilterType;

  GaussianBlurFilterType::Pointer blurFilters[3];
  for (int i = 0; i < 3; i++)
    {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma( 3.0 );
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection( i );
    }
  blurFilters[0]->SetInput(imageReader->GetOutput());
  blurFilters[1]->SetInput(blurFilters[0]->GetOutput());
  blurFilters[2]->SetInput(blurFilters[1]->GetOutput());
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

  // register the vessel and the image
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize( 137593424 );

  double parameterScales[6] = {30.0, 30.0, 30.0, 1.0, 1.0, 1.0};
  double initialPose[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  RegistrationFilterType::Pointer registrationFilter =
    RegistrationFilterType::New();

  registrationFilter->SetFixedImage( blurFilters[2]->GetOutput() );
  registrationFilter->SetMovingSpatialObject( vesselReader->GetGroup() );
  registrationFilter->SetNumberOfIteration( 1000 );
  registrationFilter->SetLearningRate( 0.1 );
  registrationFilter->SetInitialPosition( initialPose );
  registrationFilter->SetParametersScale( parameterScales );
  registrationFilter->SetVerbose( false );
  registrationFilter->SetSampling( 50 );

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
  if ( !measuresFile.is_open() )
    {
    std::cerr << "Unable to open: " << timeFile << std::endl;
    return EXIT_FAILURE;
    }

  try
    {
    memorymeter.Start( "Registration" );
    chronometer.Start( "Registration" );

    registrationFilter->Initialize();

    // Create the Command observer and register it with the optimizer.
    //
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    observer->SetFileName( valuesFile );

    registrationFilter->GetOptimizer()->
      AddObserver( itk::IterationEvent(), observer );

    // Launch Registration
    registrationFilter->Update();

    chronometer.Stop( "Registration" );
    memorymeter.Stop( "Registration" );

    // Report the time and memory taken by the registration
    chronometer.Report( measuresFile );
    memorymeter.Report( measuresFile );

    std::cout << "Optimizer stop condition = "
              << registrationFilter->GetOptimizer()->GetStopConditionDescription()
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
