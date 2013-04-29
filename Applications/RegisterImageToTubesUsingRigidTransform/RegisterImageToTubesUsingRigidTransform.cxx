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

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkEuler3DTransform.h"
#include "itkImageToTubeRigidRegistration.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkTubeToTubeTransformFilter.h"
#include "itkVesselTubeSpatialObject.h"
#include "itkSubSampleTubeTreeSpatialObjectFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "RegisterImageToTubesUsingRigidTransformCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "SampleCLIApplication",
                                                 CLPProcessInformation );
  progressReporter.Start();

  static const unsigned int Dimension = 3;
  typedef double            FloatType;

  typedef itk::VesselTubeSpatialObject< Dimension >      TubeType;
  typedef itk::GroupSpatialObject< Dimension >           TubeNetType;
  typedef itk::SpatialObjectReader< Dimension >          TubeNetReaderType;
  typedef itk::Image< FloatType, Dimension >             ImageType;
  typedef itk::ImageFileReader< ImageType >              ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >              ImageWriterType;
  typedef itk::ImageToTubeRigidRegistration< ImageType, TubeNetType, TubeType >
                                                         RegistrationFilterType;
  typedef RegistrationFilterType::TransformType          TransformType;
  typedef itk::TubeToTubeTransformFilter< TransformType, Dimension >
                                                         TubeTransformFilterType;
  typedef itk::SubSampleTubeTreeSpatialObjectFilter< TubeNetType, TubeType >
                                                         SubSampleTubeNetFilterType;

  timeCollector.Start("Load data");
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  TubeNetReaderType::Pointer vesselReader = TubeNetReaderType::New();
  vesselReader->SetFileName( inputVessel );
  try
    {
    vesselReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading vessel: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Sub-sample data");
  SubSampleTubeNetFilterType::Pointer subSampleTubeNetFilter =
    SubSampleTubeNetFilterType::New();
  subSampleTubeNetFilter->SetInput( vesselReader->GetGroup() );
  subSampleTubeNetFilter->SetSampling( 100 );
  try
    {
    subSampleTubeNetFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Sub-sampling vessel: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  timeCollector.Stop("Sub-sample data");
  progress = 0.2;
  progressReporter.Report( progress );


  typename ImageType::Pointer currentImage = reader->GetOutput();
  if( gaussianBlurStdDev > 0.0 )
    {
    timeCollector.Start("Gaussian Blur");

    typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType >
      GaussianFilterType;
    GaussianFilterType::Pointer gaussianFilter;

    // Progress per iteration
    const double progressFraction = 0.1/Dimension;
    for( unsigned int ii = 0; ii < Dimension; ++ii )
      {
      gaussianFilter = GaussianFilterType::New();
      gaussianFilter->SetInput( currentImage );
      gaussianFilter->SetSigma( gaussianBlurStdDev );

      gaussianFilter->SetOrder(
               itk::RecursiveGaussianImageFilter<ImageType>::ZeroOrder );
      gaussianFilter->SetDirection( ii );
      tube::CLIFilterWatcher watcher( gaussianFilter,
                                      "Blur Filter 1D",
                                      CLPProcessInformation,
                                      progressFraction,
                                      progress,
                                      true );

      gaussianFilter->Update();
      currentImage = gaussianFilter->GetOutput();
      }

    timeCollector.Stop("Gaussian Blur");
    }


  timeCollector.Start("Register image to tube");

  RegistrationFilterType::Pointer  registrationFilter =
    RegistrationFilterType::New();

  registrationFilter->SetFixedImage( currentImage );
  registrationFilter->SetMovingSpatialObject( subSampleTubeNetFilter->GetOutput() );
  registrationFilter->SetNumberOfIteration( 1000 );
  registrationFilter->SetLearningRate( 0.1 );
  try
    {
    registrationFilter->Initialize();
    registrationFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Performing registration: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  progress = 0.9;
  progressReporter.Report( progress );
  timeCollector.Stop("Register image to tube");


  timeCollector.Start("Save data");
  TransformType* outputTransform =
    dynamic_cast<TransformType *>(registrationFilter->GetTransform());
  outputTransform->SetParameters( registrationFilter->GetLastTransformParameters() );

  TubeTransformFilterType::Pointer transformFilter = TubeTransformFilterType::New();
  transformFilter->SetInput( vesselReader->GetGroup() );
  transformFilter->SetScale( 1.0 );
  transformFilter->SetTransform( outputTransform );
  try
    {
    transformFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Transforming tube: Exception caught: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typedef itk::SpatialObjectToImageFilter<TubeNetType, ImageType>
                                              SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer vesselToImageFilter =
    SpatialObjectToImageFilterType::New();
  ImageType::SizeType size = currentImage->GetLargestPossibleRegion().GetSize();
  const double decimationFactor = 4.0;
  typedef ImageType::SizeType::SizeValueType SizeValueType;
  size[0] = static_cast< SizeValueType >( size[0] / decimationFactor );
  size[1] = static_cast< SizeValueType >( size[1] / decimationFactor );
  size[2] = static_cast< SizeValueType >( size[2] / decimationFactor );

  ImageType::SpacingType spacing = currentImage->GetSpacing();
  spacing[0] = spacing[0] * decimationFactor;
  spacing[1] = spacing[1] * decimationFactor;
  spacing[2] = spacing[2] * decimationFactor;
  vesselToImageFilter->SetInput( transformFilter->GetOutput() );
  vesselToImageFilter->SetSize( size );
  vesselToImageFilter->SetSpacing( spacing );
  vesselToImageFilter->SetOrigin( currentImage->GetOrigin() );
  vesselToImageFilter->SetInsideValue( 1.0 );
  vesselToImageFilter->SetOutsideValue( 0.0 );
  vesselToImageFilter->Update();

  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( vesselToImageFilter->GetOutput() );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End( );

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
