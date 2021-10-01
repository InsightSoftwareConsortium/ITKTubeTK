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

#include "RegisterSpatialObjectsToImageCLP.h"

#include "itktubeSpatialObjectToImageRegistrationHelper.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkSpatialObjectReader.h"
#include "itkSpatialObjectWriter.h"
#include "itkMultiThreaderBase.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"


template < class TPixelType, unsigned int TDimension >
int DoIt( int argc, char * argv[] )
{

  PARSE_ARGS;


  enum VerboseLevelEnum { SILENT, STANDARD, VERBOSE };
  VerboseLevelEnum verbosity = SILENT;
  if( verbosityLevel == "Standard" )
    {
    verbosity = STANDARD;
    }
  else if( verbosityLevel == "Verbose" )
    {
    verbosity = VERBOSE;
    }

  typedef typename itk::Image< TPixelType, TDimension > ImageType;

  typedef itk::ImageFileReader< ImageType  > ImageReaderType;
  typedef itk::SpatialObjectReader< TDimension, float > SpatialObjectReaderType;

  typedef typename itk::tube::SpatialObjectToImageRegistrationHelper<
    TDimension, ImageType > RegistrationType;

  typename RegistrationType::Pointer reger = RegistrationType::New();

  typedef typename RegistrationType::SpatialObjectType SpatialObjectType;

  reger->SetReportProgress( true );

  if( verbosity >= STANDARD )
    {
    std::cout << "###Loading fixed image...";
    }
  typename ImageReaderType::Pointer imgReader = ImageReaderType::New();
  imgReader->SetFileName( fixedImage );
  imgReader->Update();
  reger->SetFixedImage( imgReader->GetOutput() );
  if( verbosity >= STANDARD )
    {
    std::cout << "###DONE" << std::endl;
    }

  if( verbosity >= STANDARD )
    {
    std::cout << "###Loading moving spatial object...";
    }
  typename SpatialObjectReaderType::Pointer soReader =
    SpatialObjectReaderType::New();
  soReader->SetFileName( movingSpatialObject );
  soReader->Update();
  reger->SetMovingSpatialObject( soReader->GetOutput() );
  if( verbosity >= STANDARD )
    {
    std::cout << "###DONE" << std::endl;
    }

  if( loadTransform.size() > 1 )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Loading transform...";
      }
    reger->LoadTransform( loadTransform, invertLoadedTransform );
    if( verbosity >= STANDARD )
      {
      std::cout << "###DONE" << std::endl;
      }
    }

  if( fixedLandmarks.size() > 1 || movingLandmarks.size() > 1 )
    {
    if( initialization != "Landmarks" )
      {
      std::cout << "WARNING: Landmarks specified, but initialization "
                << "process was not told to use landmarks. " << std::endl;
      std::cout << "Changing initialization to use landmarks." << std::endl;
      reger->SetInitialMethodEnum( RegistrationType::INIT_WITH_LANDMARKS );
      }
    }
  if( skipInitialRandomSearch )
    {
    reger->SetUseEvolutionaryOptimization( false );
    }
  else
    {
    reger->SetUseEvolutionaryOptimization( true );
    }

  if( initialization == "Landmarks" )
    {
    reger->SetInitialMethodEnum( RegistrationType::INIT_WITH_LANDMARKS );
    reger->SetFixedLandmarks( fixedLandmarks );
    reger->SetMovingLandmarks( movingLandmarks );
    }
  else if( initialization == "ImageCenters" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Initialization: ImageCenters" << std::endl;
      }
    reger->SetInitialMethodEnum(
      RegistrationType::INIT_WITH_IMAGE_CENTERS );
    }
  else if( initialization == "CentersOfMass" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Initialization: CentersOfMass" << std::endl;
      }
    reger->SetInitialMethodEnum(
      RegistrationType::INIT_WITH_CENTERS_OF_MASS );
    }
  else // if( initialization == "None" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Initialization: None" << std::endl;
      }
    reger->SetInitialMethodEnum( RegistrationType::INIT_WITH_NONE );
    }

  if( registration == "None" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: None" << std::endl;
      }
    reger->SetEnableInitialRegistration( false );
    reger->SetEnableRigidRegistration( false );
    reger->SetEnableAffineRegistration( false );
    }
  else if( registration == "Initial" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: Initial" << std::endl;
      }
    reger->SetEnableInitialRegistration( true );
    reger->SetEnableRigidRegistration( false );
    reger->SetEnableAffineRegistration( false );
    }
  else if( registration == "Rigid" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: Rigid" << std::endl;
      }
    reger->SetEnableInitialRegistration( false );
    reger->SetEnableRigidRegistration( true );
    reger->SetEnableAffineRegistration( false );
    }
  else if( registration == "Affine" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: Affine" << std::endl;
      }
    reger->SetEnableInitialRegistration( false );
    reger->SetEnableRigidRegistration( false );
    reger->SetEnableAffineRegistration( true );
    }
  else if( registration == "PipelineRigid" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: PipelineRigid" << std::endl;
      }
    reger->SetEnableInitialRegistration( true );
    reger->SetEnableRigidRegistration( true );
    reger->SetEnableAffineRegistration( false );
    }
  else if( registration == "PipelineAffine" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Registration: PipelineAffine" << std::endl;
      }
    reger->SetEnableInitialRegistration( true );
    reger->SetEnableRigidRegistration( true );
    reger->SetEnableAffineRegistration( true );
    }

  if( metric == "ImageIntensityMetric" )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Metric: ImageIntensityMetric" << std::endl;
      }
    reger->SetRigidMetricMethodEnum( RegistrationType
                                     ::OptimizedRegistrationMethodType
                                     ::IMAGE_INTENSITY_METRIC );
    reger->SetAffineMetricMethodEnum( RegistrationType
                                      ::OptimizedRegistrationMethodType
                                      ::IMAGE_INTENSITY_METRIC );
    }

  typedef typename itk::ImageFileReader<
    itk::Image< unsigned char, TDimension > > ImageReader;
  typedef typename itk::ImageMaskSpatialObject< TDimension >
    ImageMaskSpatialObject;

  if( fixedImageMask != "" )
    {
    reger->SetUseFixedImageMaskObject( true );

    typename ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName( fixedImageMask );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & exception )
      {
      std::cerr << "Exception caught while loading fixed image mask."
        << std::endl;
      std::cerr << exception << std::endl;
      return EXIT_FAILURE;
      }

    typename ImageMaskSpatialObject::Pointer mask =
      ImageMaskSpatialObject::New();
    mask->SetImage( reader->GetOutput() );
    reger->SetFixedImageMaskObject( mask );

    if( verbosity >= STANDARD )
      {
      std::cout << "###useFixedImageMaskObject: true" << std::endl;
      }
    }
  else
    {
    reger->SetUseFixedImageMaskObject( false );
    if( verbosity >= STANDARD )
      {
      std::cout << "###useFixedImageMaskObject: false" << std::endl;
      }
    }


  if( movingSpatialObjectMask != "" )
    {
    reger->SetUseMovingSpatialObjectMaskObject( true );

    typename ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName( movingSpatialObjectMask );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & exception )
      {
      std::cerr << "Exception caught while loading moving image mask."
        << std::endl;
      std::cerr << exception << std::endl;
      return EXIT_FAILURE;
      }

    typename ImageMaskSpatialObject::Pointer mask =
      ImageMaskSpatialObject::New();
    mask->SetImage( reader->GetOutput() );
    reger->SetMovingSpatialObjectMaskObject( mask );

    if( verbosity >= STANDARD )
      {
      std::cout << "###useMovingSpatialObjectMaskObject: true" << std::endl;
      }
    }
  else
    {
    reger->SetUseMovingSpatialObjectMaskObject( false );
    if( verbosity >= STANDARD )
      {
      std::cout << "###useMovingSpatialObjectMaskObject: false" << std::endl;
      }
    }

  reger->SetRandomNumberSeed( randomNumberSeed );

  reger->SetRigidMaxIterations( rigidMaxIterations );
  if( verbosity >= STANDARD )
    {
    std::cout << "###RigidMaxIterations: " << rigidMaxIterations
      << std::endl;
    }

  reger->SetAffineMaxIterations( affineMaxIterations );
  if( verbosity >= STANDARD )
    {
    std::cout << "###AffineMaxIterations: " << affineMaxIterations
      << std::endl;
    }

  reger->SetRigidSamplingRatio( rigidSamplingRatio );
  if( verbosity >= STANDARD )
    {
    std::cout << "###RigidSamplingRatio: " << rigidSamplingRatio
      << std::endl;
    }
  reger->SetAffineSamplingRatio( affineSamplingRatio );
  if( verbosity >= STANDARD )
    {
    std::cout << "###AffineSamplingRatio: " << affineSamplingRatio
      << std::endl;
    }

  reger->SetExpectedOffsetMagnitude( expectedOffset );
  if( verbosity >= STANDARD )
    {
    std::cout << "###ExpectedOffsetPixelMagnitude: " << expectedOffset
      << std::endl;
    }

  reger->SetExpectedRotationMagnitude( expectedRotation );
  if( verbosity >= STANDARD )
    {
    std::cout << "###ExpectedRotationMagnitude: " << expectedRotation
      << std::endl;
    }

  reger->SetExpectedScaleMagnitude( expectedScale );
  if( verbosity >= STANDARD )
    {
    std::cout << "###ExpectedScaleMagnitude: " << expectedScale
      << std::endl;
    }

  reger->SetExpectedSkewMagnitude( expectedSkew );
  if( verbosity >= STANDARD )
    {
    std::cout << "###ExpectedSkewMagnitude: " << expectedSkew
      << std::endl;
    }

  try
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Starting registration..." << std::endl;
      }
    reger->Update();
    }
  catch( itk::ExceptionObject & exception )
    {
    std::cerr << "Exception caught during helper class registration."
      << exception << std::endl;
    std::cerr << "Current Matrix Transform = " << std::endl;
    reger->GetCurrentMatrixTransform()->Print( std::cerr, 2 );
    return EXIT_FAILURE;
    }
  catch( ... )
    {
    std::cerr << "Uncaught exception during helper class registration."
      << std::endl;
    return EXIT_FAILURE;
    }

  if( resampledSpatialObject.size() > 1 )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###Resampling..." << std::endl;
      }
    typename SpatialObjectType::ConstPointer resultSpatialObject;
    try
      {
      resultSpatialObject = reger->ResampleSpatialObject( NULL, NULL,
        resampledPortion );
      }
    catch( itk::ExceptionObject & exception )
      {
      std::cerr << "Exception caught during helper class resampling."
        << exception << std::endl;
      std::cerr << "Current Matrix Transform = " << std::endl;
      reger->GetCurrentMatrixTransform()->Print( std::cerr, 2 );
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr << "Uncaught exception during helper class resampling."
                << std::endl;
      return EXIT_FAILURE;
      }

    try
      {
      typedef itk::SpatialObjectWriter< TDimension > SOWriterType;
      typename SOWriterType::Pointer soWriter =
        SOWriterType::New();
      soWriter->SetFileName( resampledSpatialObject );
      soWriter->SetInput( resultSpatialObject );
      soWriter->Update();
      }
    catch( itk::ExceptionObject & exception )
      {
      std::cerr <<
        "Exception caught during helper class resampled image saving."
        << exception << std::endl;
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr <<
        "Uncaught exception during helper class resampled image saving."
        << std::endl;
      return EXIT_FAILURE;
      }
    }

  if( saveTransform.size() > 1 )
    {
    try
      {
      reger->SaveTransform( saveTransform );
      }
    catch( itk::ExceptionObject & exception )
      {
      std::cerr << "Exception caught during helper class transform saving."
                << exception << std::endl;
      return EXIT_FAILURE;
      }
    catch( ... )
      {
      std::cerr << "Uncaught exception during helper class saving."
                << std::endl;
      return EXIT_FAILURE;
      }
    }

  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  enum VerboseLevelEnum { SILENT, STANDARD, VERBOSE };
  VerboseLevelEnum verbosity = SILENT;
  if( verbosityLevel == "Standard" )
    {
    verbosity = STANDARD;
    }
  else if( verbosityLevel == "Verbose" )
    {
    verbosity = VERBOSE;
    }

  if( numberOfThreads != 0 )
    {
    if( verbosity >= STANDARD )
      {
      std::cout << "###numberOfThreads: " << numberOfThreads << std::endl;
      }
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads( numberOfThreads );
    }

  tube::ParseArgsAndCallDoIt( fixedImage, argc, argv );

  return EXIT_SUCCESS;
}
