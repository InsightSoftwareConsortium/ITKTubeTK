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

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// It is important to use OrientedImages
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkImageToImageAnisotropicDiffusiveDeformableRegistrationFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "AnisotropicDiffusiveDeformableRegistrationCLP.h"

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
  tube::CLIProgressReporter    progressReporter(
      "AnisotropicDiffusiveDeformableRegistration", CLPProcessInformation );
  progressReporter.Start();

  // Typedefs
  const unsigned int                                      ImageDimension = 3;
  typedef pixelT                                          FixedPixelType;
  typedef pixelT                                          MovingPixelType;
  typedef MovingPixelType                                 OutputPixelType;
  typedef double                                          VectorScalarType;
  typedef itk::Image< FixedPixelType, ImageDimension >    FixedImageType;
  typedef itk::Image< MovingPixelType, ImageDimension >   MovingImageType;
  typedef itk::Image< OutputPixelType, ImageDimension >   OutputImageType;
  typedef itk::Vector< VectorScalarType, ImageDimension > VectorType;
  typedef itk::Image< VectorType, ImageDimension >        VectorImageType;
  typedef itk::Image< double, ImageDimension >            WeightImageType;
  typedef itk::Image< VectorType, ImageDimension >        FieldType;

  //--------------------------------------------------------
  typedef itk::ImageToImageAnisotropicDiffusiveDeformableRegistrationFilter
      < FixedImageType, MovingImageType, FieldType > RegistrationType;
  typename RegistrationType::Pointer registrator = RegistrationType::New();

  timeCollector.Start( "Load data" );
  typedef itk::ImageFileReader< FixedImageType > FixedImageReaderType;
  typename FixedImageReaderType::Pointer fixedImageReader
      = FixedImageReaderType::New();
  fixedImageReader->SetFileName( fixedImageFileName.c_str() );
  try
    {
    fixedImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename FixedImageType::Pointer fixed = fixedImageReader->GetOutput();
  registrator->SetFixedImage( fixed );

  typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;
  typename MovingImageReaderType::Pointer movingImageReader
      = MovingImageReaderType::New();
  movingImageReader->SetFileName( movingImageFileName.c_str() );
  try
    {
    movingImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename MovingImageType::Pointer moving = movingImageReader->GetOutput();
  registrator->SetMovingImage( moving );

  FieldType::Pointer initField = FieldType::New();
  initField->SetSpacing( fixed->GetSpacing() );
  initField->SetOrigin( fixed->GetOrigin() );
  initField->SetLargestPossibleRegion( fixed->GetLargestPossibleRegion() );
  initField->SetRequestedRegion( fixed->GetRequestedRegion() );
  initField->SetBufferedRegion( fixed->GetBufferedRegion() );
  initField->Allocate();

  // fill initial deformation with zero vectors
  VectorType zeroVec;
  zeroVec.Fill( 0.0 );
  initField->FillBuffer( zeroVec );

  typedef itk::VectorCastImageFilter< FieldType, FieldType > CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();
  registrator->SetInitialDeformationField( caster->GetOutput() );

  if( organBoundaryFileName != "" )
    {
    vtkSmartPointer< vtkPolyDataReader > polyDataReader = vtkPolyDataReader::New();
    polyDataReader->SetFileName( organBoundaryFileName.c_str() );
    polyDataReader->Update();
    if( !polyDataReader->GetOutput() )
      {
      tube::ErrorMessage( "Reading polydata: unsuccessful" );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    registrator->SetBorderSurface( polyDataReader->GetOutput() );
    }

  if( normalVectorImageFileName != "" )
    {
    typedef itk::ImageFileReader< VectorImageType > VectorImageReaderType;
    typename VectorImageReaderType::Pointer vectorImageReader
        = VectorImageReaderType::New();
    vectorImageReader->SetFileName( normalVectorImageFileName.c_str() );
    try
      {
      vectorImageReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    registrator->SetNormalVectorImage( vectorImageReader->GetOutput() );
    }

  if( weightImageFileName != "" )
    {
    typedef itk::ImageFileReader< WeightImageType > WeightImageReaderType;
    typename WeightImageReaderType::Pointer weightImageReader
        = WeightImageReaderType::New();
    weightImageReader->SetFileName( weightImageFileName.c_str() );
    try
      {
      weightImageReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    registrator->SetWeightImage( weightImageReader->GetOutput() );
    }

  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  //-------------------------------------------------------------
  timeCollector.Start( "Registration" );
  registrator->SetNumberOfIterations( numberOfIterations );
  registrator->SetComputeRegularizationTerm( !doNotPerformRegularization );
  registrator->SetUseAnisotropicRegularization(
      !doNotUseAnisotropicRegularization );

  std::cout << "compute regularization " << registrator->GetComputeRegularizationTerm() << std::endl;
  std::cout << "perform aniso " << registrator->GetUseAnisotropicRegularization() << std::endl;

  registrator->SetTimeStep( timeStep );
  if( lambda > 0.0 )
    {
    tube::ErrorMessage( "Lambda must be negative." );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  registrator->SetLambda( lambda );

  // warp moving image
  typedef itk::WarpImageFilter< MovingImageType, MovingImageType, FieldType >
      WarperType;
  typename WarperType::Pointer warper = WarperType::New();

  typedef typename WarperType::CoordRepType CoordRepType;
  typedef itk::LinearInterpolateImageFunction<MovingImageType,CoordRepType>
      InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

  warper->SetInput( moving );
  warper->SetDeformationField( registrator->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputSpacing( fixed->GetSpacing() );
  warper->SetOutputOrigin( fixed->GetOrigin() );
  warper->SetOutputDirection( fixed->GetDirection() );
  warper->SetEdgePaddingValue( 0 );

  // Update triggers the registration
  warper->Update();

  timeCollector.Stop( "Registration" );
  progress = 0.9;
  progressReporter.Report( progress );

  // ---------------------------------------------------------
  timeCollector.Start( "Write outputs" );

  if( outputDeformationFieldFileName != "" )
    {
    typedef itk::ImageFileWriter< FieldType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( outputDeformationFieldFileName );
    fieldWriter->SetInput( registrator->GetOutput() );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

  if( outputResampledImageFileName != "" )
    {
    typedef itk::ImageFileWriter< MovingImageType > ImageWriterType;
    typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetFileName( outputResampledImageFileName );
    imageWriter->SetInput( warper->GetOutput() );
    try
      {
      imageWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

  timeCollector.Stop( "Write outputs" );
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
  return tube::ParseArgsAndCallDoIt( fixedImageFileName, argc, argv );
}
