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
#include "itkAnisotropicDiffusiveRegistrationFilter.h"
#include "itkVectorCastImageFilter.h"
#include "itkWarpImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "tubeAnisotropicDiffusiveDeformableRegistrationCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  // of your algorithm.
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

  //--------------------------------------------------------
  typedef itk::AnisotropicDiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, VectorImageType > RegistrationType;
  typename RegistrationType::Pointer registrator = RegistrationType::New();

  timeCollector.Start( "Loading input data" );
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

  VectorImageType::Pointer initField = VectorImageType::New();
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

  typedef itk::VectorCastImageFilter< VectorImageType, VectorImageType >
      CasterType;
  CasterType::Pointer caster = CasterType::New();
  caster->SetInput( initField );
  caster->InPlaceOff();
  registrator->SetInitialDeformationField( caster->GetOutput() );

  if( organBoundaryFileName != "" )
    {
    // do we have .vtk or .vtp models?
    std::string::size_type loc = organBoundaryFileName.find_last_of(".");
    if( loc == std::string::npos )
      {
      tube::ErrorMessage( "Failed to find an extension for organ boundary" );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    std::string extension = organBoundaryFileName.substr(loc);

    typename RegistrationType::BorderSurfacePointer borderSurface = NULL;

    if ( extension == std::string(".vtk") )
      {
      vtkSmartPointer< vtkPolyDataReader > polyDataReader
          = vtkPolyDataReader::New();
      polyDataReader->SetFileName( organBoundaryFileName.c_str() );
      polyDataReader->Update();
      borderSurface = polyDataReader->GetOutput();
      }
    else if ( extension == std::string(".vtp") )
      {
      vtkSmartPointer< vtkXMLPolyDataReader > polyDataReader
          = vtkXMLPolyDataReader::New();
      polyDataReader->SetFileName( organBoundaryFileName.c_str() );
      polyDataReader->Update();
      borderSurface = polyDataReader->GetOutput();
      }
    if( !borderSurface )
      {
      tube::ErrorMessage( "Reading polydata: unsuccessful" );
      timeCollector.Report();
      return EXIT_FAILURE;
      }

    // If the world coordinate system is RAS (i.e. called from 3D Slicer)
    // then the model will be in RAS space while the images will be in LPS
    // space.  It's easiest to transform the model to LPS space.
    if ( worldCoordinateSystem == "RAS" )
      {
      vtkSmartPointer< vtkTransform > RAStoLPS = vtkTransform::New();
      RAStoLPS->RotateX(180); // flip in right-left
      RAStoLPS->RotateY(180); // flip in anterior-posterior
      vtkSmartPointer< vtkTransformPolyDataFilter > transformPolyDataFilter
          = vtkTransformPolyDataFilter::New();
      transformPolyDataFilter->SetInput( borderSurface );
      transformPolyDataFilter->SetTransform( RAStoLPS );
      transformPolyDataFilter->Update();
      borderSurface = transformPolyDataFilter->GetOutput();
      if( !borderSurface )
        {
        tube::ErrorMessage( "Transforming polydata: unsuccessful" );
        timeCollector.Report();
        return EXIT_FAILURE;
        }
    }

    registrator->SetBorderSurface( borderSurface );
    }

  if( inputNormalVectorImageFileName != "" )
    {
    typedef itk::ImageFileReader< VectorImageType > VectorImageReaderType;
    typename VectorImageReaderType::Pointer vectorImageReader
        = VectorImageReaderType::New();
    vectorImageReader->SetFileName( inputNormalVectorImageFileName.c_str() );
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

  if( inputWeightImageFileName != "" )
    {
    typedef itk::ImageFileReader< WeightImageType > WeightImageReaderType;
    typename WeightImageReaderType::Pointer weightImageReader
        = WeightImageReaderType::New();
    weightImageReader->SetFileName( inputWeightImageFileName.c_str() );
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

  timeCollector.Stop( "Loading input data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  //-------------------------------------------------------------
  timeCollector.Start( "Registration" );
  registrator->SetNumberOfIterations( numberOfIterations );
  registrator->SetComputeRegularizationTerm( !doNotPerformRegularization );
  registrator->SetUseAnisotropicRegularization(
      !doNotUseAnisotropicRegularization );

  registrator->SetTimeStep( timeStep );
  if( lambda > 0.0 )
    {
    tube::ErrorMessage( "Lambda must be negative." );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  registrator->SetLambda( lambda );

  tube::CLIFilterWatcher watchRegistration(registrator,
                                           "Anisotropic Diffusive Registration",
                                           CLPProcessInformation,
                                           0.8,
                                           progress );

  // warp moving image
  typedef
      itk::WarpImageFilter< MovingImageType, MovingImageType, VectorImageType >
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
    typedef itk::ImageFileWriter< VectorImageType > FieldWriterType;
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

//  if( outputTransformFileName != "" )
//    {
//    typedef itk::ImageFileWriter< VectorImageType > GridWriterType;
//    GridWriterType::Pointer gridWriter = GridWriterType::New();
//    gridWriter->SetFileName( outputTransformFileName );
//    gridWriter->SetInput( registrator->GetOutput() );
//    try
//      {
//      gridWriter->Update();
//      }
//    catch( itk::ExceptionObject & err )
//      {
//      tube::ErrorMessage( "Writing volume: Exception caught: "
//                          + std::string(err.GetDescription()) );
//      timeCollector.Report();
//      return EXIT_FAILURE;
//      }
//    }

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

  if( outputNormalVectorImageFileName != "" )
    {
    typedef itk::ImageFileWriter< VectorImageType > VectorWriterType;
    typename VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
    vectorWriter->SetFileName( outputNormalVectorImageFileName );
    vectorWriter->SetInput( registrator->GetNormalVectorImage() );
    try
      {
      vectorWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

  if( outputWeightImageFileName != "" )
    {
    typedef itk::ImageFileWriter< WeightImageType > WeightWriterType;
    typename WeightWriterType::Pointer weightWriter = WeightWriterType::New();
    weightWriter->SetFileName( outputWeightImageFileName );
    weightWriter->SetInput( registrator->GetWeightImage() );
    try
      {
      weightWriter->Update();
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
