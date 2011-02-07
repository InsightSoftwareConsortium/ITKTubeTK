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

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "itkAnisotropicDiffusiveRegistrationFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkOrientImageFilter.h"
#include "itkTransform.h"
#include "itkTransformFileReader.h"
#include "itkWarpImageFilter.h"

#include "vtkPolyDataReader.h"
#include "vtkSmartPointer.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkXMLPolyDataReader.h"

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
  tube::CLIProgressReporter progressReporter(
      "AnisotropicDiffusiveDeformableRegistration", CLPProcessInformation );
  progressReporter.Start();

  // Typedefs seeding definition of the anisotropic diffusive registration
  // filter
  const unsigned int                                      ImageDimension = 3;
  typedef pixelT                                          FixedPixelType;
  typedef pixelT                                          MovingPixelType;
  typedef MovingPixelType                                 OutputPixelType;
  typedef itk::Image< FixedPixelType, ImageDimension >    FixedImageType;
  typedef itk::Image< MovingPixelType, ImageDimension >   MovingImageType;
  typedef double                                          VectorScalarType;
  typedef itk::Vector< VectorScalarType, ImageDimension > VectorType;
  typedef itk::Image< VectorType, ImageDimension >        VectorImageType;

  // Initialize anisotropic diffusive registration filter
  typedef itk::AnisotropicDiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, VectorImageType > RegistrationType;
  typename RegistrationType::Pointer registrator = RegistrationType::New();
  typedef typename RegistrationType::WeightImageType      WeightImageType;

  // Load the fixed image
  timeCollector.Start( "Loading fixed image" );
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
  timeCollector.Stop( "Loading fixed image" );

  // Load the moving image
  timeCollector.Start( "Loading moving image" );
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
  timeCollector.Stop( "Loading moving image" );

  // Setup the initial deformation field: the initial deformation field should
  // be in the space of the fixed image
  timeCollector.Start( "Setup initial deformation field" );
  VectorImageType::Pointer initField = VectorImageType::New();
  typename FixedImageType::Pointer templateImage
      = fixedImageReader->GetOutput();
  initField->SetOrigin( templateImage->GetOrigin() );
  initField->SetSpacing( templateImage->GetSpacing() );
  initField->SetDirection( templateImage->GetDirection() );
  initField->SetLargestPossibleRegion( templateImage->GetLargestPossibleRegion() );
  initField->SetRequestedRegion( templateImage->GetRequestedRegion() );
  initField->SetBufferedRegion( templateImage->GetBufferedRegion() );
  initField->Allocate();

  // Use the initial transform if given
  if( initialTransform != "" )
    {
    typedef itk::TransformFileReader TransformReaderType;
    TransformReaderType::Pointer transformReader = TransformReaderType::New();
    transformReader->SetFileName( initialTransform );
    try
      {
      transformReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading initial transform: Exception caught: "
                          + std::string(err.GetDescription()) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    if( transformReader->GetTransformList()->size() != 0 )
      {
      TransformReaderType::TransformType::Pointer initial
          = *( transformReader->GetTransformList()->begin() );

      // Cast to transform pointer, so that we can use TransformPoint()
      typedef itk::Transform< double, ImageDimension, ImageDimension >
          TransformType;
      typename TransformType::Pointer transform
          = dynamic_cast< TransformType* >( initial.GetPointer() );

      // For each voxel, find the displacement invoked by the given initial
      // transform.  This should work for all types of transforms (linear,
      // nonlinear, B-spline, etc) because itk::Transform is the base for each.
      // Slicer saves transforms in LPS space (see Slicer's
      // vtkMRMLTransformStorageNode::WriteData()), so we don't need to modify
      // the initial transform according to the worldCoordinateSystem
      // variable (unlike the surface model).
      if( transform )
        {
        typename TransformType::InputPointType physicalPoint;
        physicalPoint.Fill( 0 );
        typename TransformType::OutputPointType transformedPoint;
        transformedPoint.Fill( 0 );
        // Initial displacement vector
        VectorType initVector;
        initVector.Fill( 0 );
        typedef itk::ImageRegionIterator< VectorImageType >
            VectorImageRegionType;
        VectorImageRegionType initIt = VectorImageRegionType(
            initField, initField->GetLargestPossibleRegion() );
        for( initIt.GoToBegin(); !initIt.IsAtEnd(); ++initIt )
          {
          initField->TransformIndexToPhysicalPoint( initIt.GetIndex(),
                                                    physicalPoint );
          transformedPoint = transform->TransformPoint( physicalPoint );
          for( unsigned int i = 0; i < ImageDimension; i++ )
            {
            initVector[i] = transformedPoint[i] - physicalPoint[i];
            }
          initIt.Set( initVector );
          }
        }
      else
        {
        tube::ErrorMessage( "Initial transform is an unsupported type" );
        timeCollector.Report();
        return EXIT_FAILURE;
        }
      }
    }
  // If no initial transform is given, fill the initial field with zero vectors
  else
    {
    VectorType zeroVector;
    zeroVector.Fill( 0.0 );
    initField->FillBuffer( zeroVector );
    }
  timeCollector.Stop( "Setup initial deformation field" );

  // Orient to axials to avoid issue with the anisotropic diffusive registration
  // function not taking image direction into account.
  //
  // Forcing image to be axials avoids this problem. Note, that
  // reorientation only affects the internal mapping from index to
  // physical coordinates.  The reoriented data spans the same
  // physical space as the original data.  Thus, the registration
  // transform calculated on the reoriented data is also the
  // transform for the original un-reoriented data.
  //
  // We wait to orient the fixed, moving and initial transform images
  // until after we have made the initial transform image because the
  // reorientation affects the internal mapping from index to physical
  // coordinates.
  timeCollector.Start( "Orient fixed image" );
  typedef itk::OrientImageFilter< FixedImageType, FixedImageType >
      FixedOrientFilterType;
  typename FixedOrientFilterType::Pointer orientFixed
      = FixedOrientFilterType::New();
  orientFixed->UseImageDirectionOn();
  orientFixed->SetDesiredCoordinateOrientationToAxial();
  orientFixed->SetInput( fixedImageReader->GetOutput() );
  orientFixed->Update();
  timeCollector.Stop( "Orient fixed image" );

  timeCollector.Start( "Orient moving image" );
  typedef itk::OrientImageFilter< MovingImageType, MovingImageType >
      MovingOrientFilterType;
  typename MovingOrientFilterType::Pointer orientMoving
      = MovingOrientFilterType::New();
  orientMoving->UseImageDirectionOn();
  orientMoving->SetDesiredCoordinateOrientationToAxial();
  orientMoving->SetInput( movingImageReader->GetOutput() );
  orientMoving->Update();
  timeCollector.Stop( "Orient moving image" );

  timeCollector.Start( "Orient initial deformation field" );
  typedef itk::OrientImageFilter< VectorImageType, VectorImageType >
      VectorOrientFilterType;
  typename VectorOrientFilterType::Pointer orientInitField
      = VectorOrientFilterType::New();
  orientInitField->UseImageDirectionOn();
  orientInitField->SetDesiredCoordinateOrientationToAxial();
  orientInitField->SetInput( initField );
  orientInitField->Update();
  timeCollector.Stop( "Orient initial deformation field" );

  // Read the organ boundary
  if( organBoundaryFileName != "" )
    {
    timeCollector.Start( "Loading organ boundary" );
    // Do we have .vtk or .vtp models?
    std::string::size_type loc = organBoundaryFileName.find_last_of(".");
    if( loc == std::string::npos )
      {
      tube::ErrorMessage( "Failed to find an extension for organ boundary" );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    std::string extension = organBoundaryFileName.substr(loc);
    typename RegistrationType::BorderSurfacePointer borderSurface = NULL;
    if( extension == std::string(".vtk") )
      {
      vtkSmartPointer< vtkPolyDataReader > polyDataReader
          = vtkPolyDataReader::New();
      polyDataReader->SetFileName( organBoundaryFileName.c_str() );
      polyDataReader->Update();
      borderSurface = polyDataReader->GetOutput();
      }
    else if( extension == std::string(".vtp") )
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
    if( worldCoordinateSystem == "RAS" )
      {
      vtkSmartPointer< vtkTransform > RAStoLPS = vtkTransform::New();
      RAStoLPS->RotateZ(180); // flip in superior-inferior
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
    timeCollector.Stop( "Loading organ boundary" );
    }

  // Read normal vector image
  if( inputNormalVectorImageFileName != "" )
    {
    timeCollector.Start( "Loading normal vector image" );
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
    typename VectorOrientFilterType::Pointer orientNormals
        = VectorOrientFilterType::New();
    orientNormals->UseImageDirectionOn();
    orientNormals->SetDesiredCoordinateOrientationToAxial();
    orientNormals->SetInput( vectorImageReader->GetOutput() );
    orientNormals->Update();
    registrator->SetNormalVectorImage( orientNormals->GetOutput() );
    timeCollector.Stop( "Loading normal vector image" );
    }

  // Read weight image
  typedef itk::OrientImageFilter< WeightImageType, WeightImageType >
      WeightOrientFilterType;
  if( inputWeightImageFileName != "" )
    {
    timeCollector.Start( "Loading weight image" );
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
    typename WeightOrientFilterType::Pointer orientWeight
        = WeightOrientFilterType::New();
    orientWeight->UseImageDirectionOn();
    orientWeight->SetDesiredCoordinateOrientationToAxial();
    orientWeight->SetInput( weightImageReader->GetOutput() );
    orientWeight->Update();
    registrator->SetWeightImage( orientWeight->GetOutput() );
    timeCollector.Stop( "Loading weight image" );
    }

  // Error checking on lambda
  if( lambda > 0.0 )
    {
    tube::ErrorMessage( "Lambda must be negative." );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  // Report progress from reading input data
  double progress = 0.1;
  progressReporter.Report( progress );

  // Setup the registration
  timeCollector.Start( "Register" );
  registrator->SetFixedImage( orientFixed->GetOutput() );
  registrator->SetMovingImage( orientMoving->GetOutput() );
  registrator->SetInitialDeformationField( orientInitField->GetOutput() );
  registrator->SetNumberOfIterations( numberOfIterations );
  registrator->SetTimeStep( timeStep );
  registrator->SetComputeRegularizationTerm( !doNotPerformRegularization );
  registrator->SetUseAnisotropicRegularization(
      !doNotUseAnisotropicRegularization );
  registrator->SetLambda( lambda );

  // Watch the registration's progress
  tube::CLIFilterWatcher watchRegistration(registrator,
                                           "Anisotropic Diffusive Registration",
                                           CLPProcessInformation,
                                           0.8,
                                           progress );

  // Setup the warper: output parameters are derived from the fixed image
  // because the warped moving image should look like the fixed image
  // (ex. consider when the fixed image is outside of the extent of the moving
  // image - the transformed moving image must be in the space of the fixed
  // image)
  typedef itk::WarpImageFilter< MovingImageType,
                                MovingImageType,
                                VectorImageType > WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  typedef itk::LinearInterpolateImageFunction
      < MovingImageType, typename WarperType::CoordRepType > InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  warper->SetInput( movingImageReader->GetOutput() );
  warper->SetDeformationField( registrator->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputParametersFromImage( fixedImageReader->GetOutput() );
  warper->SetEdgePaddingValue( 0 );

  // Update triggers the registration and the warping
  warper->Update();
  timeCollector.Stop( "Register" );

  // Report progress from doing the registration
  progress = 0.9;
  progressReporter.Report( progress );

  // Reorient the registration's output deformation field: deformation field
  // is in the space of the fixed image
  typename VectorOrientFilterType::Pointer orientOutput
      = VectorOrientFilterType::New();
  orientOutput->UseImageDirectionOn();
  orientOutput->SetDesiredCoordinateDirection(
      fixedImageReader->GetOutput()->GetDirection() );
  orientOutput->SetInput( registrator->GetOutput() );

  // Write the deformation field
  if( outputDeformationFieldFileName != "" )
    {
    timeCollector.Start( "Write deformation field" );
    typedef itk::ImageFileWriter< VectorImageType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( outputDeformationFieldFileName );
    fieldWriter->SetInput( orientOutput->GetOutput() );
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
    timeCollector.Stop( "Write deformation field" );
    }

  // Output the transformation gridTransform (commented out for now, as not yet
  // supported in 3D Slicer)
//  if( outputTransformFileName != "" )
//    {
//    timeCollector.Start( "Write output transform" );
//    typedef itk::ImageFileWriter< VectorImageType > GridWriterType;
//    GridWriterType::Pointer gridWriter = GridWriterType::New();
//    gridWriter->SetFileName( outputTransformFileName );
//    gridWriter->SetInput( orientOutput->GetOutput() );
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
//    timeCollector.Stop( "Write output transform" );
//    }

  // Write the resampled moving image (in the space of the fixed image)
  if( outputResampledImageFileName != "" )
    {
    timeCollector.Start( "Write resampled moving image" );
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
    timeCollector.Stop( "Write resampled moving image" );
    }

  // Write the normal vector image: in the space of the fixed image
  if( outputNormalVectorImageFileName != "" )
    {
    timeCollector.Start( "Write normal vector image" );
    typename VectorOrientFilterType::Pointer orientNormals
        = VectorOrientFilterType::New();
    orientNormals->UseImageDirectionOn();
    orientNormals->SetDesiredCoordinateDirection(
        fixedImageReader->GetOutput()->GetDirection() );
    orientNormals->SetInput( registrator->GetNormalVectorImage() );
    orientNormals->Update();
    typedef itk::ImageFileWriter< VectorImageType > VectorWriterType;
    typename VectorWriterType::Pointer vectorWriter = VectorWriterType::New();
    vectorWriter->SetFileName( outputNormalVectorImageFileName );
    vectorWriter->SetInput( orientNormals->GetOutput() );
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
    timeCollector.Stop( "Write normal vector image" );
    }

  // Write the weight image: in the space of the fixed image
  if( outputWeightImageFileName != "" )
    {
    timeCollector.Start( "Write weight image" );
    typename WeightOrientFilterType::Pointer orientWeight
        = WeightOrientFilterType::New();
    orientWeight->UseImageDirectionOn();
    orientWeight->SetDesiredCoordinateDirection(
        fixedImageReader->GetOutput()->GetDirection() );
    orientWeight->SetInput( registrator->GetWeightImage() );
    orientWeight->Update();
    typedef itk::ImageFileWriter< WeightImageType > WeightWriterType;
    typename WeightWriterType::Pointer weightWriter = WeightWriterType::New();
    weightWriter->SetFileName( outputWeightImageFileName );
    weightWriter->SetInput( orientWeight->GetOutput() );
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
    timeCollector.Stop( "Write weight image" );
    }

  // Report progress from writing the outputs
  progress = 1.0;
  progressReporter.Report( progress );

  // Clean up, we're done
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
