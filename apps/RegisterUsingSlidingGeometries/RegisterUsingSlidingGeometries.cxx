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

#include "itktubeAnisotropicDiffusiveRegistrationFilter.h"
#include "itktubeAnisotropicDiffusiveSparseRegistrationFilter.h"
#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMultiResolutionPDEDeformableRegistration.h>
#include <itkOrientImageFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkTransformFileReader.h>
#include <itkWarpImageFilter.h>

#include <vtkPolyDataReader.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkXMLPolyDataReader.h>

#include "RegisterUsingSlidingGeometriesCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

// Read and orient an image
template< class TImage >
bool ReadAndOrientImageAxial( TImage & outputImage, std::string fileName )
{
  typedef typename TImage::ObjectType ObjectType;

  typedef itk::ImageFileReader< ObjectType > FileReaderType;
  typename FileReaderType::Pointer imageReader = FileReaderType::New();
  imageReader->SetFileName( fileName.c_str() );
  try
    {
    imageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    return false;
    }

  typedef itk::OrientImageFilter< ObjectType, ObjectType > OrientFilterType;
  typename OrientFilterType::Pointer orient = OrientFilterType::New();
  orient->UseImageDirectionOn();
  orient->SetDesiredCoordinateOrientationToAxial();
  orient->SetInput( imageReader->GetOutput() );
  orient->Update();
  outputImage = orient->GetOutput();
  if( outputImage )
    {
    return true;
    }
  else
    {
    return false;
    }
}

// Reorient and write an image
template< class TImage >
bool ReorientAndWriteImage( TImage * inputImage,
                            typename TImage::DirectionType dir,
                            std::string fileName )
{
  typedef itk::OrientImageFilter< TImage, TImage > OrientFilterType;
  typename OrientFilterType::Pointer orient = OrientFilterType::New();
  orient->UseImageDirectionOn();
  orient->SetDesiredCoordinateDirection( dir );
  orient->SetInput( inputImage );
  orient->Update();

  typedef itk::ImageFileWriter< TImage > FileWriterType;
  typename FileWriterType::Pointer imageWriter = FileWriterType::New();
  imageWriter->SetFileName( fileName );
  imageWriter->SetUseCompression( true );
  imageWriter->SetInput( orient->GetOutput() );
  try
    {
    imageWriter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    return false;
    }
  return true;
}

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  // of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  bool reportProgress = false;
  tube::CLIProgressReporter progressReporter(
      "AnisotropicDiffusiveDeformableRegistration", CLPProcessInformation );
  if( reportProgress )
    {
    progressReporter.Start();
    }

  // Typedefs seeding definition of the anisotropic diffusive registration
  // filter
  const unsigned int                                      ImageDimension = 3;
  typedef TPixel                                          FixedPixelType;
  typedef TPixel                                          MovingPixelType;
  typedef itk::Image< FixedPixelType, ImageDimension >    FixedImageType;
  typedef itk::Image< MovingPixelType, ImageDimension >   MovingImageType;
  typedef double                                          VectorScalarType;
  typedef itk::Vector< VectorScalarType, ImageDimension > VectorType;
  typedef itk::Image< VectorType, ImageDimension >        VectorImageType;

  // Initialize the registration filter
  typedef itk::tube::DiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, VectorImageType >
      DiffusiveRegistrationFilterType;
  typedef itk::tube::AnisotropicDiffusiveRegistrationFilter
      < FixedImageType, MovingImageType, VectorImageType >
      AnisotropicDiffusiveRegistrationFilterType;
  typedef itk::tube::AnisotropicDiffusiveSparseRegistrationFilter
      < FixedImageType, MovingImageType, VectorImageType >
      AnisotropicDiffusiveSparseRegistrationFilterType;

  typename DiffusiveRegistrationFilterType::Pointer registrator = nullptr;
  typename AnisotropicDiffusiveRegistrationFilterType::Pointer
      anisotropicRegistrator = nullptr;
  typename AnisotropicDiffusiveSparseRegistrationFilterType::Pointer
      sparseAnisotropicRegistrator = nullptr;
  bool haveAnisotropicRegistrator = false;
  if( doNotPerformRegularization || doNotUseAnisotropicRegularization )
    {
    registrator = DiffusiveRegistrationFilterType::New();
    }
  else
    {
    haveAnisotropicRegistrator = true;
    if( anisotropicRegistrationType == "SlidingOrgan" )
      {
      registrator = AnisotropicDiffusiveRegistrationFilterType::New();
      anisotropicRegistrator
          = dynamic_cast< AnisotropicDiffusiveRegistrationFilterType * >(
              registrator.GetPointer() );
      }
    else if( anisotropicRegistrationType == "SparseSlidingOrgan" )
      {
      registrator = AnisotropicDiffusiveSparseRegistrationFilterType::New();
      sparseAnisotropicRegistrator
          = dynamic_cast< AnisotropicDiffusiveSparseRegistrationFilterType * >(
              registrator.GetPointer() );
      }
    else
      {
      tube::ErrorMessage( "Unknown registration type: "
                          + anisotropicRegistrationType );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    }

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
                        + std::string( err.GetDescription() ) );
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
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Loading moving image" );

  // Load the stopping criterion mask if provided
  if( stoppingCriterionMaskImageFileName != "" )
    {
    timeCollector.Start( "Loading stopping criterion mask" );
    typename DiffusiveRegistrationFilterType::StoppingCriterionMaskPointer
        stoppingCriterionMaskImage = nullptr;
    if( ReadAndOrientImageAxial( stoppingCriterionMaskImage,
                                 stoppingCriterionMaskImageFileName ) )
      {
      registrator->SetStoppingCriterionMask( stoppingCriterionMaskImage );
      }
    else
      {
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Loading stopping criterion mask" );
    }

  // Setup the initial deformation field: the initial deformation field should
  // be in the space of the fixed image
  timeCollector.Start( "Setup initial deformation field" );
  VectorImageType::Pointer initField = VectorImageType::New();

  // Preferably use an "initial transform" image, if given
  if( initialTransformImageFileName != "" )
    {
    typedef itk::ImageFileReader< VectorImageType > InitFieldReaderType;
    typename InitFieldReaderType::Pointer initFieldImageReader
        = InitFieldReaderType::New();
    initFieldImageReader->SetFileName( initialTransformImageFileName.c_str() );
    try
      {
      initFieldImageReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading initial transform image: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    initField = initFieldImageReader->GetOutput();
    }
  // If an "initial transform" transform is given, or if there is no
  // initial transform given
  else
    {
    // Setup the fixed image as a template
    typename FixedImageType::Pointer templateImage
        = fixedImageReader->GetOutput();
    initField->SetOrigin( templateImage->GetOrigin() );
    initField->SetSpacing( templateImage->GetSpacing() );
    initField->SetDirection( templateImage->GetDirection() );
    initField->SetLargestPossibleRegion(
      templateImage->GetLargestPossibleRegion() );
    initField->SetRequestedRegion( templateImage->GetRequestedRegion() );
    initField->SetBufferedRegion( templateImage->GetBufferedRegion() );
    initField->Allocate();

    // Use the "initial transform" transform if given
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
                            + std::string( err.GetDescription() ) );
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

        // For each voxel, find the displacement invoked by the given
        // initial transform.  This should work for all types of
        // transforms ( linear, nonlinear, B-spline, etc ) because
        // itk::Transform is the base for each.
        // Slicer saves transforms in LPS space ( see Slicer's
        // vtkMRMLTransformStorageNode::WriteData() ), so we don't
        // need to modify the initial transform according to the
        // worldCoordinateSystem variable ( unlike the surface model ).
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
    // If no initial transform is given, fill the initial field with
    // zero vectors
    else
      {
      VectorType zeroVector;
      zeroVector.Fill( 0.0 );
      initField->FillBuffer( zeroVector );
      }
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

  // Read the organ boundary if we are using the anisotropic regularizer
  if( haveAnisotropicRegistrator && organBoundaryFileName != "" )
    {
    timeCollector.Start( "Loading organ boundary" );
    // Do we have .vtk or .vtp models?
    std::string::size_type loc = organBoundaryFileName.find_last_of( "." );
    if( loc == std::string::npos )
      {
      tube::ErrorMessage( "Failed to find an extension for organ boundary" );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    std::string extension = organBoundaryFileName.substr( loc );
    typename AnisotropicDiffusiveRegistrationFilterType::BorderSurfacePointer
        borderSurface = NULL;
    if( extension == std::string( ".vtk" ) )
      {
      vtkSmartPointer< vtkPolyDataReader > polyDataReader
          = vtkSmartPointer< vtkPolyDataReader >::New();
      polyDataReader->SetFileName( organBoundaryFileName.c_str() );
      polyDataReader->Update();
      borderSurface = polyDataReader->GetOutput();
      }
    else if( extension == std::string( ".vtp" ) )
      {
      vtkSmartPointer< vtkXMLPolyDataReader > polyDataReader
          = vtkSmartPointer< vtkXMLPolyDataReader >::New();
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

    // If the world coordinate system is RAS ( i.e. called from 3D Slicer )
    // then the model will be in RAS space while the images will be in LPS
    // space.  It's easiest to transform the model to LPS space.
    if( worldCoordinateSystem == "RAS" )
      {
      vtkSmartPointer< vtkTransform > RAStoLPS =
        vtkSmartPointer< vtkTransform >::New();
      RAStoLPS->RotateZ( 180 ); // flip in superior-inferior
      vtkSmartPointer< vtkTransformPolyDataFilter > transformPolyDataFilter
          = vtkSmartPointer< vtkTransformPolyDataFilter >::New();
#if VTK_MAJOR_VERSION > 5
      transformPolyDataFilter->SetInputData( borderSurface );
#else
      transformPolyDataFilter->SetInput( borderSurface );
#endif
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
    if( anisotropicRegistrator )
      {
      anisotropicRegistrator->SetBorderSurface( borderSurface );
      }
    else if( sparseAnisotropicRegistrator )
      {
      sparseAnisotropicRegistrator->SetBorderSurface( borderSurface );
      }
    timeCollector.Stop( "Loading organ boundary" );
    }

  typename AnisotropicDiffusiveSparseRegistrationFilterType::TubeListPointer
    tubeList = nullptr;

  // Read tube spatial object if we are using the sparse anisotropic regularizer
  if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
    {
    timeCollector.Start( "Loading tube list" );
    typedef itk::SpatialObjectReader< ImageDimension > TubeReaderType;
    TubeReaderType::Pointer tubeReader = TubeReaderType::New();
    tubeReader->SetFileName( tubeSpatialObjectFileName.c_str() );
    try
      {
      tubeReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading tube list: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    typedef itk::GroupSpatialObject< ImageDimension > GroupType;
    typename GroupType::Pointer group = tubeReader->GetGroup();
    tubeList = group->GetChildren();
    sparseAnisotropicRegistrator->SetTubeList( tubeList );
    timeCollector.Stop( "Loading tube list" );
    }

  // Read normal vector image if we are using the anisotropic regularizer
  if( haveAnisotropicRegistrator && inputNormalVectorImageFileName != "" )
    {
    timeCollector.Start( "Loading normal vector image" );
    if( anisotropicRegistrator )
      {
      typename AnisotropicDiffusiveRegistrationFilterType::NormalVectorImageType
          ::Pointer normalImage = nullptr;
      if( ReadAndOrientImageAxial( normalImage,
                                   inputNormalVectorImageFileName ) )
        {
        anisotropicRegistrator->SetNormalVectorImage( normalImage );
        }
      else
        {
        timeCollector.Report();
        return EXIT_FAILURE;
        }
      }
    else if( sparseAnisotropicRegistrator )
      {
      typename AnisotropicDiffusiveSparseRegistrationFilterType
          ::NormalMatrixImageType::Pointer normalImage = nullptr;
      if( ReadAndOrientImageAxial( normalImage,
                                   inputNormalVectorImageFileName ) )
        {
        sparseAnisotropicRegistrator->SetNormalMatrixImage( normalImage );
        }
      else
        {
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        return EXIT_FAILURE;
        }
      }
    timeCollector.Stop( "Loading normal vector image" );
    }

  // Read weight regularizations image if we are using the anisotropic
  // registrator
  if( haveAnisotropicRegistrator
      && inputWeightRegularizationsImageFileName != "" )
    {
    timeCollector.Start( "Loading weight regularizations image" );
    if( anisotropicRegistrator )
      {
      typename AnisotropicDiffusiveRegistrationFilterType::WeightImageType
          ::Pointer weightImage = nullptr;
      if( ReadAndOrientImageAxial( weightImage,
                                   inputWeightRegularizationsImageFileName ) )
        {
        anisotropicRegistrator->SetWeightImage( weightImage );
        }
      else
        {
        timeCollector.Report();
        return EXIT_FAILURE;
        }
      }
    else if( sparseAnisotropicRegistrator )
      {
      typename AnisotropicDiffusiveSparseRegistrationFilterType
          ::WeightComponentImageType::Pointer weightImage = nullptr;
      if( ReadAndOrientImageAxial( weightImage,
                                   inputWeightRegularizationsImageFileName ) )
        {
        sparseAnisotropicRegistrator->SetWeightRegularizationsImage(
            weightImage );
        }
      else
        {
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        return EXIT_FAILURE;
        }
      }
    timeCollector.Stop( "Loading weight regularizations image" );
    }

  // Read the weight structures image ( in the space of the fixed image ) if we
  // are using the sparse anisotropic regularization
  if( sparseAnisotropicRegistrator && inputWeightStructuresImageFileName != "" )
    {
    timeCollector.Start( "Loading weight structures image" );
    typename AnisotropicDiffusiveSparseRegistrationFilterType
        ::WeightMatrixImageType::Pointer weightImage = nullptr;
    if( ReadAndOrientImageAxial( weightImage,
                                 inputWeightStructuresImageFileName ) )
      {
      sparseAnisotropicRegistrator->SetWeightStructuresImage( weightImage );
      }
    else
      {
      timeCollector.Report();
      if( tubeSpatialObjectFileName != "" )
        {
        delete tubeList;
        }
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Loading weight structures image" );
    }

  // Error checking on lambda
  if( lambda < 0.0 )
    {
    tube::ErrorMessage( "Lambda must be positive." );
    timeCollector.Report();
    if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
      {
      delete tubeList;
      }
    return EXIT_FAILURE;
    }

  // Error checking on gamma
  if( gamma < 0.0 && gamma != -1.0 )
    {
    tube::ErrorMessage( "Gamma must be positive." );
    timeCollector.Report();
    if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
      {
      delete tubeList;
      }
    return EXIT_FAILURE;
    }

  // Error checking on number of iterations
  if( numberOfIterations.size() <= 0 )
    {
    tube::ErrorMessage( "You must provide a list of number of iterations." );
    timeCollector.Report();
    if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
      {
      delete tubeList;
      }
    return EXIT_FAILURE;
    }

  // Error checking on regularization weightings
  if( regularizationWeightings.size() <= 0 )
    {
    tube::ErrorMessage( "You must provide regularization weightings." );
    timeCollector.Report();
    if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
      {
      delete tubeList;
      }
    return EXIT_FAILURE;
    }

  // Report progress from reading input data
  double progress = 0.1;
  if( reportProgress )
    {
    progressReporter.Report( progress );
    }

  // Start the registration
  timeCollector.Start( "Register" );

  // Setup the anisotropic registrator
  registrator->SetTimeStep( timeStep );
  registrator->SetComputeRegularizationTerm( !doNotPerformRegularization );
  registrator->SetComputeIntensityDistanceTerm(
    !doNotComputeIntensityDistanceTerm );
  if( anisotropicRegistrator )
    {
    anisotropicRegistrator->SetLambda( lambda );
    anisotropicRegistrator->SetGamma( gamma );
    }
  if( sparseAnisotropicRegistrator )
    {
    sparseAnisotropicRegistrator->SetLambda( lambda );
    sparseAnisotropicRegistrator->SetGamma( gamma );
    }
  registrator->SetMaximumRMSError( maximumRMSError );
  registrator->SetRegularizationWeightings( regularizationWeightings );
  registrator->SetBackgroundIntensity( backgroundIntensity );
  registrator->SetStoppingCriterionEvaluationPeriod(
    static_cast<unsigned int>( stoppingCriterionPeriod ) );
  registrator->SetStoppingCriterionMaxTotalEnergyChange(
    maximumTotalEnergyChange );

  // Setup the multiresolution PDE filter - we use the recursive pyramid because
  // we don't want the deformation field to undergo Gaussian smoothing on the
  // last iteration, which undermines the efforts we are making with anisotropic
  // diffusion

  // Setup the levels, iterations and max error of Gaussian kernel
  int numberOfLevels = numberOfIterations.size();
  unsigned int * iterations = new unsigned int [ numberOfLevels ];
  std::copy( numberOfIterations.begin(), numberOfIterations.end(), iterations );
  // This is the maximum error for the multiresolution pyramid smoother, not
  // the registration filter
  double maximumError = 0.01;

  // Setup the multiresolution pyramids
  typedef TPixel MultiResolutionRealType;
  typedef itk::Image< MultiResolutionRealType, ImageDimension >
      MultiResolutionRealImageType;
  typedef itk::RecursiveMultiResolutionPyramidImageFilter
      < FixedImageType, MultiResolutionRealImageType > FixedImagePyramidType;
  typename FixedImagePyramidType::Pointer fixedImagePyramid =
      FixedImagePyramidType::New();
  fixedImagePyramid->SetNumberOfLevels( numberOfLevels );
  fixedImagePyramid->SetMaximumError( maximumError );
  fixedImagePyramid->UseShrinkImageFilterOff();

  typedef itk::RecursiveMultiResolutionPyramidImageFilter
      < MovingImageType, MultiResolutionRealImageType > MovingImagePyramidType;
  typename MovingImagePyramidType::Pointer movingImagePyramid
      = MovingImagePyramidType::New();
  movingImagePyramid->SetNumberOfLevels( numberOfLevels );
  movingImagePyramid->SetMaximumError( maximumError );
  movingImagePyramid->UseShrinkImageFilterOff();

  // Setup the diffusion filter to work with multiresolution registration
  registrator->SetHighResolutionTemplate( orientFixed->GetOutput() );

  // Setup the multiresolution registrator
  typedef itk::MultiResolutionPDEDeformableRegistration
      < FixedImageType,
      MovingImageType,
      VectorImageType,
      MultiResolutionRealType > MultiResolutionRegistrationFilterType;
  typename MultiResolutionRegistrationFilterType::Pointer multires
      = MultiResolutionRegistrationFilterType::New();
  multires->SetRegistrationFilter( registrator );
  multires->SetFixedImage( orientFixed->GetOutput() );
  multires->SetMovingImage( orientMoving->GetOutput() );
  multires->SetArbitraryInitialDisplacementField(
    orientInitField->GetOutput() );
  multires->SetFixedImagePyramid( fixedImagePyramid );
  multires->SetMovingImagePyramid( movingImagePyramid );
  multires->SetNumberOfLevels( numberOfLevels );
  multires->SetNumberOfIterations( iterations );

  // Watch the registration's progress
  if( reportProgress )
    {
    tube::CLIFilterWatcher watchRegistration(
          registrator, "Anisotropic Diffusive Registration",
          CLPProcessInformation, 0.8, progress );
    }

  // Setup the warper: output parameters are derived from the fixed image
  // because the warped moving image should look like the fixed image
  // ( ex. consider when the fixed image is outside of the extent of the moving
  // image - the transformed moving image must be in the space of the fixed
  // image )
  typedef itk::WarpImageFilter< MovingImageType,
                                MovingImageType,
                                VectorImageType > WarperType;
  typename WarperType::Pointer warper = WarperType::New();
  typedef itk::LinearInterpolateImageFunction
      < MovingImageType, typename WarperType::CoordRepType > InterpolatorType;
  typename InterpolatorType::Pointer interpolator = InterpolatorType::New();
  warper->SetInput( movingImageReader->GetOutput() );
  warper->SetDisplacementField( multires->GetOutput() );
  warper->SetInterpolator( interpolator );
  warper->SetOutputParametersFromImage( fixedImageReader->GetOutput() );
  warper->SetEdgePaddingValue( backgroundIntensity );

  // Update triggers the registration and the warping
  warper->Update();
  timeCollector.Stop( "Register" );

  // Report progress from doing the registration
  progress = 0.9;
  if( reportProgress )
    {
    progressReporter.Report( progress );
    }

  // Reorient the registration's output deformation field: deformation field
  // is in the space of the fixed image
  typename VectorOrientFilterType::Pointer orientOutput
      = VectorOrientFilterType::New();
  orientOutput->UseImageDirectionOn();
  orientOutput->SetDesiredCoordinateDirection(
      fixedImageReader->GetOutput()->GetDirection() );
  orientOutput->SetInput( multires->GetOutput() );

  // Write the deformation field
  if( outputDeformationFieldFileName != "" )
    {
    timeCollector.Start( "Write deformation field" );
    typedef itk::ImageFileWriter< VectorImageType > FieldWriterType;
    typename FieldWriterType::Pointer fieldWriter = FieldWriterType::New();
    fieldWriter->SetFileName( outputDeformationFieldFileName );
    fieldWriter->SetUseCompression( true );
    fieldWriter->SetInput( orientOutput->GetOutput() );
    try
      {
      fieldWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
        {
        delete tubeList;
        }
      delete [] iterations;
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Write deformation field" );
    }

  // Write the resampled moving image ( in the space of the fixed image )
  if( outputResampledImageFileName != "" )
    {
    timeCollector.Start( "Write resampled moving image" );
    typedef itk::ImageFileWriter< MovingImageType > ImageWriterType;
    typename ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetFileName( outputResampledImageFileName );
    imageWriter->SetUseCompression( true );
    imageWriter->SetInput( warper->GetOutput() );
    try
      {
      imageWriter->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Writing volume: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
        {
        delete tubeList;
        }
      delete [] iterations;
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Write resampled moving image" );
    }

  // Write the normal vector image ( in the space of the fixed image ) if we are
  // using the anisotropic regularization
  if( haveAnisotropicRegistrator && outputNormalVectorImageFileName != "" )
    {
    timeCollector.Start( "Write normal vector image" );
    if( anisotropicRegistrator )
      {
      if( anisotropicRegistrator->GetHighResolutionNormalVectorImage() )
        {
        if( !ReorientAndWriteImage(
            anisotropicRegistrator->GetHighResolutionNormalVectorImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputNormalVectorImageFileName ) )
          {
          timeCollector.Report();
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      else
        {
        if( !ReorientAndWriteImage(
            anisotropicRegistrator->GetNormalVectorImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputNormalVectorImageFileName ) )
          {
          timeCollector.Report();
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      }
    else if( sparseAnisotropicRegistrator )
      {
      // We are giving back the 0th normal vector image, since Slicer cannot
      // accept the normal matrix image
      bool haveHighRes = false;
      if( sparseAnisotropicRegistrator->GetHighResolutionNormalMatrixImage() )
        {
        haveHighRes = true;
        }

      // Get the extension for the normal matrix
      std::string::size_type loc
          = outputNormalVectorImageFileName.find_last_of( "." );
      if( loc == std::string::npos )
        {
        tube::ErrorMessage( "Failed to find an extension for normal matrix" );
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        delete [] iterations;
        return EXIT_FAILURE;
        }

      std::string base = outputNormalVectorImageFileName.substr( 0, loc );
      std::string extension = outputNormalVectorImageFileName.substr( loc );
      std::string outputFileName;
      std::stringstream out;

      // Write out the normal matrix image
      out.clear();
      out.str( "" );
      out << base << "Matrix" << extension;
      outputFileName = out.str();
      if( !ReorientAndWriteImage(
          sparseAnisotropicRegistrator->GetHighResolutionNormalMatrixImage(),
          fixedImageReader->GetOutput()->GetDirection(),
          outputFileName ) )
        {
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        delete [] iterations;
        return EXIT_FAILURE;
        }

      // Write out the other normal vector images
      typedef typename AnisotropicDiffusiveSparseRegistrationFilterType
          ::NormalVectorImageType NormalVectorImageType;
      typedef typename AnisotropicDiffusiveSparseRegistrationFilterType
          ::NormalVectorImagePointer NormalVectorImagePointer;
      NormalVectorImagePointer normalImage = NormalVectorImageType::New();
      for( unsigned int i = 0; i < ImageDimension; i++ )
        {
        sparseAnisotropicRegistrator->GetHighResolutionNormalVectorImage(
            normalImage, i, haveHighRes );
        out.clear();
        out.str( "" );
        out << base << i << extension;
        if( i == 0 )
          {
          outputFileName = outputNormalVectorImageFileName;
          }
        else
          {
          outputFileName = out.str();
          }
        if( !ReorientAndWriteImage(
            normalImage.GetPointer(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputFileName ) )
          {
          timeCollector.Report();
          if( tubeSpatialObjectFileName != "" )
            {
            delete tubeList;
            }
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      }
    timeCollector.Stop( "Write normal vector image" );
    }

  // Write the weight regularizations image ( in the space of the fixed
  // image ) if we are using the anisotropic regularization
  if( haveAnisotropicRegistrator
      && outputWeightRegularizationsImageFileName != "" )
    {
    timeCollector.Start( "Write weight regularizations image" );
    if( anisotropicRegistrator )
      {
      if( anisotropicRegistrator->GetHighResolutionWeightImage() )
        {
        if( !ReorientAndWriteImage(
            anisotropicRegistrator->GetHighResolutionWeightImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputWeightRegularizationsImageFileName ) )
          {
          timeCollector.Report();
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      else
        {
        if( !ReorientAndWriteImage(
            anisotropicRegistrator->GetWeightImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputWeightRegularizationsImageFileName ) )
          {
          timeCollector.Report();
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      }
    else if( sparseAnisotropicRegistrator )
      {
      if( sparseAnisotropicRegistrator
            ->GetHighResolutionWeightRegularizationsImage() )
        {
        if( !ReorientAndWriteImage(
            sparseAnisotropicRegistrator
              ->GetHighResolutionWeightRegularizationsImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputWeightRegularizationsImageFileName ) )
          {
          timeCollector.Report();
          if( tubeSpatialObjectFileName != "" )
            {
            delete tubeList;
            }
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      else
        {
        if( !ReorientAndWriteImage(
            sparseAnisotropicRegistrator->GetWeightRegularizationsImage(),
            fixedImageReader->GetOutput()->GetDirection(),
            outputWeightRegularizationsImageFileName ) )
          {
          timeCollector.Report();
          if( tubeSpatialObjectFileName != "" )
            {
            delete tubeList;
            }
          delete [] iterations;
          return EXIT_FAILURE;
          }
        }
      }
    timeCollector.Stop( "Write weight regularizations image" );
    }

  // Write the weight matrix image ( in the space of the fixed image ) if we are
  // using the sparse anisotropic regularization
  if( sparseAnisotropicRegistrator
      && outputWeightStructuresImageFileName != "" )
    {
    timeCollector.Start( "Write weight structures image" );
    if( sparseAnisotropicRegistrator->GetHighResolutionWeightStructuresImage() )
      {
      if( !ReorientAndWriteImage(
          sparseAnisotropicRegistrator
            ->GetHighResolutionWeightStructuresImage(),
          fixedImageReader->GetOutput()->GetDirection(),
          outputWeightStructuresImageFileName ) )
        {
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        delete [] iterations;
        return EXIT_FAILURE;
        }
      }
    else
      {
      if( !ReorientAndWriteImage(
          sparseAnisotropicRegistrator->GetWeightStructuresImage(),
          fixedImageReader->GetOutput()->GetDirection(),
          outputWeightStructuresImageFileName ) )
        {
        timeCollector.Report();
        if( tubeSpatialObjectFileName != "" )
          {
          delete tubeList;
          }
        delete [] iterations;
        return EXIT_FAILURE;
        }
      }
    timeCollector.Stop( "Write weight structures image" );
    }

  // Report progress from writing the outputs
  progress = 1.0;
  if( reportProgress )
    {
    progressReporter.Report( progress );
    }

  // Clean up, we're done
  delete [] iterations;

  if( sparseAnisotropicRegistrator && tubeSpatialObjectFileName != "" )
    {
    delete tubeList;
    }

  if( reportProgress )
    {
    progressReporter.End();
    }
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( fixedImageFileName, argc, argv );
}
