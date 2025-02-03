/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "RegisterImagesCLP.h"

#include "itkImageToImageRegistrationHelper.h"

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[]);

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"


template <class TPixelType, unsigned int TDimension>
int
DoIt(int argc, char * argv[])
{

  PARSE_ARGS;

  enum VerboseLevelEnum
  {
    SILENT,
    STANDARD,
    VERBOSE
  };
  VerboseLevelEnum verbosity = SILENT;
  if (verbosityLevel == "Standard")
  {
    verbosity = STANDARD;
  }
  else if (verbosityLevel == "Verbose")
  {
    verbosity = VERBOSE;
  }

  typedef typename itk::Image<TPixelType, TDimension> ImageType;

  typedef typename itk::ImageToImageRegistrationHelper<ImageType> RegistrationType;

  typename RegistrationType::Pointer reger = RegistrationType::New();

  reger->SetReportProgress(true);

  if (verbosity >= STANDARD)
  {
    std::cout << "###Loading fixed image...";
  }
  reger->LoadFixedImage(fixedImage);
  if (verbosity >= STANDARD)
  {
    std::cout << "###DONE" << std::endl;
  }

  if (verbosity >= STANDARD)
  {
    std::cout << "###Loading moving image...";
  }
  reger->LoadMovingImage(movingImage);
  if (verbosity >= STANDARD)
  {
    std::cout << "###DONE" << std::endl;
  }

  if (loadTransform.size() > 1)
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Loading transform...";
    }
    reger->LoadTransform(loadTransform, invertLoadedTransform);
    if (verbosity >= STANDARD)
    {
      std::cout << "###DONE" << std::endl;
    }
  }

  if (skipInitialRandomSearch)
  {
    reger->SetUseEvolutionaryOptimization(false);
  }
  else
  {
    reger->SetUseEvolutionaryOptimization(true);
  }

  if (fixedLandmarks.size() > 0 || movingLandmarks.size() > 0)
  {
    if (initialization != "Landmarks")
    {
      std::cout << "WARNING: Landmarks specified, but initialization "
                << "process was not told to use landmarks. " << std::endl;
      std::cout << "Changing initialization to use landmarks." << std::endl;
    }
    if (verbosity >= STANDARD)
    {
      std::cout << "###Initialization: Landmarks" << std::endl;
    }
    initialization = "Landmarks";
    reger->SetInitialMethodEnum(RegistrationType::INIT_WITH_LANDMARKS);
    reger->SetFixedLandmarks(fixedLandmarks);
    reger->SetMovingLandmarks(movingLandmarks);
  }
  else if (initialization == "ImageCenters")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Initialization: ImageCenters" << std::endl;
    }
    reger->SetInitialMethodEnum(RegistrationType::INIT_WITH_IMAGE_CENTERS);
  }
  else if (initialization == "SecondMoments")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Initialization: SecondMoments" << std::endl;
    }
    reger->SetInitialMethodEnum(RegistrationType::INIT_WITH_SECOND_MOMENTS);
  }
  else if (initialization == "CentersOfMass")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Initialization: CentersOfMass" << std::endl;
    }
    reger->SetInitialMethodEnum(RegistrationType::INIT_WITH_CENTERS_OF_MASS);
  }
  else // if( initialization == "None" )
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Initialization: None" << std::endl;
    }
    reger->SetInitialMethodEnum(RegistrationType::INIT_WITH_NONE);
  }

  if (registration == "None")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: None" << std::endl;
    }
    reger->SetEnableInitialRegistration(false);
    reger->SetEnableRigidRegistration(false);
    reger->SetEnableAffineRegistration(false);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "Initial")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: Initial" << std::endl;
    }
    reger->SetEnableInitialRegistration(true);
    reger->SetEnableRigidRegistration(false);
    reger->SetEnableAffineRegistration(false);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "Rigid")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: Rigid" << std::endl;
    }
    reger->SetEnableInitialRegistration(false);
    reger->SetEnableRigidRegistration(true);
    reger->SetEnableAffineRegistration(false);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "Affine")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: Affine" << std::endl;
    }
    reger->SetEnableInitialRegistration(false);
    reger->SetEnableRigidRegistration(false);
    reger->SetEnableAffineRegistration(true);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "BSpline")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: BSpline" << std::endl;
    }
    reger->SetEnableInitialRegistration(false);
    reger->SetEnableRigidRegistration(false);
    reger->SetEnableAffineRegistration(false);
    reger->SetEnableBSplineRegistration(true);
  }
  else if (registration == "PipelineRigid")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: PipelineRigid" << std::endl;
    }
    reger->SetEnableInitialRegistration(true);
    reger->SetEnableRigidRegistration(true);
    reger->SetEnableAffineRegistration(false);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "PipelineAffine")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: PipelineAffine" << std::endl;
    }
    reger->SetEnableInitialRegistration(true);
    reger->SetEnableRigidRegistration(true);
    reger->SetEnableAffineRegistration(true);
    reger->SetEnableBSplineRegistration(false);
  }
  else if (registration == "PipelineBSpline")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Registration: PipelineBSpline" << std::endl;
    }
    reger->SetEnableInitialRegistration(true);
    reger->SetEnableRigidRegistration(true);
    reger->SetEnableAffineRegistration(true);
    reger->SetEnableBSplineRegistration(true);
  }

  if (metric == "NormCorr")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Metric: NormalizedCorrelation" << std::endl;
    }
    reger->SetRigidMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::NORMALIZED_CORRELATION_METRIC);
    reger->SetAffineMetricMethodEnum(
      RegistrationType ::OptimizedRegistrationMethodType ::NORMALIZED_CORRELATION_METRIC);
    reger->SetBSplineMetricMethodEnum(
      RegistrationType ::OptimizedRegistrationMethodType ::NORMALIZED_CORRELATION_METRIC);
  }
  else if (metric == "MeanSqrd")
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Metric: MeanSquared" << std::endl;
    }
    reger->SetRigidMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MEAN_SQUARED_ERROR_METRIC);
    reger->SetAffineMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MEAN_SQUARED_ERROR_METRIC);
    reger->SetBSplineMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MEAN_SQUARED_ERROR_METRIC);
  }
  else // if( metric == "MattesMI" )
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Metric: MattesMutualInformation" << std::endl;
    }
    reger->SetRigidMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MATTES_MI_METRIC);
    reger->SetAffineMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MATTES_MI_METRIC);
    reger->SetBSplineMetricMethodEnum(RegistrationType ::OptimizedRegistrationMethodType ::MATTES_MI_METRIC);
  }

  reger->SetSampleFromOverlap(sampleFromOverlap);
  if (verbosity >= STANDARD)
  {
    std::cout << "###sampleFromOverlap: " << sampleFromOverlap << std::endl;
  }

  typedef typename itk::ImageFileReader<itk::Image<unsigned char, TDimension>> ImageReader;
  typedef typename itk::ImageMaskSpatialObject<TDimension>                     ImageMaskSpatialObject;

  if (fixedImageMask != "")
  {
    reger->SetUseFixedImageMaskObject(true);

    typename ImageReader::Pointer reader = ImageReader::New();
    reader->SetFileName(fixedImageMask);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & exception)
    {
      std::cerr << "Exception caught while loading fixed image mask." << std::endl;
      std::cerr << exception << std::endl;
      return EXIT_FAILURE;
    }

    typename ImageMaskSpatialObject::Pointer mask = ImageMaskSpatialObject::New();
    mask->SetImage(reader->GetOutput());
    reger->SetFixedImageMaskObject(mask);

    if (verbosity >= STANDARD)
    {
      std::cout << "###useFixedImageMaskObject: true" << std::endl;
    }
  }
  else
  {
    reger->SetUseFixedImageMaskObject(false);
    if (verbosity >= STANDARD)
    {
      std::cout << "###useFixedImageMaskObject: false" << std::endl;
    }
  }

  // reger->SetSampleIntensityPortion( sampleIntensityPortion );
  // if ( verbosity >= STANDARD )
  //   {
  //   std::cout << "###sampleIntensityPortion: "
  //     << sampleIntensityPortion << std::endl;
  //   }
  // if( regionOfInterest.size() == 2*TDimension )
  //   {
  //   reger->SetRegionOfInterest( regionOfInterest );
  //   if ( verbosity >= STANDARD )
  //     {
  //     std::cout << "###regionOfInterest: ";
  //     std::cout << "    ###point1: ";
  //     for( unsigned int i=0; i< TDimension; i++ )
  //       {
  //       std::cout << regionOfInterest[i] << " ";
  //       }
  //     std::cout << "    ###point2: ";
  //     for( unsigned int i=0; i< TDimension; i++ )
  //       {
  //       std::cout << regionOfInterest[i+TDimension] << " ";
  //       }
  //     std::cout << std::endl;
  //     }
  //   }
  // else if( regionOfInterest.size() > 0 )
  //   {
  //   std::cerr <<
  //     "Error: region of interest does not contain two bounding points"
  //     << std::endl;
  //   return EXIT_FAILURE;
  //   }

  reger->SetMinimizeMemory(minimizeMemory);
  if (verbosity >= STANDARD)
  {
    std::cout << "###MinimizeMemory: " << minimizeMemory << std::endl;
  }

  reger->SetRandomNumberSeed(randomNumberSeed);

  reger->SetRigidMaxIterations(rigidMaxIterations);
  if (verbosity >= STANDARD)
  {
    std::cout << "###RigidMaxIterations: " << rigidMaxIterations << std::endl;
  }

  reger->SetAffineMaxIterations(affineMaxIterations);
  if (verbosity >= STANDARD)
  {
    std::cout << "###AffineMaxIterations: " << affineMaxIterations << std::endl;
  }

  reger->SetBSplineMaxIterations(bsplineMaxIterations);
  if (verbosity >= STANDARD)
  {
    std::cout << "###BSplineMaxIterations: " << bsplineMaxIterations << std::endl;
  }

  reger->SetRigidSamplingRatio(rigidSamplingRatio);
  if (verbosity >= STANDARD)
  {
    std::cout << "###RigidSamplingRatio: " << rigidSamplingRatio << std::endl;
  }
  reger->SetAffineSamplingRatio(affineSamplingRatio);
  if (verbosity >= STANDARD)
  {
    std::cout << "###AffineSamplingRatio: " << affineSamplingRatio << std::endl;
  }
  reger->SetBSplineSamplingRatio(bsplineSamplingRatio);
  if (verbosity >= STANDARD)
  {
    std::cout << "###BSplineSamplingRatio: " << bsplineSamplingRatio << std::endl;
  }

  /** not sure */
  if (interpolation == "NearestNeighbor")
  {
    reger->SetRigidInterpolationMethodEnum(
      RegistrationType ::OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION);
    reger->SetAffineInterpolationMethodEnum(
      RegistrationType ::OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION);
    reger->SetBSplineInterpolationMethodEnum(
      RegistrationType ::OptimizedRegistrationMethodType::NEAREST_NEIGHBOR_INTERPOLATION);
  }
  else if (interpolation == "Linear")
  {
    reger->SetRigidInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION);
    reger->SetAffineInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION);
    reger->SetBSplineInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION);
  }
  else if (interpolation == "BSpline")
  {
    reger->SetRigidInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION);
    reger->SetAffineInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION);
    reger->SetBSplineInterpolationMethodEnum(RegistrationType ::OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION);
  }
  if (verbosity >= STANDARD)
  {
    std::cout << "###RigidInterpolationMethod: " << interpolation << std::endl;
  }
  if (verbosity >= STANDARD)
  {
    std::cout << "###AffineInterpolationMethod: " << interpolation << std::endl;
  }
  if (verbosity >= STANDARD)
  {
    std::cout << "###BSplineInterpolationMethod: " << interpolation << std::endl;
  }

  reger->SetExpectedOffsetMagnitude(expectedOffset);
  if (verbosity >= STANDARD)
  {
    std::cout << "###ExpectedOffsetPixelMagnitude: " << expectedOffset << std::endl;
  }

  reger->SetExpectedRotationMagnitude(expectedRotation);
  if (verbosity >= STANDARD)
  {
    std::cout << "###ExpectedRotationMagnitude: " << expectedRotation << std::endl;
  }

  reger->SetExpectedScaleMagnitude(expectedScale);
  if (verbosity >= STANDARD)
  {
    std::cout << "###ExpectedScaleMagnitude: " << expectedScale << std::endl;
  }

  reger->SetExpectedSkewMagnitude(expectedSkew);
  if (verbosity >= STANDARD)
  {
    std::cout << "###ExpectedSkewMagnitude: " << expectedSkew << std::endl;
  }

  reger->SetBSplineControlPointPixelSpacing(controlPointSpacing);
  if (verbosity >= STANDARD)
  {
    std::cout << "###ExpectedBSplineControlPointPixelSpacing: " << controlPointSpacing << std::endl;
  }

  try
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Starting registration..." << std::endl;
    }
    reger->Update();
  }
  catch (itk::ExceptionObject & exception)
  {
    std::cerr << "Exception caught during helper class registration." << exception << std::endl;
    std::cerr << "Current Matrix Transform = " << std::endl;
    reger->GetCurrentMatrixTransform()->Print(std::cerr, 2);
    return EXIT_FAILURE;
  }
  catch (...)
  {
    std::cerr << "Uncaught exception during helper class registration." << std::endl;
    return EXIT_FAILURE;
  }

  if (resampledImage.size() > 1)
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###Resampling..." << std::endl;
    }
    typename ImageType::ConstPointer resultImage;
    try
    {
      if (interpolation == "NearestNeighbor")
      {
        resultImage =
          reger->ResampleImage(RegistrationType ::OptimizedRegistrationMethodType ::NEAREST_NEIGHBOR_INTERPOLATION,
                               NULL,
                               NULL,
                               NULL,
                               0,
                               resampledImagePortion);
      }
      else if (interpolation == "Linear")
      {
        resultImage = reger->ResampleImage(RegistrationType ::OptimizedRegistrationMethodType::LINEAR_INTERPOLATION,
                                           NULL,
                                           NULL,
                                           NULL,
                                           0,
                                           resampledImagePortion);
      }
      else if (interpolation == "BSpline")
      {
        resultImage = reger->ResampleImage(RegistrationType ::OptimizedRegistrationMethodType::BSPLINE_INTERPOLATION,
                                           NULL,
                                           NULL,
                                           NULL,
                                           0,
                                           resampledImagePortion);
      }
    }
    catch (itk::ExceptionObject & exception)
    {
      std::cerr << "Exception caught during helper class resampling." << exception << std::endl;
      std::cerr << "Current Matrix Transform = " << std::endl;
      reger->GetCurrentMatrixTransform()->Print(std::cerr, 2);
      return EXIT_FAILURE;
    }
    catch (...)
    {
      std::cerr << "Uncaught exception during helper class resampling." << std::endl;
      return EXIT_FAILURE;
    }

    try
    {
      reger->SaveImage(resampledImage, resultImage);
    }
    catch (itk::ExceptionObject & exception)
    {
      std::cerr << "Exception caught during helper class resampled image saving." << exception << std::endl;
      return EXIT_FAILURE;
    }
    catch (...)
    {
      std::cerr << "Uncaught exception during helper class resampled image saving." << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (saveTransform.size() > 1)
  {
    try
    {
      reger->SaveTransform(saveTransform);
    }
    catch (itk::ExceptionObject & exception)
    {
      std::cerr << "Exception caught during helper class transform saving." << exception << std::endl;
      return EXIT_FAILURE;
    }
    catch (...)
    {
      std::cerr << "Uncaught exception during helper class saving." << std::endl;
      return EXIT_FAILURE;
    }
  }

  if (saveDisplacementField.size() > 1)
  {
    try
    {
      reger->SaveDisplacementField(saveDisplacementField);
    }
    catch (itk::ExceptionObject & exception)
    {
      std::cerr << "Exception caught during helper class transform saving." << exception << std::endl;
      return EXIT_FAILURE;
    }
    catch (...)
    {
      std::cerr << "Uncaught exception during helper class saving." << std::endl;
      return EXIT_FAILURE;
    }
  }
  return EXIT_SUCCESS;
}

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  enum VerboseLevelEnum
  {
    SILENT,
    STANDARD,
    VERBOSE
  };
  VerboseLevelEnum verbosity = SILENT;
  if (verbosityLevel == "Standard")
  {
    verbosity = STANDARD;
  }
  else if (verbosityLevel == "Verbose")
  {
    verbosity = VERBOSE;
  }

  if (numberOfThreads != 0)
  {
    if (verbosity >= STANDARD)
    {
      std::cout << "###numberOfThreads: " << numberOfThreads << std::endl;
    }
    itk::MultiThreaderBase::SetGlobalDefaultNumberOfThreads(numberOfThreads);
  }

  tube::ParseArgsAndCallDoIt(fixedImage, argc, argv);

  return EXIT_SUCCESS;
}
