/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// TubeTKITK includes
#include "tubeEnhanceTubesUsingDiffusion.h"
#include "tubeMessage.h"

// TubeTK includes
#include "EnhanceTubesUsingDiffusionCLP.h"

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[]);

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template <class TPixel, unsigned int VDimension>
int
DoIt(int argc, char * argv[])
{
  PARSE_ARGS;

  typedef TPixel                                                  PixelType;
  typedef itk::Image<PixelType, VDimension>                       ImageType;
  typedef itk::ImageFileReader<ImageType>                         ReaderType;
  typedef tube::EnhanceTubesUsingDiffusion<PixelType, VDimension> FilterType;

  // read input image
  typename ReaderType::Pointer reader = ReaderType::New();

  try
  {
    reader->SetFileName(inputVolume.c_str());
    reader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    tube::ErrorMessage("Error reading input image: " + std::string(err.GetDescription()));
    return EXIT_FAILURE;
  }

  // compute vesselness
  typename FilterType::Pointer filter = FilterType::New();

  filter->SetInput(reader->GetOutput());

  filter->SetMinSigma(minSigma);
  filter->SetMaxSigma(maxSigma);
  filter->SetNumSigmaSteps(numSigmaSteps);
  filter->SetRecalculateTubeness(recalculateTubeness);
  filter->SetBeta(beta);
  filter->SetGamma(gamma);
  filter->SetEpsilon(epsilon);
  filter->SetOmega(omega);
  filter->SetSensitivity(sensitivity);
  filter->SetTimeStep(timeStep);
  filter->SetIterations(numIterations);

  filter->Update();

  // write output image
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  typename ImageWriterType::Pointer       writer = ImageWriterType::New();

  try
  {
    writer->SetFileName(outputVolume.c_str());
    writer->SetInput(filter->GetOutput());
    writer->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    tube::ErrorMessage("Error writing output vesselness image: " + std::string(err.GetDescription()));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt(inputVolume, argc, argv);
}
