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

#include "itktubeSpatialObjectToImageRegistrationHelper.h"
#include "itktubeSubSampleSpatialObjectFilter.h"

#include <itkImageFileReader.h>
#include <itkSpatialObjectReader.h>

#include <itkMemoryProbesCollectorBase.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>
#include <itkRecursiveGaussianImageFilter.h>
#include <itkTimeProbesCollectorBase.h>

#include <itkGradientDescentOptimizer.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>

//  The following section of code implements a Command observer
//  used to monitor the evolution of the registration process.
//
class CommandIterationUpdate : public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  itkNewMacro(Self);

protected:
  CommandIterationUpdate(void) {}
  ~CommandIterationUpdate(void)
  {
    if (measuresFileStream.is_open())
    {
      measuresFileStream.close();
    }
  }

  std::ofstream measuresFileStream;

public:
  typedef itk::OnePlusOneEvolutionaryOptimizer EvoOptimizerType;
  typedef itk::GradientDescentOptimizer        GradOptimizerType;

  void
  Execute(itk::Object * caller, const itk::EventObject & event) override
  {
    Execute((const itk::Object *)caller, event);
  }

  void
  Execute(const itk::Object * object, const itk::EventObject & event) override
  {
    if (!(itk::IterationEvent().CheckEvent(&event)))
    {
      return;
    }

    if (measuresFileStream.is_open())
    {
      typename EvoOptimizerType::ConstPointer  evoOptimizer = dynamic_cast<const EvoOptimizerType *>(object);
      typename GradOptimizerType::ConstPointer gradOptimizer = dynamic_cast<const GradOptimizerType *>(object);
      if (evoOptimizer.IsNotNull())
      {
        measuresFileStream << evoOptimizer->GetCurrentIteration() << ";";
        measuresFileStream << evoOptimizer->GetValue() << std::endl;
      }
      else if (gradOptimizer.IsNotNull())
      {
        measuresFileStream << gradOptimizer->GetCurrentIteration() << ";";
        measuresFileStream << gradOptimizer->GetValue() << std::endl;
      }
    }
  }

  void
  SetFileName(char * fileName)
  {
    measuresFileStream.open(fileName, std::ios::trunc);
  }

}; // End class CommandIterationUpdate

int
itktubeSpatialObjectToImageRegistrationPerformanceTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Missing Parameters: " << argv[0] << " Input_Image " << "Input_group " << "Output_Results"
              << std::endl;
    return EXIT_FAILURE;
  }

  typedef itk::SpatialObjectReader<3>                                     GroupReaderType;
  typedef itk::Image<double, 3>                                           ImageType;
  typedef itk::ImageFileReader<ImageType>                                 ImageReaderType;
  typedef itk::tube::SpatialObjectToImageRegistrationHelper<3, ImageType> RegistrationHelperType;

  // read image
  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(argv[1]);
  try
  {
    imageReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // Gaussian blur the original input image to increase the likelihood of vessel
  // spatial object overlapping with the vessel image at their initial alignment.
  // this enlarges the convergence zone.
  typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType> GaussianBlurFilterType;

  GaussianBlurFilterType::Pointer blurFilters[3];
  for (int i = 0; i < 3; i++)
  {
    blurFilters[i] = GaussianBlurFilterType::New();
    blurFilters[i]->SetSigma(3.0);
    blurFilters[i]->SetZeroOrder();
    blurFilters[i]->SetDirection(i);
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
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // read group
  GroupReaderType::Pointer groupReader = GroupReaderType::New();
  groupReader->SetFileName(argv[2]);
  try
  {
    groupReader->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // subsample points in group
  typedef itk::tube::SubSampleSpatialObjectFilter<> SubSampleFilterType;
  SubSampleFilterType::Pointer                      subSampleFilter = SubSampleFilterType::New();
  subSampleFilter->SetInput(groupReader->GetGroup());
  subSampleFilter->SetSampling(30);
  try
  {
    subSampleFilter->Update();
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // register the group and the image
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randGenerator =
    itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  randGenerator->Initialize(137593424);


  RegistrationHelperType::Pointer registrationHelper = RegistrationHelperType::New();

  registrationHelper->SetFixedImage(blurFilters[2]->GetOutput());
  registrationHelper->SetMovingSpatialObject(subSampleFilter->GetOutput());

  // Set Optimizer parameters.
  registrationHelper->SetExpectedOffsetMagnitude(3);
  registrationHelper->SetRigidMaxIterations(1000);


  // Add a time probe
  itk::TimeProbesCollectorBase   chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  char timeFile[250] = "";
  char valuesFile[250] = "";
  std::strcat(timeFile, argv[3]);
  std::strcat(valuesFile, argv[3]);
  std::strcat(timeFile, "TimeMemory.txt");
  std::strcat(valuesFile, "Values.txt");

  // Create stream to record the measure
  std::ofstream measuresFile;
  measuresFile.open(timeFile);
  if (!measuresFile.is_open())
  {
    std::cerr << "Unable to open: " << timeFile << std::endl;
    return EXIT_FAILURE;
  }

  try
  {
    memorymeter.Start("Registration");
    chronometer.Start("Registration");

    registrationHelper->Initialize();

    // Create the Command observer and register it with the optimizer.
    //
    CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
    observer->SetFileName(valuesFile);

    registrationHelper->SetObserver(observer);

    // Launch Registration
    registrationHelper->Update();

    chronometer.Stop("Registration");
    memorymeter.Stop("Registration");

    // Report the time and memory taken by the registration
    chronometer.Report(measuresFile);
    memorymeter.Report(measuresFile);
  }
  catch (itk::ExceptionObject & err)
  {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
  }

  measuresFile.close();
  return EXIT_SUCCESS;
}
