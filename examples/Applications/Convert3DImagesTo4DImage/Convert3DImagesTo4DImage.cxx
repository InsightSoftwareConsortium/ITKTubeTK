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

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"

#include "tubeMessage.h"
#include <tubeWrite4DImageFrom3DImages.h>

#include <itkTimeProbesCollectorBase.h>

#include <itkExtractImageFilter.h>

#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <metaUtils.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "Convert3DImagesTo4DImageCLP.h"

// Your code should be within the DoIt function...
template <class TPixel>
int
DoIt(int argc, char * argv[])
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter("Convert3DImagesTo4DImage", CLPProcessInformation);
  progressReporter.Start();

  typedef TPixel                               InputPixelType;
  typedef itk::Image<InputPixelType, 3>        InputImageType;
  typedef itk::ImageFileReader<InputImageType> ReaderType;

  typedef tube::Write4DImageFrom3DImages<InputImageType> FilterType;

  typename FilterType::Pointer filter = FilterType::New();
  unsigned int                 num3DImages = argc - 2;
  filter->SetNumberOfInputImages(num3DImages);

  timeCollector.Start("Load data");
  double progress = 0.0;
  progressReporter.Report(progress);
  for (unsigned int i = 0; i < num3DImages; ++i)
  {
    typename ReaderType::Pointer reader = ReaderType::New();
    std::cout << "Reading #" << i << " : " << inputImageFileNames[i] << std::endl;
    reader->SetFileName(inputImageFileNames[i]);
    try
    {
      reader->Update();
    }
    catch (itk::ExceptionObject & err)
    {
      tube::ErrorMessage("Reading volume: Exception caught: " + std::string(err.GetDescription()));
      timeCollector.Report();
      return EXIT_FAILURE;
    }
    filter->SetNthInputImage(i, reader->GetOutput());

    progress = i * 0.75 / num3DImages;
    progressReporter.Report(progress);
  }
  timeCollector.Stop("Load data");

  filter->SetFileName(outputImageFileName);
  filter->Write();

  progress = 1.0;
  progressReporter.Report(progress);
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

namespace tube
{

// Get the component type and dimension of the image.
void
GetImageInformation(const std::string &                 fileName,
                    itk::ImageIOBase::IOComponentEnum & componentType,
                    unsigned int &                      dimension)
{
  typedef itk::ImageIOBase    ImageIOType;
  typedef itk::ImageIOFactory ImageIOFactoryType;

  ImageIOType::Pointer imageIO = ImageIOFactoryType::CreateImageIO(fileName.c_str(), itk::IOFileModeEnum::ReadMode);

  if (imageIO)
  {
    // Read the metadata from the image file.
    imageIO->SetFileName(fileName);
    imageIO->ReadImageInformation();

    componentType = imageIO->GetComponentType();
    dimension = imageIO->GetNumberOfDimensions();
  }
  else
  {
    tubeErrorMacro(<< "No ImageIO was found.");
  }
}

} // namespace tube

// Main
int
main(int argc, char * argv[])
{
  PARSE_ARGS;

  typedef itk::ImageIOBase             ImageIOType;
  typedef ImageIOType::IOComponentEnum IOComponentType;

  IOComponentType componentType = IOComponentType::UNKNOWNCOMPONENTTYPE;
  unsigned int    dimension = 0;
  try
  {
    tube::GetImageInformation(inputImageFileNames[0], componentType, dimension);

    if (dimension == 3)
    {
      switch (componentType)
      {
        case IOComponentType::UCHAR:
          return DoIt<unsigned char>(argc, argv);
        case IOComponentType::CHAR:
          return DoIt<char>(argc, argv);
        case IOComponentType::USHORT:
          return DoIt<unsigned short>(argc, argv);
        case IOComponentType::SHORT:
          return DoIt<short>(argc, argv);
        case IOComponentType::FLOAT:
          return DoIt<float>(argc, argv);
        case IOComponentType::DOUBLE:
          return DoIt<double>(argc, argv);
        case IOComponentType::INT:
          return DoIt<int>(argc, argv);
        case IOComponentType::UINT:
          return DoIt<unsigned int>(argc, argv);
        case IOComponentType::UNKNOWNCOMPONENTTYPE:
        default:
          tubeErrorMacro(<< "Unknown component type.");
          return EXIT_FAILURE;
      }
    }
    tubeErrorMacro(<< "Dimension size of " << dimension << " not supported");
    return EXIT_FAILURE;
  }
  catch (itk::ExceptionObject & ex)
  {
    tubeErrorMacro(<< "ITK exception caught. " << ex);
    return EXIT_FAILURE;
  }
  catch (...)
  {
    tubeErrorMacro(<< "Exception caught.");
    return EXIT_FAILURE;
  }

  return EXIT_FAILURE;
}
