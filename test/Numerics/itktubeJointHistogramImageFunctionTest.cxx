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

#include "itktubeJointHistogramImageFunction.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int
itktubeJointHistogramImageFunctionTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " Image1 Image2 mask <1=linearize> outZImage [meanHist] [stdDevHist]" << std::endl;
    return EXIT_FAILURE;
  }

  // Define the dimension of the images
  enum
  {
    Dimension = 2
  };

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension> ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;

  // Declare the type for the Filter
  typedef itk::tube::JointHistogramImageFunction<ImageType> FunctionType;

  typedef itk::ImageFileWriter<FunctionType::HistogramType> HistoWriterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[1]);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during input read:\n" << e;
    return EXIT_FAILURE;
  }

  ImageType::Pointer inputImage1 = reader->GetOutput();

  // Read the mask
  reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during mask read:\n" << e;
    return EXIT_FAILURE;
  }
  ImageType::Pointer inputImage2 = reader->GetOutput();

  // Read the mask
  reader = ReaderType::New();
  reader->SetFileName(argv[3]);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during mask read:\n" << e;
    return EXIT_FAILURE;
  }
  ImageType::Pointer maskImage = reader->GetOutput();

  FunctionType::Pointer func = FunctionType::New();
  func->SetInputImage(inputImage1);
  func->SetInputMask(inputImage2);

  if (*argv[4] == '1')
  {
    func->SetForceDiagonalHistogram(true);
  }

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation(inputImage1);
  outputImage->SetRegions(inputImage1->GetLargestPossibleRegion());
  outputImage->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> outIter(outputImage, outputImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> maskIter(maskImage, maskImage->GetLargestPossibleRegion());

  // Precompute
  while (!outIter.IsAtEnd())
  {
    if (maskIter.Get() != 0)
    {
      func->PrecomputeAtIndex(outIter.GetIndex());
    }
    ++maskIter;
    ++outIter;
  }

  func->ComputeMeanAndStandardDeviation();

  if (argc > 4)
  {
    HistoWriterType::Pointer writerMean = HistoWriterType::New();
    writerMean->SetFileName(argv[6]);
    writerMean->SetUseCompression(true);
    writerMean->SetInput(func->GetMeanHistogram());
    try
    {
      writerMean->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught during mean write:\n" << e;
      return EXIT_FAILURE;
    }
  }

  if (argc > 5)
  {
    HistoWriterType::Pointer writerStdDev = HistoWriterType::New();
    writerStdDev->SetFileName(argv[7]);
    writerStdDev->SetUseCompression(true);
    writerStdDev->SetInput(func->GetStandardDeviationHistogram());
    try
    {
      writerStdDev->Update();
    }
    catch (itk::ExceptionObject & e)
    {
      std::cerr << "Exception caught during std. dev. write:\n" << e;
      return EXIT_FAILURE;
    }
  }

  // Evaluate
  outIter.GoToBegin();
  maskIter.GoToBegin();
  while (!outIter.IsAtEnd())
  {
    if (maskIter.Get() != 0)
    {
      double tf = func->EvaluateAtIndex(outIter.GetIndex());
      if (outIter.GetIndex()[0] == outIter.GetIndex()[1])
      {
        std::cout << "i = " << outIter.GetIndex()[0] << " : " << tf << std::endl;
      }
      outIter.Set(tf);
    }
    else
    {
      outIter.Set(0);
    }
    ++maskIter;
    ++outIter;
  }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[5]);
  writer->SetUseCompression(true);
  writer->SetInput(outputImage);
  try
  {
    writer->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during write:\n" << e;
    return EXIT_FAILURE;
  }


  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
