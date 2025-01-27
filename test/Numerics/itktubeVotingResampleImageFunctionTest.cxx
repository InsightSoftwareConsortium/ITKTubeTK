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

#include "itktubeVotingResampleImageFunction.h"

#include <itkAffineTransform.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>

int
itktubeVotingResampleImageFunctionTest(int argc, char * argv[])
{
  if (argc != 4)
  {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " testNumber inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
  }

  const unsigned int                       Dimension = 2;
  typedef unsigned char                    PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;

  typedef ImageType::IndexType                                            ImageIndexType;
  typedef ImageType::SizeType                                             ImageSizeType;
  typedef double                                                          CoordRepType;
  typedef itk::AffineTransform<CoordRepType, Dimension>                   AffineTransformType;
  typedef itk::tube::VotingResampleImageFunction<ImageType, CoordRepType> InterpolatorType;

  typedef itk::ImageFileReader<ImageType> ImageReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageWriterType;


  // Create and configure an image
  ImageReaderType::Pointer imReader = ImageReaderType::New();
  imReader->SetFileName(argv[2]);
  imReader->Update();

  int testNumber = std::atoi(argv[1]);

  double        scale = 1;
  ImageSizeType size = imReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  if (testNumber == 0)
  {
    scale = 2;
    size[0] *= 2;
    size[1] *= 2;
  }
  else if (testNumber == 1)
  {
    scale = 1.5;
    size[0] *= 1.5;
    size[1] *= 1.5;
  }
  else if (testNumber == 2)
  {
    scale = 0.5;
    size[0] *= 0.5;
    size[1] *= 0.5;
  }
  // Create an affine transformation
  AffineTransformType::Pointer aff = AffineTransformType::New();
  aff->Scale(1.0 / scale);

  // Create a linear interpolation image function
  InterpolatorType::Pointer interp = InterpolatorType::New();
  interp->SetInputImage(imReader->GetOutput());

  // Create and configure a resampling filter
  itk::ResampleImageFilter<ImageType, ImageType>::Pointer resample;
  resample = itk::ResampleImageFilter<ImageType, ImageType>::New();
  resample->SetInput(imReader->GetOutput());
  resample->SetSize(size);
  resample->SetTransform(aff);
  resample->SetInterpolator(interp);

  ImageIndexType index;
  index.Fill(0);
  resample->SetOutputStartIndex(index);

  ImageType::PointType origin;
  origin.Fill(0.0);
  resample->SetOutputOrigin(origin);

  ImageType::SpacingType spacing;
  spacing.Fill(1.0);
  resample->SetOutputSpacing(spacing);

  // Run the resampling filter
  resample->Update();

  ImageWriterType::Pointer imWriter = ImageWriterType::New();
  imWriter->SetFileName(argv[3]);
  imWriter->SetInput(resample->GetOutput());
  imWriter->Update();

  return EXIT_SUCCESS;
}
