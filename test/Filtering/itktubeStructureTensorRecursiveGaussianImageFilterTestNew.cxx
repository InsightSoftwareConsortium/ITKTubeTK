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

#include "itktubeStructureTensorRecursiveGaussianImageFilter.h"

#include <itkImageDuplicator.h>
#include <itkImageFileWriter.h>

int
itktubeStructureTensorRecursiveGaussianImageFilterTestNew(int argc, char * argv[])
{
  // Define image
  enum
  {
    Dimension = 3
  };
  typedef double                           PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;

  // Set the value of sigma if specified in command line
  double sigma = 0.1;
  if (argc > 4)
  {
    sigma = std::atof(argv[4]);
  }

  double threshold = 0.05;
  if (argc > 5)
  {
    threshold = std::atof(argv[5]);
  }

  // Create test image I = f( x )
  ImageType::IndexType start;
  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  ImageType::SizeType size;
  size[0] = 100;
  size[1] = 100;
  size[2] = 100;

  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  ImageType::SpacingType spacing;
  spacing[0] = 0.1;
  spacing[1] = 0.1;
  spacing[2] = 0.1;

  ImageType::PointType origin;
  origin[0] = -5.0;
  origin[1] = -5.0;
  origin[2] = -5.0;

  ImageType::Pointer inImage = ImageType::New();
  inImage->SetRegions(region);
  inImage->SetSpacing(spacing);
  inImage->SetOrigin(origin);
  inImage->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> it(inImage, inImage->GetLargestPossibleRegion());

  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
  {
    ImageType::IndexType index = it.GetIndex();
    ImageType::PointType point;
    inImage->TransformIndexToPhysicalPoint(index, point);
    it.Set(static_cast<ImageType::PixelType>(std::sin(point[0]))); // sinx
  }

  typedef itk::ImageFileWriter<ImageType> WriterType;
  WriterType::Pointer                     writer = WriterType::New();
  writer->SetFileName("sin.mha");
  writer->SetInput(inImage);
  writer->Update();

  // Create a copy of the image
  typedef itk::ImageDuplicator<ImageType> ImageDuplicatorType;

  ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage(inImage);
  duplicator->Update();

  // Compute the outer ( diadic ) product of the gradient
  ImageType::Pointer                           prodImage = duplicator->GetOutput();
  itk::ImageRegionIteratorWithIndex<ImageType> pt(prodImage, prodImage->GetLargestPossibleRegion());

  for (pt.GoToBegin(); !pt.IsAtEnd(); ++pt)
  {
    ImageType::IndexType index = pt.GetIndex();
    ImageType::PointType point;
    prodImage->TransformIndexToPhysicalPoint(index, point);
    double cosValue = std::cos(point[0]);
    pt.Set(static_cast<ImageType::PixelType>(cosValue * cosValue)); // cosx*cosx
  }

  // Compute reference image as convolution of Gaussian and product of gradient
  typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType> GaussianFilterType;

  ImageType::Pointer refImage = prodImage;
  for (unsigned int i = 0; i < Dimension; i++)
  {
    GaussianFilterType::Pointer gaussian = GaussianFilterType::New();
    gaussian->SetInput(refImage);
    gaussian->SetDirection(i);
    gaussian->SetSigma(sigma);
    gaussian->Update();
    refImage = gaussian->GetOutput();
  }

  WriterType::Pointer writer1 = WriterType::New();
  writer1->SetFileName("cos2.mha");
  writer1->SetInput(refImage);
  writer1->Update();

  // Declare the type for the filter
  typedef itk::tube::StructureTensorRecursiveGaussianImageFilter<ImageType> StructureTensorFilterType;

  // Create a  Filter
  StructureTensorFilterType::Pointer filter = StructureTensorFilterType::New();

  // Connect the input images
  filter->SetInput(inImage);

  // set sigma
  filter->SetSigma(sigma);

  // Execute the filter
  filter->Update();

  // Define output type
  typedef StructureTensorFilterType::OutputImageType TensorImageType;
  typedef StructureTensorFilterType::OutputPixelType TensorImagePixelType;

  TensorImageType::Pointer outImage = filter->GetOutput();

  typedef itk::ImageFileWriter<TensorImageType> TensorWriterType;
  TensorWriterType::Pointer                     writer2 = TensorWriterType::New();
  writer2->SetFileName("tensor.nrrd");
  writer2->SetInput(outImage);
  writer2->Update();

  itk::ImageRegionIteratorWithIndex<TensorImageType> ot(outImage, outImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType>       rt(refImage, refImage->GetLargestPossibleRegion());

  for (ot.GoToBegin(), rt.GoToBegin(); !ot.IsAtEnd(); ++ot, ++rt)
  {
    TensorImageType::IndexType index = ot.GetIndex();

    bool boundaryPixel = false;
    for (unsigned int i = 0; i < 3; i++)
    {
      if (index[i] < 10 || index[i] >= 90)
      {
        boundaryPixel = true;
        break;
      }
    }
    if (boundaryPixel)
    {
      continue;
    }

    ImageType::PointType point;
    inImage->TransformIndexToPhysicalPoint(index, point);

    TensorImagePixelType tensor = ot.Get();

    for (unsigned int i = 0; i < Dimension; i++)
    {
      for (unsigned int j = i; j < Dimension; j++)
      {
        double refValue = 0.0;
        if (i == 0 && j == 0)
        {
          refValue = rt.Get();
        }
        if (tensor(i, j) - refValue > threshold || tensor(i, j) - refValue < -threshold)
        {
          std::cout << "Found mismatch at pixel of index[" << index[0] << " " << index[1] << " " << index[2]
                    << "] or point [" << point[0] << " " << point[1] << " " << point[2] << "]. Tensor element [" << i
                    << " " << j << "] is expected to be " << refValue << " but computed value is " << tensor(i, j)
                    << std::endl;
          return EXIT_FAILURE;
        }
      }
    }
  }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
