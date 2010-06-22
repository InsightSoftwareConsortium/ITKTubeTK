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


#include <itkImage.h>
#include <itkStructureTensorRecursiveGaussianImageFilter.h>
#include <itkSymmetricSecondRankTensor.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkSymmetricEigenAnalysisImageFilter.h>
#include <itkSymmetricEigenVectorAnalysisImageFilter.h>
#include <itkMatrix.h>
#include <itkVectorImage.h>
#include <itkVariableLengthVector.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMath.h>
#include <itkImageDuplicator.h>
#include <itkRecursiveGaussianImageFilter.h>

int itkStructureTensorRecursiveGaussianImageFilterTestNew(int argc, char* argv []  ) 
{
  // Define image
  const unsigned int Dimension = 3;
  typedef double PixelType;
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Set the value of sigma if specificed in command line
  double sigma = 0.1;
  if (argc > 4)
    {
    sigma = atof( argv[4] );
    }

  // Create test image I = f(x)
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
  
  itk::ImageRegionIteratorWithIndex<ImageType> it(inImage, 
                                           inImage->GetLargestPossibleRegion());
  
  for (it.GoToBegin(); !it.IsAtEnd(); ++it)
    {
    ImageType::IndexType index = it.GetIndex();
    ImageType::PointType point;
    inImage->TransformIndexToPhysicalPoint(index, point);
    it.Set(static_cast<ImageType::PixelType>(sin(point[0])));
    }

  typedef itk::ImageFileWriter<ImageType>        WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("sin.mha");
  writer->SetInput(inImage);
  writer->Update();
  
  // Create a copy of the image
  typedef itk::ImageDuplicator<ImageType>        ImageDuplicatorType;

  ImageDuplicatorType::Pointer duplicator = ImageDuplicatorType::New();
  duplicator->SetInputImage(inImage);
  duplicator->Update();

  // Blur the duplicated image
  typedef itk::RecursiveGaussianImageFilter<ImageType, ImageType>  
                                                 GaussianFilterType;
  
  ImageType::Pointer blurImage = duplicator->GetOutput();
  for (unsigned int i = 0; i < Dimension; i++)
    {
    GaussianFilterType::Pointer gaussian = GaussianFilterType::New();
    gaussian->SetInput(inImage);
    gaussian->SetDirection(i);
    gaussian->SetSigma(sigma);
    gaussian->Update();
    blurImage = gaussian->GetOutput();
    }
  
  // Declare the type for the filter
  typedef itk::StructureTensorRecursiveGaussianImageFilter<ImageType>  
                                                     StructureTensorFilterType;
            
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
  
  itk::ImageRegionIteratorWithIndex<TensorImageType> ot(outImage, 
                                           outImage->GetLargestPossibleRegion());
  itk::ImageRegionIteratorWithIndex<ImageType> bt(blurImage, 
                                           blurImage->GetLargestPossibleRegion());
  
  const double eps = 0.051;
  for (ot.GoToBegin(), bt.GoToBegin(); !ot.IsAtEnd(); ++ot, ++bt)
    {
    TensorImageType::IndexType index = ot.GetIndex();

    bool boundaryPixel = false;
    for (unsigned int i = 0; i < 3; i++)
      {
      if (index[i] < 20 || index[i] >= 80)
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

    // second x derivative of sin(x) is -sin(x), 
    // other second derivatives are all zero
    double minusSinX = -bt.Get();
    TensorImagePixelType upperMat = ot.Get();

    if (fabs(upperMat(0,0) - minusSinX) > eps)
      {
      std::cout << "The second X derivative of sin(x) image at index [" << 
        index[0] << " " << index[1] << " " << index[2] << "], point [" << 
        point[0] << " " << point[1] << " " << point[2] << "] is expected to be " 
        << minusSinX << " but the result is " << upperMat(0,0) << std::endl;
      }
    
    for (unsigned int i = 1; i < Dimension; i++)
      {
      if (upperMat(i, i) > eps || upperMat(i, i) < -eps)
        {
        std::cout << "Other than the second X derivative, other second derivatives "
          "of sin(x) image should be zero. This is violated at index [" << 
          index[0] << " " << index[1] << " " << index[2] << "], point [" << 
          point[0] << " " << point[1] << " " << point[2] << "]" << std::endl;
        }
      }
    for (unsigned int i = 0; i < Dimension-1; i++)
      {
      for (unsigned int j = i+1; j < Dimension; j++)
        {
        if (upperMat(i, i) > eps || upperMat(i, i) < -eps)
          {
          std::cout << "Other than the second X derivative, other second derivatives "
            "of sin(x) image should be zero. This is violated at index [" << 
            index[0] << " " << index[1] << " " << index[2] << "], point [" << 
            point[0] << " " << point[1] << " " << point[2] << "]" << std::endl;
          }
        }
      }
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}




