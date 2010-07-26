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
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include "../itkBlur3DImageFunction.h"

int itkBlur3DImageFunctionTest(int argc, char * argv[])
  {
  if(argc != 2)
    {
    std::cout << 
      "Usage: itkBlur3DImageFunctionTest <outputFilename>" << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::Image<float, 3>   ImageType;
  typedef ImageType::SizeType    ImageSizeType;
  typedef ImageType::SpacingType ImageSpacingType;

  ImageType::Pointer im = ImageType::New();
  ImageSizeType imSize;
  imSize[0] = 100;
  imSize[1] = 100;
  imSize[2] = 50;
  im->SetRegions(imSize);
  ImageSpacingType imSpacing;
  imSpacing[0] = 1;
  imSpacing[1] = 1;
  imSpacing[2] = 2;
  im->SetSpacing(imSpacing);

  im->Allocate();

  ImageType::IndexType index;
  index[0] = 50;
  index[1] = 50;
  index[2] = 25;
  im->SetPixel(index, 100);

  typedef itk::Blur3DImageFunction<ImageType> ImageOpType;
  ImageOpType::Pointer imOp = ImageOpType::New();

  imOp->SetInputImage(im);
  imOp->SetScale(3);

  ImageType::Pointer imOut = ImageType::New();
  imOut->SetRegions(imSize);
  imOut->SetSpacing(imSpacing);
  imOut->Allocate();

  itk::ImageRegionIteratorWithIndex<ImageType> itOut(imOut,
    imOut->GetLargestPossibleRegion());
  itOut.GoToBegin();
  while(!itOut.IsAtEnd())
    {
    itOut.Set(imOp->EvaluateAtIndex(itOut.GetIndex()));
    ++itOut;
    }

  typedef itk::ImageFileWriter<ImageType> ImageWriterType;
  ImageWriterType::Pointer imWriter = ImageWriterType::New();
  imWriter->SetFileName(argv[1]);
  imWriter->SetInput(imOut);
  imWriter->Update();

  return EXIT_SUCCESS;
  }

