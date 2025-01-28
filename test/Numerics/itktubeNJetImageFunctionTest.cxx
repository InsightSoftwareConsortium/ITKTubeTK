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

#include "itktubeNJetImageFunction.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int
itktubeNJetImageFunctionTest(int argc, char * argv[])
{
  if (argc < 4)
  {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " function inputImage outputImage" << std::endl;
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
  typedef itk::tube::NJetImageFunction<ImageType> FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(argv[2]);
  try
  {
    reader->Update();
  }
  catch (itk::ExceptionObject & e)
  {
    std::cerr << "Exception caught during read:\n" << e;
    return EXIT_FAILURE;
  }

  ImageType::Pointer inputImage = reader->GetOutput();

  typedef itk::tube::NJetImageFunction<ImageType> FunctionType;
  FunctionType::Pointer                           func = FunctionType::New();
  func->SetInputImage(inputImage);

  ImageType::Pointer outputImage = ImageType::New();
  outputImage->CopyInformation(inputImage);
  outputImage->SetRegions(inputImage->GetLargestPossibleRegion());
  outputImage->Allocate();

  bool error = false;

  func->ComputeStatistics();
  func->SetUseProjection(false);
  double minI = 126;
  double maxI = 160;
  if (func->GetMin() != minI)
  {
    error = true;
    std::cerr << "Min = " << func->GetMin() << " != " << minI << std::endl;
  }
  if (func->GetMax() != maxI)
  {
    error = true;
    std::cerr << "Max = " << func->GetMax() << " != " << maxI << std::endl;
  }

  double scale = 1;

  FunctionType::VectorType v1, v2, d;
  v1.Fill(0);
  v1[0] = 1;
  v2.Fill(0);
  v2[1] = 1;
  double                   val;
  double                   val2;
  FunctionType::MatrixType h;
  FunctionType::MatrixType h2;

  int function = std::atoi(argv[1]);

  itk::ImageRegionIteratorWithIndex<ImageType> outIter(outputImage, outputImage->GetLargestPossibleRegion());
  ImageType::PointType                         pnt;
  outIter.GoToBegin();
  while (!outIter.IsAtEnd())
  {
    inputImage->TransformIndexToPhysicalPoint(outIter.GetIndex(), pnt);
    switch (function)
    {
      case 0:
      {
        outIter.Set(func->Evaluate(pnt, scale));
        break;
      }
      case 1:
      {
        outIter.Set(func->Evaluate(pnt, v1, scale));
        break;
      }
      case 2:
      {
        outIter.Set(func->Evaluate(pnt, v1, v2, scale));
        break;
      }
      case 3:
      {
        outIter.Set(func->EvaluateAtIndex(outIter.GetIndex(), scale));
        break;
      }
      case 4:
      {
        outIter.Set(func->EvaluateAtIndex(outIter.GetIndex(), v1, scale));
        break;
      }
      case 5:
      {
        outIter.Set(func->EvaluateAtIndex(outIter.GetIndex(), v1, v2, scale));
        break;
      }
      case 6:
      {
        func->Derivative(pnt, scale, d);
        outIter.Set(d[0]);
        break;
      }
      case 7:
      {
        func->Derivative(pnt, v1, scale, d);
        outIter.Set(d[0]);
        break;
      }
      case 8:
      {
        func->Derivative(pnt, v1, v2, scale, d);
        outIter.Set(d[0]);
        break;
      }
      case 9:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), scale, d);
        outIter.Set(d[1]);
        break;
      }
      case 10:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), v1, scale, d);
        outIter.Set(d[1]);
        break;
      }
      case 11:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), v1, v2, scale, d);
        outIter.Set(d[1]);
        break;
      }
      case 12:
      {
        func->Derivative(pnt, scale, d);
        val = func->GetMostRecentIntensity();
        outIter.Set(val);
        break;
      }
      case 13:
      {
        func->Derivative(pnt, v1, scale, d);
        val = func->GetMostRecentIntensity();
        outIter.Set(d[0]);
        break;
      }
      case 14:
      {
        func->Derivative(pnt, v1, v2, scale, d);
        val = func->GetMostRecentIntensity();
        outIter.Set(val);
        break;
      }
      case 15:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), scale, d);
        val = func->GetMostRecentIntensity();
        outIter.Set(d[1]);
        break;
      }
      case 16:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), v1, scale, d);
        val = func->GetMostRecentIntensity();
        outIter.Set(d[1]);
        break;
      }
      case 17:
      {
        func->DerivativeAtIndex(outIter.GetIndex(), v1, v2, scale, d);
        outIter.Set(d[1]);
        break;
      }
      case 18:
      {
        val = func->Jet(pnt, d, h, scale);
        outIter.Set(h[0][0]);
        break;
      }
      case 19:
      {
        val = func->JetAtIndex(outIter.GetIndex(), d, h, scale);
        outIter.Set(h[1][1]);
        break;
      }
      case 20:
      {
        val = func->Ridgeness(pnt, scale);
        outIter.Set(val);
        break;
      }
      case 21:
      {
        func->Ridgeness(pnt, scale);
        val = func->GetMostRecentRidgeLevelness();
        outIter.Set(val);
        break;
      }
      case 22:
      {
        func->Ridgeness(pnt, scale);
        val = func->GetMostRecentRidgeRoundness();
        outIter.Set(val);
        break;
      }
      case 23:
      {
        func->Ridgeness(pnt, scale);
        val = func->GetMostRecentRidgeCurvature();
        outIter.Set(val);
        break;
      }
      case 24:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), v1, v2, scale);
        outIter.Set(val);
        break;
      }
      case 25:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), v1, v2, scale);
        d = func->GetMostRecentDerivative();
        outIter.Set(val + d[1]);
        break;
      }
      case 26:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), v1, v2, scale);
        h = func->GetMostRecentHessian();
        outIter.Set(h[0][0]);
        break;
      }
      case 27:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), scale);
        val2 = func->GetMostRecentRidgeness();
        if (val != val2)
        {
          error = true;
          std::cerr << "Ridgeness and recent ridgeness differ at " << outIter.GetIndex() << std::endl;
          outIter.Set(val - val2);
        }
        else
        {
          outIter.Set(0);
        }
        break;
      }
      case 28:
      {
        func->RidgenessAtIndex(outIter.GetIndex(), scale);
        h = func->GetMostRecentHessian();
        val2 = func->HessianAtIndex(outIter.GetIndex(), scale, h2);
        val = std::fabs(h[0][0] - h2[0][0]) + std::fabs(h[0][1] - h2[0][1]) + std::fabs(h[1][0] - h2[1][0]) +
              std::fabs(h[1][1] - h2[1][1]);
        if (val > 0.0001)
        {
          error = true;
          std::cerr << "Ridgeness Hessian and Hessian differ at " << outIter.GetIndex() << std::endl;
          outIter.Set(val);
        }
        else
        {
          outIter.Set(0);
        }
        break;
      }
      case 29:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), scale);
        d = func->GetMostRecentDerivative();
        if (d.GetNorm() != 0)
        {
          d.Normalize();
        }
        else
        {
          d.Fill(0);
        }
        outIter.Set(d[0]);
        break;
      }
      case 30:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), v1, scale);
        d = func->GetMostRecentDerivative();
        if (d.GetNorm() != 0)
        {
          d.Normalize();
        }
        else
        {
          d.Fill(0);
        }
        outIter.Set(val);
        break;
      }
      case 31:
      {
        val = func->RidgenessAtIndex(outIter.GetIndex(), v1, v2, scale);
        d = func->GetMostRecentDerivative();
        if (d.GetNorm() != 0)
        {
          d.Normalize();
        }
        else
        {
          d.Fill(0);
        }
        outIter.Set(d[1]);
        break;
      }
      case 32:
      {
        func->Hessian(pnt, scale, h);
        outIter.Set(h[0][0]);
        break;
      }
      case 33:
      {
        func->Hessian(pnt, v1, scale, h);
        outIter.Set(h[0][0]);
        break;
      }
      case 34:
      {
        func->Hessian(pnt, v1, v2, scale, h);
        outIter.Set(h[0][0]);
        break;
      }
      case 35:
      {
        func->HessianAtIndex(outIter.GetIndex(), scale, h);
        outIter.Set(h[0][1]);
        break;
      }
      case 36:
      {
        func->Hessian(pnt, scale, h);
        outIter.Set(h[0][1]);
        break;
      }
      case 37:
      {
        func->HessianAtIndex(outIter.GetIndex(), v1, scale, h);
        outIter.Set(h[0][0] + h[1][1]);
        break;
      }
      case 38:
      {
        func->HessianAtIndex(outIter.GetIndex(), v1, v2, scale, h);
        outIter.Set(h[0][0] + h[1][1]);
        break;
      }
    }
    ++outIter;
  }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(argv[3]);
  writer->SetInput(outputImage);
  writer->SetUseCompression(true);

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
  if (!error)
  {
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}
