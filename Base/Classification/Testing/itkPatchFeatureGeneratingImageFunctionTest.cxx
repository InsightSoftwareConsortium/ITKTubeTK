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

#include <vector>

#include <itkImage.h>
#include <itkFilterWatcher.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkPatchFeatureGeneratingImageFunction.h>

int itkPatchFeatureGeneratingImageFunctionTest(int argc, char* argv [] ) 
{
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << 
      " inputImage maskImage outputImage [meanHisto] [stdDevHisto]"
      << std::endl;
    return EXIT_FAILURE;
    }
  
  // Define the dimension of the images
  const unsigned int Dimension = 2;

  // Define the pixel type
  typedef float PixelType;
  
  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  // Declare the type for the Filter
  typedef itk::PatchFeatureGeneratingImageFunction< ImageType > FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }

  ImageType::Pointer inputImage = reader->GetOutput();

  FunctionType::Pointer func = FunctionType::New();
  func->SetInputImage( inputImage );
  func->SetWidth(20);

  itk::ImageRegionIteratorWithIndex< ImageType > outIter( inputImage,
    inputImage->GetLargestPossibleRegion() );
  ImageType::PointType pnt;

  // Evaluate
  outIter.GoToBegin();
  while( !outIter.IsAtEnd() )
    {
    inputImage->TransformIndexToPhysicalPoint( outIter.GetIndex(), pnt);
    std::vector<double> tf = func->Evaluate( pnt );
    if( outIter.GetIndex()[0] == outIter.GetIndex()[1] )
      {
      std::cout << "tf = " << " : " << tf[0] << " " << tf[1] << " " << tf[2] 
                << std::endl;
      }
    ++outIter;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
