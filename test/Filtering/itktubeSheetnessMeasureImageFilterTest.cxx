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

#include "itktubeSheetnessMeasureImageFilter.h"

#include <itkHessianRecursiveGaussianImageFilter.h>
#include <itkImageFileWriter.h>

int itktubeSheetnessMeasureImageFilterTest( int ,char *[] )
{
  std::cout << "itktubeSheetnessMeasureImageFilterTest running..." << std::endl;

  // Define the dimension of the images
  const unsigned int myDimension = 3;

  // Declare the types of the images
  using myImageType = itk::Image<float, myDimension>;

  // Declare the type of the index to access images
  using myIndexType = itk::Index<myDimension>;

  // Declare the type of the size
  using mySizeType = itk::Size<myDimension>;

  // Declare the type of the Region
  using myRegionType = itk::ImageRegion<myDimension>;

  // Create the image
  myImageType::Pointer inputImage  = myImageType::New();

  std::cout << "Creating a synthetic image" << std::endl;

  // Define their size, and start index
  mySizeType size;
  size[0] = 8;
  size[1] = 8;
  size[2] = 8;

  myIndexType start;
  start.Fill( 0 );

  myRegionType region;
  region.SetIndex( start );
  region.SetSize( size );

  // Initialize Image A
  inputImage->SetLargestPossibleRegion( region );
  inputImage->SetBufferedRegion( region );
  inputImage->SetRequestedRegion( region );
  inputImage->Allocate();

  // Declare Iterator type for the input image
  using myIteratorType = itk::ImageRegionIteratorWithIndex<myImageType>;

  // Create one iterator for the Input Image A ( this is a light object )
  myIteratorType it( inputImage, inputImage->GetRequestedRegion() );

  // Initialize the content of Image A
  while( !it.IsAtEnd() )
    {
    it.Set( 0.0 );
    ++it;
    }

  size[0] = 4;
  size[1] = 4;
  size[2] = 2;

  start[0] = 2;
  start[1] = 2;
  start[2] = 2;

  // Create one iterator for an internal region
  region.SetSize( size );
  region.SetIndex( start );
  myIteratorType itb( inputImage, region );

  // Initialize the content the internal region
  while( !itb.IsAtEnd() )
    {
    itb.Set( 100.0 );
    ++itb;
    }

  std::cout << "Finished creating a synthetic image" << std::endl;

  std::cout << "Writing out the synthetic image" << std::endl;
  using InputImageWriterType = itk::ImageFileWriter<myImageType>;
  InputImageWriterType::Pointer inputImageWriter= InputImageWriterType::New();
  inputImageWriter->SetFileName( "SyntheticImageForSheetnessTest.mha" );
  inputImageWriter->SetInput( inputImage );
  inputImageWriter->Update();

  // Declare the type for the Hessian filter
  using myHessianFilterType = itk::HessianRecursiveGaussianImageFilter<
                                            myImageType >;

  // Declare the type for the sheetness measure filter
  using mySheetnessFilterType = itk::tube::SheetnessMeasureImageFilter< float >;

  using mySheetnessImageType = mySheetnessFilterType::OutputImageType;


  // Create a Hessian Filter
  myHessianFilterType::Pointer filterHessian = myHessianFilterType::New();

  // Create a sheetness Filter
  mySheetnessFilterType::Pointer filterSheetness = mySheetnessFilterType::New();


  // Connect the input images
  filterHessian->SetInput( inputImage );
  filterSheetness->SetInput( filterHessian->GetOutput() );

  // Select the value of Sigma
  filterHessian->SetSigma( 0.5 );


  // Execute the filter
  std::cout << "Generate sheetness measure" << std::endl;
  filterSheetness->Update();

  // Get the Smart Pointer to the Filter Output
  // It is important to do it AFTER the filter is Updated
  // Because the object connected to the output may be changed
  // by another during GenerateData() call
  mySheetnessImageType::Pointer outputImage = filterSheetness->GetOutput();

  //Write out the sheetness image
  //Define output type
  using SheetnessImageType = mySheetnessFilterType::OutputImageType;

  std::cout << "Write out the sheetness image" << std::endl;
  using SheetnessImageWriterType = itk::ImageFileWriter<SheetnessImageType>;
  SheetnessImageWriterType::Pointer writer= SheetnessImageWriterType::New();
  writer->SetFileName( "SheetnessImage.mha" );
  writer->SetInput( outputImage );
  writer->Update();

  return EXIT_SUCCESS;
}
