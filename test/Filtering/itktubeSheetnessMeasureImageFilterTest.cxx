/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

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
  typedef itk::Image<float, myDimension>           myImageType;

  // Declare the type of the index to access images
  typedef itk::Index<myDimension>             myIndexType;

  // Declare the type of the size
  typedef itk::Size<myDimension>              mySizeType;

  // Declare the type of the Region
  typedef itk::ImageRegion<myDimension>        myRegionType;

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
  typedef itk::ImageRegionIteratorWithIndex<myImageType>  myIteratorType;

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
  typedef itk::ImageFileWriter<myImageType>     InputImageWriterType;
  InputImageWriterType::Pointer inputImageWriter= InputImageWriterType::New();
  inputImageWriter->SetFileName( "SyntheticImageForSheetnessTest.mha" );
  inputImageWriter->SetInput( inputImage );
  inputImageWriter->Update();

  // Declare the type for the Hessian filter
  typedef itk::HessianRecursiveGaussianImageFilter<
                                            myImageType >  myHessianFilterType;

  // Declare the type for the sheetness measure filter
  typedef itk::tube::SheetnessMeasureImageFilter< float >  mySheetnessFilterType;

  typedef mySheetnessFilterType::OutputImageType mySheetnessImageType;


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
  typedef mySheetnessFilterType::OutputImageType SheetnessImageType;

  std::cout << "Write out the sheetness image" << std::endl;
  typedef itk::ImageFileWriter<SheetnessImageType>     SheetnessImageWriterType;
  SheetnessImageWriterType::Pointer writer= SheetnessImageWriterType::New();
  writer->SetFileName( "SheetnessImage.mha" );
  writer->SetInput( outputImage );
  writer->Update();

  return EXIT_SUCCESS;
}
