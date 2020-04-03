/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#include "itktubeImageRegionMomentsCalculator.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

int itktubeImageRegionMomentsCalculatorTest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0]
      << " inputImage outputImage"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 2 };

  // Define the pixel type
  typedef float PixelType;

  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;


  // Declare the type for the Filter
  typedef itk::tube::ImageRegionMomentsCalculator< ImageType > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }

  ImageType::Pointer inputImage = reader->GetOutput();

  FilterType::Pointer filter = FilterType::New();
  filter->SetImage( inputImage );
  filter->Compute();

  std::cout << "Mass = " << filter->GetTotalMass() << std::endl;
  if( filter->GetTotalMass() != 117300 )
    {
    return EXIT_FAILURE;
    }
  std::cout << "First = " << filter->GetFirstMoments() << std::endl;
  if( std::fabs( filter->GetFirstMoments()[0] - 71.5565 ) > 0.001 )
    {
    return EXIT_FAILURE;
    }
  std::cout << "Second = " << filter->GetSecondMoments() << std::endl;
  std::cout << "CoG = " << filter->GetCenterOfGravity() << std::endl;
  std::cout << "Moments = " << filter->GetCentralMoments() << std::endl;
  std::cout << "Principal Moments = " << filter->GetPrincipalMoments()
            << std::endl;
  std::cout << "Principal Axes = " << filter->GetPrincipalAxes()
            << std::endl;
  std::cout << "Principal Axes Transform = "
            << filter->GetPrincipalAxesToPhysicalAxesTransform()
            << std::endl;
  std::cout << "Physical Axes Transform = "
            << filter->GetPhysicalAxesToPrincipalAxesTransform()
            << std::endl;

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
