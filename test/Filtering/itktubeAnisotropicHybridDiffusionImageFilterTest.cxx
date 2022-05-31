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
#include "itktubeAnisotropicHybridDiffusionImageFilter.h"

#include <itkImageFileReader.h>

int itktubeAnisotropicHybridDiffusionImageFilterTest( int argc,
  char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
      << argv[0]
      << " Input_Image"
      << " HDCS_Enhanced_Output_Image "
      << " [Sigma] [Sigma Outer] [EED contrast] [CED contrast]"
      << " [Hybrid contrast] "
      << " [Alpha] [Time step size] [ Number of Iterations]"
      << std::endl;
    return EXIT_FAILURE;
    }

  // Define the dimension of the images
  enum { Dimension = 3 };
  typedef double      InputPixelType;

  // Declare the types of the images
  typedef itk::Image< InputPixelType, Dimension>           InputImageType;
  typedef itk::Image< InputPixelType, Dimension>           OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >      ImageReaderType;

  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] );

  std::cout << "Reading input image : " << argv[1] << std::endl;
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // Declare the anisotropic diffusion hybrid enhancement filter
  typedef itk::tube::AnisotropicHybridDiffusionImageFilter< InputImageType,
                                            OutputImageType>  HybridFilterType;

  // Create a hybrid diffusion Filter
  HybridFilterType::Pointer HybridFilter =
                                      HybridFilterType::New();
  //HybridFilter->DebugOn();

  HybridFilter->SetInput( reader->GetOutput() );

  //Set/Get hybrid filter parameters
  //Sigma
  if( argc > 3 )
    {
    double sigma = std::atof( argv[3] );
    HybridFilter->SetSigma( sigma );
    }
  if( argc > 4 )
    {
    double sigmaOuter = std::atof( argv[4] );
    HybridFilter->SetSigmaOuter( sigmaOuter );
    }
  //Contrast EED
  if( argc > 5 )
    {
    double contrastEED = std::atof( argv[5] );
    HybridFilter->SetContrastParameterLambdaEED( contrastEED );
    }
  //Contrast CED
  if( argc > 6 )
    {
    double contrastCED = std::atof( argv[6] );
    HybridFilter->SetContrastParameterLambdaCED( contrastCED );
    }
  //Contrast Hybrid
  if( argc > 7 )
    {
    double contrastHybrid = std::atof( argv[7] );
    HybridFilter->SetContrastParameterLambdaHybrid( contrastHybrid );
    }
  //alpha
  if( argc > 8 )
    {
    double alpha = std::atof( argv[8] );
    HybridFilter->SetAlpha( alpha );
    }
  // time step size
  if( argc > 9 )
    {
    double timestep = std::atof( argv[9] );
    HybridFilter->SetTimeStep( timestep );
    }

  // Number Of iterations
  if( argc > 10 )
    {
    int numberOfIterations= std::atoi( argv[10] );
    HybridFilter->SetNumberOfIterations( numberOfIterations );
    }

  HybridFilter->Print( std::cout );
  std::cout << "Enhancing .........: " << argv[1] << std::endl;

  try
    {
    HybridFilter->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << "Writing out the enhanced image to " <<  argv[2] << std::endl;

  typedef itk::ImageFileWriter< OutputImageType  >      ImageWriterType;
  ImageWriterType::Pointer writer = ImageWriterType::New();

  writer->SetFileName( argv[2] );
  writer->SetInput ( HybridFilter->GetOutput() );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
