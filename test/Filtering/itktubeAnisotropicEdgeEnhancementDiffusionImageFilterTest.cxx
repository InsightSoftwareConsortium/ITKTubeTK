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
#include "itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h"

#include <itkImageFileReader.h>

int itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest( int argc, char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image"
              << " Edge_Enhanced_Output_Image [ScaleParameter] [ContrastParameter] "
              << " [TimeStep] [NumberOfIterations]" << std::endl;
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


  // Declare the anisotropic diffusion edge enhancement filter
  typedef itk::tube::AnisotropicEdgeEnhancementDiffusionImageFilter<
    InputImageType, OutputImageType>  EdgeEnhancementFilterType;

  // Create a edge enhancement Filter
  EdgeEnhancementFilterType::Pointer EdgeEnhancementFilter =
    EdgeEnhancementFilterType::New();
  EdgeEnhancementFilter->SetInput( reader->GetOutput() );

  //Set/Get VED parameters

  //Set scale/sigma value
  if( argc > 3 )
    {
    double scaleParameter = std::atof( argv[3] );
    EdgeEnhancementFilter->SetSigma( scaleParameter );
    }
  //Set contrast parameter
  if( argc > 4 )
    {
    double contrastParameter = std::atof( argv[4] );
    EdgeEnhancementFilter->SetContrastParameterLambdaE( contrastParameter );
    }
  //Set time step
  if( argc > 5 )
    {
    double timeStep = std::atof( argv[5] );
    EdgeEnhancementFilter->SetTimeStep( timeStep );
    }
  //Set number of iterations
  if( argc > 6 )
    {
    double numberOfIterations = std::atoi( argv[6] );
    EdgeEnhancementFilter->SetNumberOfIterations( numberOfIterations );
    }

  std::cout << "Enhancing .........: " << argv[1] << std::endl;

  EdgeEnhancementFilter->Print( std::cout );

  try
    {
    EdgeEnhancementFilter->Update();
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
  writer->SetInput ( EdgeEnhancementFilter->GetOutput() );

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
