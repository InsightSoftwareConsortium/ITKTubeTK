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
#include "itktubeAnisotropicCoherenceEnhancingDiffusionImageFilter.h"

#include <itkImageFileReader.h>

int itktubeAnisotropicCoherenceEnhancingDiffusionImageFilterTest( int argc,
  char * argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image"
              << " Coherence_Enhanced_Output_Image [Sigma]"
              << " [Alpha] [SigmaOuter] [ContrastParameter]"
              << "[TimeStep] [NumberOfIterations]" << std::endl;
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
  typedef itk::tube::AnisotropicCoherenceEnhancingDiffusionImageFilter<
    InputImageType, OutputImageType>  CoherenceEnhancingFilterType;

  // Create a edge enhancement Filter
  CoherenceEnhancingFilterType::Pointer CoherenceEnhancingFilter =
                                      CoherenceEnhancingFilterType::New();
  //CoherenceEnhancingFilter->DebugOn();

  CoherenceEnhancingFilter->SetInput( reader->GetOutput() );

  //Set/Get Coherence Enhancing Diffusion parameters
  //
  //Set sigma
  if( argc > 3 )
    {
    double sigma = std::atof( argv[3] );
    std::cout << "Set sigma parameter value \t" << sigma << std::endl;
    CoherenceEnhancingFilter->SetSigma( sigma );
    }

  //set sigma outer
  if( argc > 4 )
    {
    double sigmaOuter = std::atof( argv[4] );
    std::cout << "Set outer sigma parameter value \t" << sigmaOuter << std::endl;
    CoherenceEnhancingFilter->SetSigmaOuter( sigmaOuter );
    }

  //set alpha
  if( argc > 5 )
    {
    double alpha = std::atof( argv[5] );
    std::cout << "Set alpha parameter value \t" << alpha << std::endl;
    CoherenceEnhancingFilter->SetAlpha( alpha );
    }

  //Set contrast parameter
  if( argc > 6 )
    {
    double contrastParamter = std::atof( argv[6] );
    std::cout << "Set contrast parameter value \t" << contrastParamter << std::endl;
    CoherenceEnhancingFilter->SetContrastParameterLambdaC( contrastParamter );
    }

  //Set time step
  if( argc > 7 )
    {
    double timeStep = std::atof( argv[7] );
    CoherenceEnhancingFilter->SetTimeStep( timeStep );
    }

  //Set number of iterations
  if( argc > 8 )
    {
    double numberOfIterations = std::atoi( argv[8] );
    CoherenceEnhancingFilter->SetNumberOfIterations( numberOfIterations );
    }

  CoherenceEnhancingFilter->Print ( std::cout );
  std::cout << "Enhancing .........: " << argv[1] << std::endl;

  try
    {
    CoherenceEnhancingFilter->Update();
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
  writer->SetInput ( CoherenceEnhancingFilter->GetOutput() );

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
