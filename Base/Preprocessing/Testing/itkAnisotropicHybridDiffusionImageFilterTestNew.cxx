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
#include <itkAnisotropicHybridDiffusionImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkMath.h>
#include <itkImageDuplicator.h>
#include <itkMersenneTwisterRandomVariateGenerator.h>

int itkAnisotropicHybridDiffusionImageFilterTestNew(int argc, char* argv []  ) 
{
  // Define the dimension of the images
  const unsigned int Dimension = 3;
  typedef double      InputPixelType;
  typedef double      OutputPixelType;

  // Declare the types of the images
  typedef itk::Image< InputPixelType, Dimension>       InputImageType;
  typedef itk::Image< InputPixelType, Dimension>       OutputImageType;

  typedef itk::ImageFileReader< InputImageType  >      ImageReaderType;

  // Read input image
  ImageReaderType::Pointer   reader = ImageReaderType::New();
  reader->SetFileName ( argv[1] ); 

  std::cout << "Reading input image : " << argv[1] << std::endl;
  try
    {
    reader->Update();
    }
  catch ( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }

  // Declare the anisotropic diffusion edge enhancement filter
  typedef itk::AnisotropicHybridDiffusionImageFilter< InputImageType,
                                            OutputImageType>  HybridFilterType;

  // Create a edge enhancement Filter
  HybridFilterType::Pointer HybridFilter = 
                                      HybridFilterType::New();
  
  //HybridFilter->DebugOn();

  HybridFilter->SetInput( reader->GetOutput() );

  //Set/Get hybrid filter parameters
  
  //Sigma 
  if( argc > 3 )
    {
    double sigma = atof( argv[3] );
    HybridFilter->SetSigma( sigma );
    }
  
  //Contrast EED 
  if( argc > 4 )
    {
    double contrastEED = atof( argv[4] );
    HybridFilter->SetContrastParameterLambdaEED( contrastEED );
    }
 
  //Contrast CED 
  if( argc > 5 )
    {
    double contrastCED = atof( argv[5] );
    HybridFilter->SetContrastParameterLambdaCED( contrastCED );
    }
 
  //Contrast Hybrid 
  if( argc > 6 )
    {
    double contrastHybrid = atof( argv[6] );
    HybridFilter->SetContrastParameterLambdaHybrid( contrastHybrid );
    }
 
  //alpha 
  if( argc > 7 )
    {
    double alpha = atof( argv[7] );
    HybridFilter->SetAlpha( alpha );
    }
  
  // timestep 
  if( argc > 8 )
    {
    double timestep = atof( argv[8] );
    HybridFilter->SetTimeStep( timestep );
    }

  // Number Of iterations
  if( argc > 9 )
    {
    int numberOfIterations= atoi( argv[9] );
    HybridFilter->SetNumberOfIterations( numberOfIterations );
    }

 
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

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}




