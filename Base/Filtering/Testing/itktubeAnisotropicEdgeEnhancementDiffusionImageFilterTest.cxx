/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itktubeAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#include "itktubeAnisotropicEdgeEnhancementDiffusionImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

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
  typedef double      OutputPixelType;

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
  typedef itk::tube::AnisotropicEdgeEnhancementDiffusionImageFilter< InputImageType,
                                            OutputImageType>  EdgeEnhancementFilterType;

  // Create a edge enhancement Filter
  EdgeEnhancementFilterType::Pointer EdgeEnhancementFilter =
                                      EdgeEnhancementFilterType::New();
  EdgeEnhancementFilter->SetInput( reader->GetOutput() );

  //Set/Get VED parameters

  //Set scale/sigma value
  if( argc > 3 )
  {
  double scaleParameter = std::atof(argv[3]);
  EdgeEnhancementFilter->SetSigma( scaleParameter );
  }
  //Set contrast parameter
  if( argc > 4 )
  {
  double contrastParameter = std::atof(argv[4]);
  EdgeEnhancementFilter->SetContrastParameterLambdaE( contrastParameter );
  }
  //Set time step
  if( argc > 5 )
  {
  double timeStep = std::atof(argv[5]);
  EdgeEnhancementFilter->SetTimeStep( timeStep );
  }
  //Set number of iterations
  if( argc > 6 )
  {
  double numberOfIterations = std::atoi(argv[6]);
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
