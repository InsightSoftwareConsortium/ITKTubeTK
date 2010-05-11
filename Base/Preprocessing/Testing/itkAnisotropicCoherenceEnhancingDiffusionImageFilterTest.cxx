/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2007/06/20 16:03:23 $
  Version:   $Revision: 1.5 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include "itkAnisotropicCoherenceEnhancingDiffusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int itkAnisotropicCoherenceEnhancingDiffusionImageFilterTest(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Image"
              << " Edge_Enhanced_Output_Image [Alpha] [ContrastParameter]"<< std::endl; 
    return EXIT_FAILURE;
    }
 
 
  // Define the dimension of the images
  const unsigned int Dimension = 3;
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
  catch ( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // Declare the anisotropic diffusion edge enhancement filter
  typedef itk::AnisotropicCoherenceEnhancingDiffusionImageFilter< InputImageType,
                                            OutputImageType>  CoherenceEnhancingFilterType;

  // Create a edge enhancement Filter
  CoherenceEnhancingFilterType::Pointer CoherenceEnhancingFilter = 
                                      CoherenceEnhancingFilterType::New();
  
  //CoherenceEnhancingFilter->DebugOn();

  CoherenceEnhancingFilter->SetInput( reader->GetOutput() );

  //Set/Get Coherence Enhancing Diffusion parameters
  //Set alpha
  if( argc > 3 ) 
    {
    double alpha = atof( argv[3] );
    std::cout << "Set alpha parameter value \t" << alpha << std::endl;
    CoherenceEnhancingFilter->SetAlpha( alpha );
    }

  //Set contrast parameter
  if( argc > 4 ) 
    {
    double contrastParamter = atof( argv[4] );
    std::cout << "Set contrast parameter value \t" << contrastParamter << std::endl;
    CoherenceEnhancingFilter->SetContrastParameterLambdaC( contrastParamter );
    }

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

