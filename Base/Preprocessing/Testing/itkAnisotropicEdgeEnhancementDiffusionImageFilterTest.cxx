/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicEdgeEnhancementDiffusionImageFilterTest.cxx,v $
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


#include "itkAnisotropicEdgeEnhancementDiffusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int itkAnisotropicEdgeEnhancementDiffusionImageFilterTest(int argc, char* argv [] )
{
  if ( argc < 3 )
    {
    std::cerr << "Missing Parameters: " 
              << argv[0]
              << " Input_Image"
              << " Edge_Enhanced_Output_Image "<< std::endl; 
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
  typedef itk::AnisotropicEdgeEnhancementDiffusionImageFilter< InputImageType,
                                            OutputImageType>  EdgeEnhancementFilterType;

  // Create a edge enhancement Filter
  EdgeEnhancementFilterType::Pointer EdgeEnhancementFilter = 
                                      EdgeEnhancementFilterType::New();
  
 //EdgeEnhancementFilter->DebugOn();

  EdgeEnhancementFilter->SetInput( reader->GetOutput() );

  //Test Set/Get VED parameters
 
  EdgeEnhancementFilter->SetSensitivity( 4.0 );
  EdgeEnhancementFilter->SetWStrength( 24.0 );
  EdgeEnhancementFilter->SetEpsilon( 0.01 );

  if ( fabs( EdgeEnhancementFilter->GetSensitivity() - 4.0 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  if ( fabs( EdgeEnhancementFilter->GetWStrength() - 24.0 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  if ( fabs( EdgeEnhancementFilter->GetEpsilon() - 0.01 ) > 0.01 )
    {
    std::cerr << "Error Set/Get Sensitivity" << std::endl;
    return EXIT_FAILURE;
    }

  EdgeEnhancementFilter->SetSensitivity( 5.0 );
  EdgeEnhancementFilter->SetWStrength( 25.0 );
  EdgeEnhancementFilter->SetEpsilon( 10e-2 );

  std::cout << "Enhancing vessels.........: " << argv[1] << std::endl;

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

