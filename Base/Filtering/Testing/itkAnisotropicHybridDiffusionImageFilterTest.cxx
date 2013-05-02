/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkAnisotropicHybridDiffusionImageFilterTest.cxx,v $
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


#include "itkAnisotropicHybridDiffusionImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

int itkAnisotropicHybridDiffusionImageFilterTest(int argc, char* argv[] )
{
  if( argc < 3 )
    {
    std::cerr << "Missing Parameters: "
              << argv[0]
              << " Input_Image"
              << " HDCS_Enhanced_Output_Image "
              << " [Sigma] [Sigma Outer] [EED contrast] [CED contrast] [Hybrid contrast] "
              << " [Alpha] [Time step size] [ Number of Iterations]"
              << std::endl;
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
  catch( itk::ExceptionObject &err )
    {
    std::cerr << "Exception thrown: " << err << std::endl;
    return EXIT_FAILURE;
    }


  // Declare the anisotropic diffusion hybrid enhancement filter
  typedef itk::AnisotropicHybridDiffusionImageFilter< InputImageType,
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
    double sigma = atof( argv[3] );
    HybridFilter->SetSigma( sigma );
    }
  if( argc > 4 )
    {
    double sigmaOuter = atof( argv[4] );
    HybridFilter->SetSigmaOuter( sigmaOuter );
    }
  //Contrast EED
  if( argc > 5 )
    {
    double contrastEED = atof( argv[5] );
    HybridFilter->SetContrastParameterLambdaEED( contrastEED );
    }
  //Contrast CED
  if( argc > 6 )
    {
    double contrastCED = atof( argv[6] );
    HybridFilter->SetContrastParameterLambdaCED( contrastCED );
    }
  //Contrast Hybrid
  if( argc > 7 )
    {
    double contrastHybrid = atof( argv[7] );
    HybridFilter->SetContrastParameterLambdaHybrid( contrastHybrid );
    }
  //alpha
  if( argc > 8 )
    {
    double alpha = atof( argv[8] );
    HybridFilter->SetAlpha( alpha );
    }
  // time step size
  if( argc > 9 )
    {
    double timestep = atof( argv[9] );
    HybridFilter->SetTimeStep( timestep );
    }

  // Number Of iterations
  if( argc > 10 )
    {
    int numberOfIterations= atoi( argv[10] );
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
