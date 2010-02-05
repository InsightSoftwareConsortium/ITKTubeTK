/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkTopHatImageFilterTest.cxx,v $
  Language:  C++
  Date:      $Date: 2005-09-01 13:41:14 $
  Version:   $Revision: 1.1 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif


#include <itkImage.h>
#include "itkFilterWatcher.h"
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkVesselEnhancingDiffusion2DImageFilter.h>

int itkVesselEnhancingDiffusion2DImageFilterTest(int argc, char* argv [] ) 
{
  if( argc != 3 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << "  inputImage outputImage" << std::endl;
    return EXIT_FAILURE;
    }
  
  // Define the dimension of the images
  const unsigned int Dimension = 2;

  // Define the pixel type
  typedef short PixelType;
  
  // Declare the types of the images
  typedef itk::Image<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  
 
  // Declare the type for the Filter
  typedef itk::VesselEnhancingDiffusion2DImageFilter<
                           PixelType, Dimension > FilterType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  
  FilterType::Pointer filter = FilterType::New();
  FilterWatcher watcher(filter, "filter");

  // Connect the pipeline
  filter->SetInput( reader->GetOutput() );
  filter->SetDefaultPars( ); // duplicates assignments given below
  filter->SetTimeStep( 0.25 );
  filter->SetIterations( 20 ); // Default is 200
  filter->SetRecalculateVesselness( 10 ); // Default is 100
  filter->SetBeta( 0.5 );
  filter->SetGamma( 5.0 );
  filter->SetEpsilon( 0.01 );
  filter->SetOmega( 25.0 );
  filter->SetSensitivity( 20.0 );
  std::vector< float > scales;
  scales.resize(2);
  scales[0] = 1;
  scales[1] = 3;
  filter->SetScales( scales );
  filter->SetDarkObjectLightBackground( true );
  filter->SetVerbose( true );
  filter->Update();
  writer->SetInput( filter->GetOutput() );
  
  // Execute the filter
  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during pipeline Update\n"  << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point

  return EXIT_SUCCESS;

}




