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

#include "itktubeCVTImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#define numCentroids 4

int itktubeCVTImageFilterTest( int argc, char * argv[] )
{
  if( argc != 3 )
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
  typedef itk::ImageFileWriter< ImageType > WriterType;


  // Declare the type for the Filter
  typedef itk::tube::CVTImageFilter< ImageType > FilterType;

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
  filter->SetInput( inputImage );
  filter->SetNumberOfCentroids( numCentroids );
  filter->SetInitialSamplingMethod( FilterType::CVT_GRID );
  filter->SetNumberOfSamples( 1000 );
  filter->SetNumberOfIterations( 500 );
  filter->SetNumberOfIterationsPerBatch( 10 );
  filter->SetNumberOfSamplesPerBatch( 100 );
  filter->SetBatchSamplingMethod( FilterType::CVT_RANDOM );
  filter->SetSeed( 1 );
  filter->Update();

  double val[numCentroids];
  for( unsigned int i=0; i<numCentroids; i++ )
    {
    val[i] = 0;
    }

  ImageType::Pointer outImage = filter->GetOutput();

  itk::ImageRegionIterator< ImageType > inItr( inputImage,
    inputImage->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > outItr( outImage,
    outImage->GetLargestPossibleRegion() );
  while( !inItr.IsAtEnd() )
    {
    val[ ( int )( outItr.Get() )-1 ] += inItr.Get();
    ++inItr;
    ++outItr;
    }

  double mean = 0;
  for( unsigned int i=0; i<numCentroids; i++ )
    {
    mean += val[i];
    }
  mean /= numCentroids;
  std::cout << "Mean val = " << mean << std::endl;

  bool valid = true;
  for( unsigned int i=0; i<numCentroids; i++ )
    {
    double tf = vnl_math_abs( ( val[i]-mean )/mean );
    std::cout << "val[" << i << "] = " << val[i] << " ( "
      << ( int )( tf*100 ) << "% diff from mean )" << std::endl;
    if( tf > 0.15 )
      {
      std::cout << "  Error: not within tolerance of mean." << std::endl;
      valid = false;
      }
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( argv[2] );
  writer->SetUseCompression( true );
  writer->SetInput( filter->GetOutput() );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject& e )
    {
    std::cerr << "Exception caught during write:\n"  << e;
    return EXIT_FAILURE;
    }

  // All objects should be automatically destroyed at this point
  if( valid )
    {
    return EXIT_SUCCESS;
    }
  else
    {
    return EXIT_FAILURE;
    }
}
