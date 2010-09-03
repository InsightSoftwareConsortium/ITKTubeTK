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

#include <vector>

#include <itkImage.h>
#include <itkFilterWatcher.h>
#include <itkExceptionObject.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itkStandardFeatureGeneratingImageFunction.h>

int itkStandardFeatureGeneratingImageFunctionTest(int argc, char* argv [] ) 
{
  if( argc < 8 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " inputImage maskImage addMeanHist addStdevHist "
      "subMeanHist subStdevHist nomMeanHist nomStdevHist" << std::endl;
    return EXIT_FAILURE;
    }
  
  // Define the dimension of the images
  const unsigned int Dimension = 2;

  // Define the pixel type
  typedef float PixelType;
  
  // Declare the types of the images
  typedef itk::OrientedImage<PixelType, Dimension>  ImageType;

  // Declare the reader and writer
  typedef itk::ImageFileReader< ImageType > ReaderType;
  
  // Declare the type for the Filter
  typedef itk::StandardFeatureGeneratingImageFunction< ImageType > FunctionType;

  // Create the reader and writer
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during input read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer inputImage = reader->GetOutput();
  
  reader = ReaderType::New();
  reader->SetFileName( argv[2] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during prior read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer priorImage = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[3] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during addMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer addMeanHist = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[4] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during addStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer addStdevHist = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[4] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during subMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer subMeanHist = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[5] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during subStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer subStdevHist = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[6] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during nomMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer nomMeanHist = reader->GetOutput();

  reader = ReaderType::New();
  reader->SetFileName( argv[7] );
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during nomStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  ImageType::Pointer nomStdevHist = reader->GetOutput();

  FunctionType::Pointer func = FunctionType::New();
  func->SetInputImage( inputImage );
  func->SetPriorImage( priorImage );
  func->SetMeanAddHistogram( addMeanHist );
  func->SetStdevAddHistogram( addStdevHist );
  func->SetMeanSubHistogram( subMeanHist );
  func->SetStdevSubHistogram( subStdevHist );
  func->SetMeanNormHistogram( nomMeanHist ); 
  func->SetStdevNormHistogram( nomStdevHist );
  func->SetScale( 20 );
  func->PrepFilter();

  // shrink the region
  ImageType::RegionType region = inputImage->GetLargestPossibleRegion();
  ImageType::SizeType size;
  size[0] = 2;
  size[1] = 2;
  region.SetSize( size );

  // Create truth set
  double ep = 0.0001;
  double truth[] = {-0.204853,-0.31364,-0.209311,-6.69886e-05,-0.997691,
                    18.9597,1.08886,6.41422,-0.00671102,-0.00494839,
                    -0.00119436,-7.11408e-05,-0.997691,20.3184,1.20397,
                    7.24889,-0.00373373,-9.6231e-05,-2.376e-09,-6.74754e-05,
                    -0.997691,20.1709,1.17282,7.0434,0.0737599,0.107649,
                    0.177815,-7.15384e-05,-0.997691,21.6709,1.29748,7.96332};

  // Evaluate
  itk::ImageRegionIteratorWithIndex< ImageType > outIter( inputImage, region );
  ImageType::PointType pnt;
  outIter.GoToBegin();
  unsigned int truthIndex = 0;
  while( !outIter.IsAtEnd() )
    {
    inputImage->TransformIndexToPhysicalPoint( outIter.GetIndex(), pnt);
    std::vector<double> tf = func->Evaluate( pnt );
    std::vector<double>::const_iterator itr;
    std::cout << "tf =";
    for( itr = tf.begin(); itr != tf.end(); ++itr )
      {
      double cur = *itr;
      double curTruth = truth[truthIndex];
      if( !(cur > curTruth-ep && cur < curTruth+ep) )
        {
        std::cout << "\n" << std::endl;
        std::cout << cur << " != " << curTruth << std::endl;
        return EXIT_FAILURE;
        }
      std::cout << " " << *itr; 
      ++truthIndex;
      }
    std::cout << std::endl;
    ++outIter;
    }

  // All objects should be automatically destroyed at this point
  return EXIT_SUCCESS;
}
