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
  typedef FunctionType::HistogramType                              HistogramType;
  typedef itk::ImageFileReader< HistogramType >                    HistReaderType;

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

  HistReaderType::Pointer histReader = HistReaderType::New();
  histReader->SetFileName( argv[3] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during addMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer addMeanHist = histReader->GetOutput();

  histReader = HistReaderType::New();
  histReader->SetFileName( argv[4] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during addStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer addStdevHist = histReader->GetOutput();

  histReader = HistReaderType::New();
  histReader->SetFileName( argv[4] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during subMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer subMeanHist = histReader->GetOutput();

  histReader = HistReaderType::New();
  histReader->SetFileName( argv[5] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during subStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer subStdevHist = histReader->GetOutput();

  histReader = HistReaderType::New();
  histReader->SetFileName( argv[6] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during nomMeanHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer nomMeanHist = histReader->GetOutput();

  histReader = HistReaderType::New();
  histReader->SetFileName( argv[7] );
  try
    {
    histReader->Update();
    }
  catch (itk::ExceptionObject& e)
    {
    std::cerr << "Exception caught during nomStdevHist read:\n"  << e;
    return EXIT_FAILURE;
    }
  HistogramType::Pointer nomStdevHist = histReader->GetOutput();

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
  double truth[] = {
    0.459082,0.384968,0.213724,0.0103309,0.256258,0.271751,99.1086,
    24.694,
    0.735631,0.980133,0.0328008,0.00943892,0.255945,0.272139,99.1066,
    24.692,
    0.878761,0.129995,0.000634503,0.00784074,0.256168,0.272067,99.1072,
    24.6925,
    0.496005,-0.107331,-0.177035,0.00708768,0.255856,0.272473,99.1052,
    24.6905,
    };


  // Evaluate
  itk::ImageRegionIteratorWithIndex< ImageType > outIter( inputImage, region );
  outIter.GoToBegin();
  unsigned int truthIndex = 0;
  while( !outIter.IsAtEnd() )
    {
    std::vector<double> tf = func->EvaluateAtIndex( outIter.GetIndex() );
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
