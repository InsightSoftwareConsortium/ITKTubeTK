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

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// It is important to use OrientedImages
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

#include "itkMersenneTwisterRandomVariateGenerator.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Project includes
#include "itkRidgeExtractor.h"
#include "itkJointHistogramImageFunction.h"
#include "itkStandardFeatureGeneratingImageFunction.h"
#include "itkPatchFeatureGeneratingImageFunction.h"

// ITK filters used
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkDivideImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkNeighborhoodIterator.h"

#include "floatfann.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>
#include <cassert>

// Must include CLP before including tubeCLIHleperFunctions
#include "tubeFANNClassifyCLP.h"

// No need for the DoIt nonsense (we arent' templated over input type)
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter( "FANNClassify",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef unsigned char                                 OutputPixelType;
  typedef float                                         InputPixelType;
  const unsigned int                                    Dimension = 2;
  typedef itk::Image< InputPixelType,  Dimension >      InputImageType;
  typedef itk::Image< OutputPixelType,  Dimension >     OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType >       WriterType;
  typedef itk::ImageFileReader< InputImageType >        ReaderType;
  typedef itk::ImageFileReader< OutputImageType >       MaskReaderType;

  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start("Load data");
  fann* network = fann_create_from_file(inputFile.c_str());

  ReaderType::Pointer reader = ReaderType::New();
  ReaderType::Pointer priorReader = ReaderType::New();
  MaskReaderType::Pointer centerlinesReader = MaskReaderType::New();

  reader->SetFileName(inputVolume);
  reader->Update();
  priorReader->SetFileName(inputPrior);
  priorReader->Update();
  centerlinesReader->SetFileName(inputCenterlines);
  centerlinesReader->Update();

  timeCollector.Stop("Load data");

  // Get the centerlines
  OutputImageType::Pointer centerlines = centerlinesReader->GetOutput();

  // typedefs for numerics
  typedef itk::RidgeExtractor< InputImageType >             CalculatorType;

  // typedefs for classification
  typedef itk::StandardFeatureGeneratingImageFunction< InputImageType >
                                                           StandardFunctionType;
  typedef itk::PatchFeatureGeneratingImageFunction< InputImageType >
                                                           PatchFunctionType;
  typedef StandardFunctionType::HistogramType              HistogramType;
  typedef itk::ImageFileReader<HistogramType>              HistReaderType;

  // typedefs for filters
  typedef itk::RescaleIntensityImageFilter<InputImageType,InputImageType>
                                                            RescaleType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<OutputImageType>
                                                            IterType;
  typedef itk::NeighborhoodIterator<InputImageType>         NeighborIterType;

  // Rescale the input
  InputImageType::Pointer curImage = NULL;
  if(!skipNormalizeImage)
    {
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput( reader->GetOutput() );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();
    curImage = rescale->GetOutput();
    }
  else
    {
    curImage = reader->GetOutput();
    }

  // Rescale the prior
  InputImageType::Pointer curPrior = NULL;
  if(!skipNormalizeImage)
    {
    RescaleType::Pointer rescale = RescaleType::New();
    rescale->SetInput( priorReader->GetOutput() );
    rescale->SetOutputMinimum( 0 );
    rescale->SetOutputMaximum( 1 );
    rescale->Update();
    curPrior = rescale->GetOutput();
    }
  else
    {
    curPrior = priorReader->GetOutput();
    }

  // Read the histograms
  HistReaderType::Pointer histReader = HistReaderType::New();
  histReader->SetFileName( addMeans );
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
  histReader->SetFileName( addStdDevs );
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
  histReader->SetFileName( subMeans );
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
  histReader->SetFileName( subStdDevs);
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
  histReader->SetFileName( nomMeans );
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
  histReader->SetFileName( nomStdDevs );
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

  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( curImage->GetLargestPossibleRegion() );
  output->CopyInformation( curImage );
  output->Allocate();
  output->FillBuffer( 0 );

  progress = 0.12;
  progressReporter.Report( progress );

  IterType centerlineItr( centerlines,
                          centerlines->GetLargestPossibleRegion() );

  double numCenterlinePoints = 0;
  for( centerlineItr.GoToBegin(); !centerlineItr.IsAtEnd(); ++centerlineItr )
    {
    ++numCenterlinePoints;
    }

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer
    randGen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  if( seed != -1 )
    {
    randGen->Initialize( seed );
    }
  else
    {
    randGen->Initialize();
    }

  StandardFunctionType::Pointer standGen = StandardFunctionType::New();
  PatchFunctionType::Pointer patchGen = PatchFunctionType::New();

  standGen->SetInputImage( curImage );
  standGen->SetPriorImage( curPrior );
  standGen->SetMeanAddHistogram( addMeanHist );
  standGen->SetStdevAddHistogram( addStdevHist );
  standGen->SetMeanSubHistogram( subMeanHist );
  standGen->SetStdevSubHistogram( subStdevHist );
  standGen->SetMeanNormHistogram( nomMeanHist );
  standGen->SetStdevNormHistogram( nomStdevHist );
  standGen->SetScale( scale );
  standGen->PrepFilter();

  patchGen->SetInputImage( curImage );
  patchGen->SetPriorImage( curPrior );
  patchGen->SetWidth( 5 );

  unsigned int stag = 0;
  for( centerlineItr.GoToBegin(); !centerlineItr.IsAtEnd();
       ++centerlineItr, ++stag, progress += (1/numCenterlinePoints)*0.88 )
    {

    // Staggered Progress Reporting
    if( stag == 10000 )
      {
      progressReporter.Report( progress );
      stag = 0;
      }

    if( centerlineItr.Get() > 0 &&
        static_cast<int>(randGen->GetUniformVariate( 0, subsampling )) == 0 )
      {
      InputImageType::PointType curPoint;
      InputImageType::IndexType curIndex;

      curIndex = centerlineItr.GetIndex();

      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );

      std::vector<float> features;
      std::vector<double>::const_iterator itr;

      std::vector<double> standFeats = standGen->Evaluate( curPoint );
      std::vector<double> patchFeats = patchGen->Evaluate( curPoint );

      for( itr = standFeats.begin(); itr != standFeats.end(); ++itr )
        {
        features.push_back(*itr);
        }
      for( itr = patchFeats.begin(); itr != patchFeats.end(); ++itr )
        {
        features.push_back(*itr);
        }

      float* rawData = &(features[0]);
      assert( features.size() == fann_get_num_input( network ) );
      float* pixeloutput = fann_run( network, rawData );
      unsigned int nOutputs = fann_get_num_output( network );
      assert( nOutputs == 3 );

      output->TransformPhysicalPointToIndex( curPoint, curIndex);

      // bad hard code
      unsigned char label = 1 + std::distance( pixeloutput,
                                               std::max_element( pixeloutput,
                                                                 pixeloutput +
                                                                 nOutputs ) );

      output->SetPixel(curIndex, label);

      }
    }

  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(output);
  writer->SetFileName(outputVolume);
  writer->Update();

  fann_destroy(network);

  progressReporter.End();
  timeCollector.Report();

  return EXIT_SUCCESS;
}
