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
#include "FANNClassifyCLP.h"

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
  typedef itk::JointHistogramImageFunction<InputImageType>  HistCalcType;

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

  OutputImageType::Pointer output = OutputImageType::New();
  output->SetRegions( curImage->GetLargestPossibleRegion() );
  output->CopyInformation( curImage );
  output->Allocate();
  output->FillBuffer( 0 );

  timeCollector.Start("calculators");

  // Setup the Calculators
    CalculatorType::Pointer inputCalcSmall = CalculatorType::New();
  inputCalcSmall->SetInputImage( curImage );
  double inputDataMin = inputCalcSmall->GetDataMin();
  double inputDataMax = inputCalcSmall->GetDataMin();
  inputCalcSmall->SetDataMin( inputDataMax );
  inputCalcSmall->SetDataMax( inputDataMin );

  CalculatorType::Pointer inputCalcMedium = CalculatorType::New();
  inputCalcMedium->SetInputImage( curImage );
  inputCalcMedium->SetDataMin( inputDataMax );
  inputCalcMedium->SetDataMax( inputDataMin );

  CalculatorType::Pointer inputCalcLarge = CalculatorType::New();
  inputCalcLarge->SetInputImage( curImage );
  inputCalcLarge->SetDataMin( inputDataMax );
  inputCalcLarge->SetDataMax( inputDataMin );

  CalculatorType::Pointer priorCalcSmall = CalculatorType::New();
  priorCalcSmall->SetInputImage( curPrior );
  double priorDataMin = priorCalcSmall->GetDataMin();
  double priorDataMax = priorCalcSmall->GetDataMin();
  priorCalcSmall->SetDataMin( priorDataMax );
  priorCalcSmall->SetDataMax( priorDataMin );

  CalculatorType::Pointer priorCalcMedium = CalculatorType::New();
  priorCalcMedium->SetInputImage( curPrior );
  priorCalcMedium->SetDataMin( priorDataMax );
  priorCalcMedium->SetDataMax( priorDataMin );

  CalculatorType::Pointer priorCalcLarge = CalculatorType::New();
  priorCalcLarge->SetInputImage( curPrior );
  priorCalcLarge->SetDataMin( priorDataMax );
  priorCalcLarge->SetDataMax( priorDataMin );

  // Setup the Joint-Histogram Calculators
  HistCalcType::Pointer addJHCalc = HistCalcType::New();
  addJHCalc->SetInputImage( curImage );
  addJHCalc->SetInputMask( curPrior );
  // need to read from file
  typedef HistCalcType::HistogramType HistogramType;
  typedef itk::ImageFileReader<HistogramType> HistogramReader;
  HistogramReader::Pointer histReader = HistogramReader::New();
  histReader->SetFileName(addMeans);
  histReader->Update();
  addJHCalc->SetMeanHistogram(histReader->GetOutput());
  histReader->SetFileName(addStdDevs);
  histReader->Update();
  addJHCalc->SetStandardDeviationHistogram(histReader->GetOutput());

  HistCalcType::Pointer subJHCalc = HistCalcType::New();
  subJHCalc->SetInputImage( curImage );
  subJHCalc->SetInputMask( curPrior );
  // need to read from file
  histReader->SetFileName(subMeans);
  histReader->Update();
  subJHCalc->SetMeanHistogram(histReader->GetOutput());
  histReader->SetFileName(subStdDevs);
  histReader->Update();
  subJHCalc->SetStandardDeviationHistogram(histReader->GetOutput());

  HistCalcType::Pointer nomJHCalc = HistCalcType::New();
  nomJHCalc->SetInputImage( curImage );
  nomJHCalc->SetInputMask( curPrior );
  // need to read from file
  histReader->SetFileName(nomMeans);
  histReader->Update();
  nomJHCalc->SetMeanHistogram(histReader->GetOutput());
  histReader->SetFileName(nomStdDevs);
  histReader->Update();
  nomJHCalc->SetStandardDeviationHistogram(histReader->GetOutput());

  timeCollector.Stop("calculators");
  
  // Setup the sigmas based on the wire scale
  InputPixelType sigmaMedium = (scale/2)*0.6667;
  InputPixelType sigmaSmall = 0.6667*sigmaMedium;
  InputPixelType sigmaLarge = 1.3333*sigmaMedium;

  inputCalcSmall->SetScale( sigmaSmall );
  inputCalcMedium->SetScale( sigmaMedium );
  inputCalcLarge->SetScale( sigmaLarge );

  priorCalcSmall->SetScale( sigmaSmall );
  priorCalcMedium->SetScale( sigmaMedium );
  priorCalcLarge->SetScale( sigmaLarge );

  progress = 0.12;
  progressReporter.Report( progress );

  NeighborIterType::SizeType radius  = {{ceil(scale/2),
                                                  ceil(scale/2)}};
  NeighborIterType imageItr(radius,
                            curImage,
                            curImage->GetLargestPossibleRegion() );

  NeighborIterType priorItr(radius,
                            curPrior,
                            curPrior->GetLargestPossibleRegion() );

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

  unsigned int stag = 0;
  for(  centerlineItr.GoToBegin(), imageItr.GoToBegin(), priorItr.GoToBegin();
        !centerlineItr.IsAtEnd();
        ++centerlineItr, ++imageItr, ++priorItr, ++stag,
        progress += (1/numCenterlinePoints)*0.88 )
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
      CalculatorType::IndexType curIndex;
      CalculatorType::ContinuousIndexType curContIndex;

      curIndex = centerlineItr.GetIndex();
      curContIndex = curIndex;

      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );
      // holders for features
      double v1sg, v1mg, v1lg, v2g, v3g, v1se, v1me, v1le, v2e, v3e;
      double roundness = 0;
      double curvature = 0;
      double zAdd, zSub, zNom;
      
      v1sg = priorCalcSmall->Ridgeness( curContIndex, roundness,
        curvature );
      v1mg = priorCalcMedium->Ridgeness( curContIndex, roundness,
        curvature );
      v1lg = priorCalcLarge->Ridgeness( curContIndex, roundness,
        curvature );

      v1se = inputCalcSmall->Ridgeness( curContIndex, roundness,
        curvature );
      v1me = inputCalcMedium->Ridgeness( curContIndex, roundness,
        curvature );
      v1le = inputCalcLarge->Ridgeness( curContIndex, roundness,
        curvature );
  
      double pS = priorCalcSmall->Intensity( curIndex );
      double pM = priorCalcMedium->Intensity( curIndex );
      double pL = priorCalcLarge->Intensity( curIndex );

      double iS = inputCalcSmall->Intensity( curIndex );
      double iM = inputCalcMedium->Intensity( curIndex );
      double iL = inputCalcLarge->Intensity( curIndex );

      if( pL != 0 )
        {
        v2g = ( pS - pL ) / pL;
        }
      else
        {
        v2g = 0;
        }
      if( iL != 0 )
        {
        v2e = ( iS - iL ) / iL;
        }
      else
        {
        v2e = 0;
        }
      
      v3g = pM;
      v3e = iM;
      
      zAdd = addJHCalc->Evaluate( curPoint );
      zSub = subJHCalc->Evaluate( curPoint );
      zNom = nomJHCalc->Evaluate( curPoint );
      
      // begin patch features
      NeighborIterType::NeighborhoodType imageN = imageItr.GetNeighborhood();
      NeighborIterType::NeighborhoodType priorN = priorItr.GetNeighborhood();
      
      assert(imageN.Size() == priorN.Size());
      assert(std::distance(imageN.Begin(), imageN.End()) == imageN.Size());
      
      std::vector<float> normalizedImagePatch(imageN.Size());
      std::vector<float> normalizedPriorPatch(imageN.Size());
      
      // compute the mean intensity of the two patches
      float imageMean = std::accumulate(imageN.Begin(), imageN.End(), 0.0) / 
        imageN.Size();
      float priorMean = std::accumulate(priorN.Begin(), priorN.End(), 0.0) / 
        imageN.Size();
      
      std::transform(imageN.Begin(), imageN.End(),
                     normalizedImagePatch.begin(),
                     std::bind2nd(std::minus<float>(), imageMean));
      
      std::transform(priorN.Begin(), priorN.End(),
                       normalizedPriorPatch.begin(),
                     std::bind2nd(std::minus<float>(), priorMean));                    
      
      float imageStdDev = std::sqrt( 
        std::inner_product( normalizedImagePatch.begin(), 
                            normalizedImagePatch.end(),
                            normalizedImagePatch.begin(),
                            0.0 ) / 
        ( imageN.Size() - 1 ) );
      
      float priorStdDev = std::sqrt(
        std::inner_product( normalizedPriorPatch.begin(), 
                            normalizedPriorPatch.end(),
                            normalizedPriorPatch.begin(),
                            0.0 ) /
        ( priorN.Size() - 1 ) );
      
      float imageNorm = std::sqrt(
        std::inner_product( normalizedImagePatch.begin(), 
                            normalizedImagePatch.end(),
                            normalizedImagePatch.begin(), 0.0 ) );
      
      float priorNorm = std::sqrt(
        std::inner_product( normalizedPriorPatch.begin(), 
                            normalizedPriorPatch.end(),
                            normalizedPriorPatch.begin(), 0.0 ) );
      
      if(imageNorm > 0.0)
        {
        std::transform(normalizedImagePatch.begin(), normalizedImagePatch.end(),
                       normalizedImagePatch.begin(),
                       std::bind2nd(std::divides<float>(), imageNorm));
        }
      
      if(priorNorm > 0.0)
        {
        std::transform(normalizedPriorPatch.begin(), normalizedPriorPatch.end(),
                       normalizedPriorPatch.begin(),
                       std::bind2nd(std::divides<float>(), priorNorm));
        }
      
      float crossCorrelation = std::inner_product( normalizedImagePatch.begin(),
                                                   normalizedImagePatch.end(),
                                                   normalizedPriorPatch.begin(),
                                                   0.0 );
      
      std::vector<float> differencePatch(imageN.Size());
      
      std::transform(imageN.Begin(), imageN.End(),
                     priorN.Begin(),
                     differencePatch.begin(),
                     std::minus<float>());
      
      float norm = std::inner_product( differencePatch.begin(), 
                                       differencePatch.end(),
                                       differencePatch.begin(), 0.0 );
      
      float total = std::accumulate( differencePatch.begin(), 
                                     differencePatch.end(), 0.0);
      
      float maxdiff = *std::max_element( differencePatch.begin(),
                                         differencePatch.end() );
      
      float mindiff = *std::min_element( differencePatch.begin(),
                                         differencePatch.end() );
      
      size_t q1 = static_cast<size_t>(.25*differencePatch.size());
      size_t q2 = static_cast<size_t>(.50*differencePatch.size());
      size_t q3 = static_cast<size_t>(.75*differencePatch.size());
      size_t q95 = static_cast<size_t>(.95*differencePatch.size());
      
      std::sort(differencePatch.begin(), differencePatch.end());
      
      float q1val = differencePatch[q1];
      float q2val = differencePatch[q2];
      float q3val = differencePatch[q3];
      float q95val = differencePatch[q95];

      std::vector<float> features;
      features.push_back(v1sg-v1se);
      features.push_back(v1mg-v1me);
      features.push_back(v1lg-v1le); 
      features.push_back(v2g-v2e);
      features.push_back(v3g-v3e);
      features.push_back(zAdd);
      features.push_back(zSub);
      features.push_back(zNom);
      features.push_back(imageMean-priorMean);
      features.push_back(imageStdDev-priorStdDev);
      features.push_back(crossCorrelation);
      features.push_back(mindiff);
      features.push_back(q1val);
      features.push_back(q2val);
      features.push_back(q3val);
      features.push_back(q95val);
      features.push_back(maxdiff);
      features.push_back(total);
      features.push_back(norm);
      
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
