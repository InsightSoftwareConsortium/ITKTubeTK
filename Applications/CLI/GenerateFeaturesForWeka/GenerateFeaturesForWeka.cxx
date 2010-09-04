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
#include "itkOrientedImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following three should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
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
#include "vnl/vnl_math.h"
#include "math.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "GenerateFeaturesForWekaCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;
  
  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "GenerateFeaturesForWeka",
                                              CLPProcessInformation );
  progressReporter.Start();

  // typedefs for data structures
  typedef float                                            PixelType;
  typedef itk::OrientedImage<PixelType, dimensionT>        ImageType;
  typedef itk::ImageFileReader<ImageType>                  ReaderType;
  typedef itk::ImageFileWriter<ImageType>                  WriterType;

  // typedefs for numerics
  typedef itk::RidgeExtractor< ImageType >                 CalculatorType;
  typedef itk::JointHistogramImageFunction< ImageType >    HistCalcType;

  // typedefs for classification
  typedef itk::StandardFeatureGeneratingImageFunction< ImageType > 
                                                           StandardFunctionType;
  typedef itk::PatchFeatureGeneratingImageFunction< ImageType > 
                                                           PatchFunctionType;

  // typedefs for filters
  typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>     IterType;

  // Setup the readers to load the input data (image + prior)
  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  typename ReaderType::Pointer priorReader = ReaderType::New();
  typename ReaderType::Pointer additionsReader = ReaderType::New();
  typename ReaderType::Pointer subtractionsReader = ReaderType::New();
  typename ReaderType::Pointer centerlinesReader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  priorReader->SetFileName( inputPrior.c_str() );
  additionsReader->SetFileName( inputAdditions.c_str() );
  subtractionsReader->SetFileName( inputSubtractions.c_str() );
  centerlinesReader->SetFileName( inputCenterlines.c_str() );

  // Load the input image with exception handling
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input image!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    priorReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input prior!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    additionsReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input defects!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    subtractionsReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input defects!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  // Load the input prior with exception handling
  try
    {
    centerlinesReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught while reading the input centerlines!"
              << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  // Get the centerlines
  typename ImageType::Pointer centerlines = centerlinesReader->GetOutput();
  typename ImageType::Pointer additions = additionsReader->GetOutput();
  typename ImageType::Pointer subtractions = subtractionsReader->GetOutput();

  // Rescale the input
  typename ImageType::Pointer curImage = NULL;
  if(!skipNormalizeImage)
    {
    typename RescaleType::Pointer rescale = RescaleType::New();
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
  typename ImageType::Pointer curPrior = NULL;
  if(!skipNormalizeImage)
    {
    typename RescaleType::Pointer rescale = RescaleType::New();
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


  double samplesNom = 0;
  double samplesAdd = 0;
  double samplesSub = 0;

  // Setup the sigmas based on the wire scale
  PixelType sigmaMedium = (scale/2)*0.6667;
  PixelType sigmaSmall = 0.6667*sigmaMedium;

  // Init joint histogram calculators for training
  typename HistCalcType::Pointer addJHCalc = HistCalcType::New();
  typename HistCalcType::Pointer subJHCalc = HistCalcType::New();
  typename HistCalcType::Pointer nomJHCalc = HistCalcType::New();
  addJHCalc->SetInputImage( curImage );
  addJHCalc->SetInputMask( curPrior );
  subJHCalc->SetInputImage( curImage );
  subJHCalc->SetInputMask( curPrior );
  nomJHCalc->SetInputImage( curImage );
  nomJHCalc->SetInputMask( curPrior );

  // Init add and sub calc for class checking
  typename CalculatorType::Pointer addCalc = CalculatorType::New();
  typename CalculatorType::Pointer subCalc = CalculatorType::New();
  addCalc->SetInputImage( additions );
  subCalc->SetInputImage( subtractions );
  addCalc->SetScale( sigmaSmall );
  subCalc->SetScale( sigmaSmall );

  IterType centerlineItr( centerlines, 
                          centerlines->GetLargestPossibleRegion() );

  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename ImageType::PointType curPoint;
      typename CalculatorType::IndexType curIndex;
      typename CalculatorType::ContinuousIndexType curContIndex;
  
      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );
  
      if( curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) 
            && addCalc->Intensity( curIndex ) )
          {
          ++samplesAdd;
          }
        else if( subtractions->TransformPhysicalPointToIndex( curPoint,
                                                              curIndex ) 
                 && subCalc->Intensity( curIndex ) )
          {
          ++samplesSub;
          }
        else
          {
          ++samplesNom;
          }
        }
      }
    ++centerlineItr;
    }

  double sampleMin = samplesAdd;
  if( samplesSub < sampleMin )
    {
    sampleMin = samplesSub;
    }
  if( samplesNom < sampleMin )
    {
    sampleMin = samplesNom;
    }
  double sampleRateAdd = sampleMin/samplesAdd;
  double sampleRateSub = sampleMin/samplesSub;
  double sampleRateNom = sampleMin/samplesNom;

  // train joint histograms
  samplesNom = 0;
  samplesAdd = 0;
  samplesSub = 0;
  double samplesTotal = 0;
  centerlineItr.GoToBegin();
  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename ImageType::PointType curPoint;
      typename CalculatorType::IndexType curIndex;

      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );
      
      if( curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) 
            && addCalc->Intensity( curIndex ) )
          {
          if( (int)(samplesAdd) != (int)(samplesAdd+sampleRateAdd) )
            {
            addJHCalc->Precompute( curPoint );
            }
          samplesAdd += sampleRateAdd;
          }
        else if( subtractions->TransformPhysicalPointToIndex( curPoint, 
                                                            curIndex ) 
                 && subCalc->Intensity( curIndex ) )
          {
          if( (int)(samplesSub) != (int)(samplesSub+sampleRateSub) )
            {
            subJHCalc->Precompute( curPoint );
            }
          samplesSub += sampleRateSub;
          }
        else 
          {
          if( (int)(samplesNom) != (int)(samplesNom+sampleRateNom) )
            {
            nomJHCalc->Precompute( curPoint );
            }
          samplesNom += sampleRateNom;
          }
        }
      }
    ++samplesTotal;
    ++centerlineItr;
    }

  // write out histogram files
  typedef typename HistCalcType::HistogramType          HistrogramType;
  typedef typename itk::ImageFileWriter<HistrogramType> HistogramWriter;
  typename HistogramWriter::Pointer histWriter = HistogramWriter::New();

  addJHCalc->ComputeMeanAndStandardDeviation();
  histWriter->SetInput(addJHCalc->GetMeanHistogram());
  histWriter->SetFileName(addMeans);
  histWriter->Update();

  histWriter->SetInput(addJHCalc->GetStandardDeviationHistogram());
  histWriter->SetFileName(addStdDevs);
  histWriter->Update();

  subJHCalc->ComputeMeanAndStandardDeviation();
  histWriter->SetInput(subJHCalc->GetMeanHistogram());
  histWriter->SetFileName(subMeans);
  histWriter->Update();

  histWriter->SetInput(subJHCalc->GetStandardDeviationHistogram());
  histWriter->SetFileName(subStdDevs);
  histWriter->Update();

  nomJHCalc->ComputeMeanAndStandardDeviation();
  histWriter->SetInput(nomJHCalc->GetMeanHistogram());
  histWriter->SetFileName(nomMeans);
  histWriter->Update();

  histWriter->SetInput(nomJHCalc->GetStandardDeviationHistogram());
  histWriter->SetFileName(nomStdDevs);
  histWriter->Update();
  
  // Get the outfile started
  std::ofstream output( outputFile.c_str() );
  output << "% Generated by GenerateFeaturesForWeka in TubeTK\n";
  output << "@RELATION errors\n";
  output << "\n";
  output << "@ATTRIBUTE x NUMERIC\n";
  output << "@ATTRIBUTE y NUMERIC\n";

  output << "@ATTRIBUTE RidgenessDiffSmall NUMERIC\n";
  output << "@ATTRIBUTE RidgenessDiffMedium NUMERIC\n";
  output << "@ATTRIBUTE RidgenessDiffLarge NUMERIC\n";

  output << "@ATTRIBUTE NormalizedScaledIntensityDiff NUMERIC\n";
  output << "@ATTRIBUTE ScaledIntensityDiff NUMERIC\n";

  output << "@ATTRIBUTE Z_add NUMERIC\n";
  output << "@ATTRIBUTE Z_sub NUMERIC\n";
  output << "@ATTRIBUTE Z_nom NUMERIC\n";

  output << "@ATTRIBUTE PatchMean-PriorMean NUMERIC\n";
  output << "@ATTRIBUTE PatchStdDev-PriorStdDev NUMERIC\n";
  output << "@ATTRIBUTE PatchCrossCorrelation NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceMin NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceP25 NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceP50 NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceP75 NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceP95 NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceMax NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceSum NUMERIC\n";
  output << "@ATTRIBUTE PatchDifferenceL2Norm NUMERIC\n";

  output << "@ATTRIBUTE class {typical,addition,subtraction}\n";
  output << "\n";
  output << "@DATA\n";


  samplesNom = 0;
  samplesAdd = 0;
  samplesSub = 0;
  double portion = 0.9;
  double step = portion/samplesTotal;
  unsigned int count = 0;

  typename StandardFunctionType::Pointer standGen = StandardFunctionType::New();
  typename PatchFunctionType::Pointer patchGen = PatchFunctionType::New();
  
  standGen->SetInputImage( curImage );
  standGen->SetPriorImage( curPrior );
  standGen->SetMeanAddHistogram( addJHCalc->GetMeanHistogram() );
  standGen->SetStdevAddHistogram( addJHCalc->GetStandardDeviationHistogram() );
  standGen->SetMeanSubHistogram( subJHCalc->GetMeanHistogram() );
  standGen->SetStdevSubHistogram( subJHCalc->GetStandardDeviationHistogram() );
  standGen->SetMeanNormHistogram( nomJHCalc->GetMeanHistogram() ); 
  standGen->SetStdevNormHistogram( nomJHCalc->GetStandardDeviationHistogram() );
  standGen->SetScale( scale );
  standGen->PrepFilter();

  patchGen->SetInputImage( curImage );
  patchGen->SetPriorImage( curPrior );
  patchGen->SetWidth( 5 );

  centerlineItr.GoToBegin();
  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename ImageType::PointType curPoint;
      typename CalculatorType::IndexType curIndex;
      typename CalculatorType::ContinuousIndexType curContIndex;

      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );

      // Set label value
      bool validPoint = true;
      std::string label;
      if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) &&
          addCalc->Intensity( curIndex ) &&
          curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( (int)(samplesAdd) != (int)(samplesAdd+sampleRateAdd) )
          {
          label = "addition";
          validPoint = true;
          }
        else
          {
          validPoint = false;
          }
        samplesAdd += sampleRateAdd;
        }
      else if( subtractions->TransformPhysicalPointToIndex( curPoint, 
                                                            curIndex ) &&
               subCalc->Intensity( curIndex ) &&
               curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( (int)(samplesSub) != (int)(samplesSub+sampleRateSub) )
          {
          label = "subtraction";
          validPoint = true;
          }
        else
          {
          validPoint = false;
          }
        samplesSub += sampleRateSub;
        }
      else if( curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( (int)(samplesNom) != (int)(samplesNom+sampleRateNom) )
          {
          label = "typical";
          validPoint = true;
          }
        else
          {
          validPoint = false;
          }
        samplesNom += sampleRateNom;
        }

      if( validPoint )
        {
        std::vector<double>::const_iterator itr;
        std::vector<double> standFeats = standGen->Evaluate( curPoint );
        std::vector<double> patchFeats = patchGen->Evaluate( curPoint );
        output << curIndex[0] << "," 
               << curIndex[1] << ",";
        for( itr = standFeats.begin(); itr != standFeats.end(); ++itr )
          {
          output << *itr << ",";
          }
        for( itr = patchFeats.begin(); itr != patchFeats.end(); ++itr )
          {
          output << *itr << ",";
          }
        output << label << "\n";
  
        if( count == 1000 )
          {
          progressReporter.Report( progress );
          count = 0;
          }
        else
          {
          ++count;
          }
        }
      } // end if centerline is at valid point
    progress += step;
    ++centerlineItr;
    }

  output.close();

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End( );
  
  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You May Need To update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
