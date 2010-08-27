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
#include "itkNJetImageFunction.h"
#include "itkJointHistogramImageFunction.h"
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
#include "GenerateReducedFeaturesForWekaCLP.h"

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
  tube::CLIProgressReporter progressReporter( "GenerateReducedFeaturesForWeka",
                                              CLPProcessInformation );
  progressReporter.Start();

  // typedefs for data structures
  typedef float                                                 PixelType;
  typedef itk::OrientedImage<PixelType, dimensionT>             ImageType;
  typedef itk::ImageFileReader<ImageType>                       ReaderType;
  typedef itk::ImageFileWriter<ImageType>                       WriterType;

  // typedefs for numerics
  typedef itk::NJetImageFunction<ImageType>                     CalculatorType;
  typedef itk::JointHistogramImageFunction<ImageType>           HistCalcType;

  // typedefs for filters
  typedef itk::RescaleIntensityImageFilter<ImageType,ImageType> RescaleType;

  // typedefs for iterators
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType>     IterType;
  typedef itk::NeighborhoodIterator<ImageType>                  NeighborIterType;

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

  // Setup the Calculators
  typename CalculatorType::Pointer inputCalc = CalculatorType::New();
  inputCalc->SetInputImage( curImage );
  inputCalc->SetInverseRidgeness( true );

  typename CalculatorType::Pointer priorCalc = CalculatorType::New();
  priorCalc->SetInputImage( curPrior );
  priorCalc->SetInverseRidgeness( true );

  typename CalculatorType::Pointer addCalc = CalculatorType::New();
  addCalc->SetInputImage( additions );

  typename CalculatorType::Pointer subCalc = CalculatorType::New();
  subCalc->SetInputImage( subtractions );

  // Setup the Joint-Histogram Calculators
  typename HistCalcType::Pointer addJHCalc = HistCalcType::New();
  addJHCalc->SetInputImage( curImage );
  addJHCalc->SetInputMask( curPrior );
  typename HistCalcType::Pointer subJHCalc = HistCalcType::New();
  subJHCalc->SetInputImage( curImage );
  subJHCalc->SetInputMask( curPrior );
  typename HistCalcType::Pointer nomJHCalc = HistCalcType::New();
  nomJHCalc->SetInputImage( curImage );
  nomJHCalc->SetInputMask( curPrior );

  // Setup the sigmas based on the wire scale
  PixelType sigmaMedium = (scale/2)*0.6667;
  PixelType sigmaSmall = 0.6667*sigmaMedium;
  PixelType sigmaLarge = 1.3333*sigmaMedium;

  double samplesNom = 0;
  double samplesAdd = 0;
  double samplesSub = 0;

  IterType centerlineItr( centerlines, 
                          centerlines->GetLargestPossibleRegion() );

  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename CalculatorType::PointType curPoint;
      typename ImageType::IndexType curIndex;
  
      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );
  
      if( curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) 
            && addCalc->Evaluate( curPoint, sigmaSmall ) )
          {
          ++samplesAdd;
          }
        else if( subtractions->TransformPhysicalPointToIndex( curPoint,
                                                              curIndex ) 
                 && subCalc->Evaluate( curPoint, sigmaSmall ) )
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

  std::cout << "Add = " << samplesAdd << std::endl;
  std::cout << "Sub = " << samplesSub << std::endl;
  std::cout << "Nom = " << samplesNom << std::endl;

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
      typename CalculatorType::PointType curPoint;
      typename ImageType::IndexType curIndex;
      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );
      
      if( curImage->TransformPhysicalPointToIndex( curPoint, curIndex ) )
        {
        if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) 
            && addCalc->Evaluate( curPoint, sigmaSmall ) )
          {
          if( (int)(samplesAdd) != (int)(samplesAdd+sampleRateAdd) )
            {
            addJHCalc->Precompute( curPoint );
            }
          samplesAdd += sampleRateAdd;
          }
        else if( subtractions->TransformPhysicalPointToIndex( curPoint, 
                                                            curIndex ) 
                 && subCalc->Evaluate( curPoint, sigmaSmall ) )
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
  typedef typename HistCalcType::HistogramType    HistrogramType;
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
  
  // Get the additional features ready
  std::vector<std::string> featureNames;
  std::vector<std::string> featureFilenames;
  std::vector<typename ImageType::Pointer> featureImages;
  std::vector<std::string>::const_iterator featureItr;
  bool isName = true;
  unsigned int featureCount = 0;
  for( featureItr = otherFeatures.begin(); featureItr != otherFeatures.end();
       ++featureItr )
    {
    std::string val = *featureItr;
    if( isName )
      {
      featureNames.push_back( val );
      }
    else
      {
      featureFilenames.push_back( val );
      typename ReaderType::Pointer featureReader = ReaderType::New();
      featureReader->SetFileName( val.c_str() );
      try
        {
        featureReader->Update();
        }
      catch( itk::ExceptionObject & err )
        {
        std::cerr << "ExceptionObject caught while reading feature "
                  << featureCount << std::endl;
        std::cerr << err << std::endl;
        return EXIT_FAILURE;
        }
      featureImages.push_back( featureReader->GetOutput() );
      ++featureCount;
      }
    isName = !isName;
    }


  // Get the outfile started
  std::ofstream output( outputFile.c_str() );
  output << "% Generated by GenerateReducedFeaturesForWeka in TubeTK\n";
  output << "@RELATION errors\n";
  output << "\n";
  output << "@ATTRIBUTE x NUMERIC\n";
  output << "@ATTRIBUTE y NUMERIC\n";

  output << "@ATTRIBUTE v1sg-v1se NUMERIC\n";
  output << "@ATTRIBUTE v1mg-v1me NUMERIC\n";
  output << "@ATTRIBUTE v1lg-v1le NUMERIC\n";

  output << "@ATTRIBUTE v2g-v2e NUMERIC\n";
  output << "@ATTRIBUTE v3g-v3e NUMERIC\n";

  output << "@ATTRIBUTE Z_add NUMERIC\n";
  output << "@ATTRIBUTE Z_sub NUMERIC\n";
  output << "@ATTRIBUTE Z_nom NUMERIC\n";

  // Patch Based - begin

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

  // Patch Based - END

  for( featureItr = featureNames.begin(); featureItr != featureNames.end();
       ++featureItr )
    {
    output << "@ATTRIBUTE " << *featureItr << " NUMERIC\n";
    }
  
  output << "@ATTRIBUTE class {typical,addition,subtraction}\n";
  output << "\n";
  output << "@DATA\n";


  samplesNom = 0;
  samplesAdd = 0;
  samplesSub = 0;
  double portion = 0.9;
  double step = portion/samplesTotal;
  unsigned int count = 0;

  typename NeighborIterType::SizeType radius  = {{ceil(scale/2),
                                                  ceil(scale/2)}};
  NeighborIterType imageItr(radius,
                            curImage,
                            curImage->GetLargestPossibleRegion() );

  NeighborIterType priorItr(radius,
                            curPrior,
                            curPrior->GetLargestPossibleRegion() );

  centerlineItr.GoToBegin();
  imageItr.GoToBegin();
  priorItr.GoToBegin();
  while( !centerlineItr.IsAtEnd() )
    {
    if( centerlineItr.Get() > 0 )
      {
      typename CalculatorType::PointType curPoint;
      typename ImageType::IndexType curIndex;
      centerlines->TransformIndexToPhysicalPoint( centerlineItr.GetIndex(),
                                                  curPoint );

      // Set label value
      bool validPoint = true;
      std::string label;
      if( additions->TransformPhysicalPointToIndex( curPoint, curIndex ) &&
          addCalc->Evaluate( curPoint, sigmaSmall ) &&
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
               subCalc->Evaluate( curPoint, sigmaSmall ) &&
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
        
        // holders for features
        PixelType v1sg, v1mg, v1lg, v2g, v3g, v1se, v1me, v1le, v2e, v3e;
        double zAdd, zSub, zNom;
        std::vector<PixelType> extras;
  
        v1sg = priorCalc->Ridgeness( curPoint, sigmaSmall );
        v1mg = priorCalc->Ridgeness( curPoint, sigmaMedium );
        v1lg = priorCalc->Ridgeness( curPoint, sigmaLarge );

        v1se = inputCalc->Ridgeness( curPoint, sigmaSmall );
        v1me = inputCalc->Ridgeness( curPoint, sigmaMedium );
        v1le = inputCalc->Ridgeness( curPoint, sigmaLarge );
  
        double pS = priorCalc->Evaluate( curPoint, sigmaSmall );
        double pM = priorCalc->Evaluate( curPoint, sigmaMedium );
        double pL = priorCalc->Evaluate( curPoint, sigmaLarge );

        double iS = inputCalc->Evaluate( curPoint, sigmaSmall );
        double iM = inputCalc->Evaluate( curPoint, sigmaMedium );
        double iL = inputCalc->Evaluate( curPoint, sigmaLarge );

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

        typename std::vector<typename ImageType::Pointer>::const_iterator 
          featureImageItr;
        for( featureImageItr = featureImages.begin(); 
             featureImageItr != featureImages.end();
             ++featureImageItr )
          {
          typename ImageType::IndexType featIndex; 
          (*featureImageItr)->TransformPhysicalPointToIndex( curPoint, 
                                                             featIndex );
          PixelType featPix = (*featureImageItr)->GetPixel( featIndex );
          extras.push_back( featPix );
          }
  
        // begin patch features
        typename NeighborIterType::NeighborhoodType imageN = imageItr.GetNeighborhood();
        typename NeighborIterType::NeighborhoodType priorN = priorItr.GetNeighborhood();

        assert(imageN.Size() == priorN.Size());
        assert(std::distance(imageN.Begin(), imageN.End()) == imageN.Size());

        std::vector<float> normalizedImagePatch(imageN.Size());
        std::vector<float> normalizedPriorPatch(imageN.Size());

        // compute the mean intensity of the two patches
        float imageMean = std::accumulate(imageN.Begin(), imageN.End(), 0.0) / imageN.Size();
        float priorMean = std::accumulate(priorN.Begin(), priorN.End(), 0.0) / imageN.Size();
      
        std::transform(imageN.Begin(), imageN.End(),
                       normalizedImagePatch.begin(),
                       std::bind2nd(std::minus<float>(), imageMean));

        std::transform(priorN.Begin(), priorN.End(),
                       normalizedPriorPatch.begin(),
                       std::bind2nd(std::minus<float>(), priorMean));                    
      
        float imageStdDev = std::sqrt(std::inner_product(normalizedImagePatch.begin(), normalizedImagePatch.end(),
                                                         normalizedImagePatch.begin(),
                                                         0.0) /
                                      (imageN.Size() - 1));

        float priorStdDev = std::sqrt(std::inner_product(normalizedPriorPatch.begin(), normalizedPriorPatch.end(),
                                                         normalizedPriorPatch.begin(),
                                                         0.0) /
                                      (priorN.Size() - 1));

        float imageNorm = std::sqrt(std::inner_product(normalizedImagePatch.begin(), normalizedImagePatch.end(),
                                                       normalizedImagePatch.begin(), 0.0));

        float priorNorm = std::sqrt(std::inner_product(normalizedPriorPatch.begin(), normalizedPriorPatch.end(),
                                                       normalizedPriorPatch.begin(), 0.0));

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

        float crossCorrelation = std::inner_product(normalizedImagePatch.begin(), normalizedImagePatch.end(),
                                                    normalizedPriorPatch.begin(), 0.0);

        std::vector<float> differencePatch(imageN.Size());

        std::transform(imageN.Begin(), imageN.End(),
                       priorN.Begin(),
                       differencePatch.begin(),
                       std::minus<float>());

        float norm = std::inner_product(differencePatch.begin(), differencePatch.end(),
                                        differencePatch.begin(), 0.0);

        float total = std::accumulate(differencePatch.begin(), differencePatch.end(), 0.0);

        float maxdiff = *std::max_element(differencePatch.begin(),
                                          differencePatch.end());

        float mindiff = *std::min_element(differencePatch.begin(),
                                          differencePatch.end());

        size_t q1 = static_cast<size_t>(.25*differencePatch.size());
        size_t q2 = static_cast<size_t>(.50*differencePatch.size());
        size_t q3 = static_cast<size_t>(.75*differencePatch.size());
        size_t q95 = static_cast<size_t>(.95*differencePatch.size());

        std::sort(differencePatch.begin(), differencePatch.end());

        float q1val = differencePatch[q1];
        float q2val = differencePatch[q2];
        float q3val = differencePatch[q3];
        float q95val = differencePatch[q95];

        // end patch features

        output << curIndex[0] << "," << curIndex[1] << ","
               << v1sg-v1se << "," << v1mg-v1me << "," << v1lg-v1le << "," 
               << v2g-v2e << "," << v3g-v3e << ","
               << zAdd << "," << zSub << "," << zNom << ","
               << imageMean-priorMean << "," 
               << imageStdDev-priorStdDev << ","
               << crossCorrelation << ","
               << mindiff << ","
               << q1val << ","
               << q2val << ","
               << q3val << ","
               << q95val << ","
               << maxdiff << ","
               << total << ","
               << norm << ",";
        std::vector<PixelType>::const_iterator extraItr;
        for( extraItr = extras.begin(); 
             extraItr != extras.end(); ++extraItr )
          {
          output << *extraItr << ",";
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
    ++priorItr;
    ++imageItr;
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
