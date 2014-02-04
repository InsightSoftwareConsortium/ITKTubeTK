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

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkTimeProbesCollectorBase.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageRegionIteratorWithIndex.h>

#include <itktubeRidgeExtractor.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ComputeSegmentTubesParametersCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class ImageT >
int WriteOutputImage( std::string & fileName, typename ImageT::Pointer
  & image )
{
  typedef itk::ImageFileWriter< ImageT  >  WriterType;

  typename WriterType::Pointer writer = WriterType::New();

  writer->SetInput( image );
  writer->SetFileName( fileName.c_str() );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing " + fileName + " : Exception caught: "
      + std::string(err.GetDescription()) );
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter    progressReporter(
    "ComputeSegmentTubesParameters", CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                    InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;

  typedef float                                     OutputPixelType;
  typedef itk::Image< OutputPixelType, VDimension > OutputImageType;

  typedef itk::tube::RidgeExtractor< InputImageType > RidgeFuncType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename InputImageType::Pointer inImage = reader->GetOutput();

  if( scale < inImage->GetSpacing()[0] * 0.3 )
    {
    scale = inImage->GetSpacing()[0] * 0.3;
    tubeWarningMacro( << "Reseting scale to " << scale );
    }

  timeCollector.Start("Compute ridgeness images");

  typename OutputImageType::RegionType region =
    inImage->GetLargestPossibleRegion();
  typename OutputImageType::SpacingType spacing =
    inImage->GetSpacing();
  if( supersample != 1 )
    {
    typename OutputImageType::RegionType::SizeType size = region.GetSize();
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      size[i] *= supersample;
      spacing[i] /= supersample;
      }
    region.SetSize( size );
    }

  typename OutputImageType::Pointer outImageRidgeness =
    OutputImageType::New();
  outImageRidgeness->CopyInformation( inImage );
  outImageRidgeness->SetSpacing( spacing );
  outImageRidgeness->SetRegions( region );
  outImageRidgeness->Allocate();

  typename OutputImageType::Pointer outImageRoundness =
    OutputImageType::New();
  outImageRoundness->CopyInformation( inImage );
  outImageRoundness->SetSpacing( spacing );
  outImageRoundness->SetRegions( region );
  outImageRoundness->Allocate();

  typename OutputImageType::Pointer outImageCurvature =
    OutputImageType::New();
  outImageCurvature->CopyInformation( inImage );
  outImageCurvature->SetSpacing( spacing );
  outImageCurvature->SetRegions( region );
  outImageCurvature->Allocate();

  typename OutputImageType::Pointer outImageLinearity =
    OutputImageType::New();
  outImageLinearity->CopyInformation( inImage );
  outImageLinearity->SetSpacing( spacing );
  outImageLinearity->SetRegions( region );
  outImageLinearity->Allocate();

  itk::ImageRegionIteratorWithIndex< OutputImageType > itR(
    outImageRidgeness, region );
  itk::ImageRegionIteratorWithIndex< OutputImageType > itO(
    outImageRoundness, region );
  itk::ImageRegionIteratorWithIndex< OutputImageType > itC(
    outImageCurvature, region );
  itk::ImageRegionIteratorWithIndex< OutputImageType > itL(
    outImageLinearity, region );

  typename RidgeFuncType::Pointer ridgeFunc = RidgeFuncType::New();
  ridgeFunc->SetInputImage( inImage );
  ridgeFunc->SetScale( scale );

  double supersampleFactor = 1.0 / supersample;
  double ridgeness = 0;
  double roundness = 0;
  double curvature = 0;
  double linearity = 0;
  typename RidgeFuncType::ContinuousIndexType cIndx;
  while( !itR.IsAtEnd() )
    {
    if( supersample != 1 )
      {
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        cIndx[i] = itR.GetIndex()[i] * supersampleFactor;
        }
      }
    else
      {
      for( unsigned int i = 0; i < VDimension; ++i )
        {
        cIndx[i] = itR.GetIndex()[i];
        }
      }
    ridgeness = ridgeFunc->Ridgeness( cIndx, roundness, curvature,
      linearity );
    itR.Set( (OutputPixelType) ridgeness );
    itO.Set( (OutputPixelType) roundness );
    itC.Set( (OutputPixelType) curvature );
    itL.Set( (OutputPixelType) linearity );
    ++itR;
    ++itO;
    ++itC;
    ++itL;
    }

  timeCollector.Stop("Compute ridgeness images");

  timeCollector.Start("Save data");

  std::string outName = outputImagesBaseFileName + ".ridge.mha";
  int result = WriteOutputImage< OutputImageType >( outName,
    outImageRidgeness );

  outName = outputImagesBaseFileName + ".round.mha";
  result += WriteOutputImage< OutputImageType >( outName,
    outImageRoundness );

  outName = outputImagesBaseFileName + ".curve.mha";
  result += WriteOutputImage< OutputImageType >( outName,
    outImageCurvature );

  outName = outputImagesBaseFileName + ".line.mha";
  result += WriteOutputImage< OutputImageType >( outName,
    outImageLinearity );

  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return result;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );
}
