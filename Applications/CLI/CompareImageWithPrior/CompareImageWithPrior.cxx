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
#include "itkImageRegionConstIteratorWithIndex.h"

// The following three should be used in every CLI application
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Local Includes
#include "tubeSubImageGenerator.h"

// Includes specific to this CLI application
#include "itkRecursiveGaussianImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "CompareImageWithPriorCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "SampleCLIApplication",
                                                 CLPProcessInformation );
  progressReporter.Start();

  // typedefs for inputs
  typedef float                                              PixelType;
  typedef itk::Image< PixelType,  dimensionT >               ImageType;
  typedef itk::ImageFileReader< ImageType >                  ReaderType;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType > FullItrType;

  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  typename ReaderType::Pointer priorReader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  priorReader->SetFileName( priorVolume.c_str() );
  try
    {
    reader->Update();
    priorReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Load data");
  double progress = 0.1;
  progressReporter.Report( progress );

  typename ImageType::Pointer curImage = reader->GetOutput();
  typename ImageType::Pointer curPrior = priorReader->GetOutput();

  
  timeCollector.Start("Pull Sections");
  
  typename ImageType::RegionType region;
  typename ImageType::SizeType size;
  typename ImageType::IndexType start;
  
  std::vector<int> roiSize = std::vector<int>(dimensionT);
  
  size = curImage->GetLargestPossibleRegion().GetSize();

  for( unsigned int i = 0; i < dimensionT; ++i )
    {
    size[i] -= regionRadius;
    start[i] = regionRadius/2;
    roiSize[i] = regionRadius;
    }
  region.SetSize(size);
  region.SetIndex(start);

  typename ImageType::Pointer outImage = ImageType::New();
  outImage->SetRegions(region);
  outImage->Allocate();
  outImage->FillBuffer(0);
  
  FullItrType imageItr( curImage, region);
  imageItr.GoToBegin();
  while( !imageItr.IsAtEnd() )
    {
    bool doCalculation = true;
    bool useThreshold;
    bool useSkip;
    if( threshold != 0 )
      {
      useThreshold = true;
      useSkip = false;
      }
    else if( skipSize != -1 )
      {
      useThreshold = false;
      useSkip = true;
      }
    else
      {
      useThreshold = false;
      useSkip = false;
      }

    if( useThreshold && imageItr.Get() > threshold )
      {
      doCalculation = false;
      }

    typename ImageType::IndexType curIndex;
    std::vector<int> roiCenter;
    if( doCalculation )
      {
      curIndex = imageItr.GetIndex();    
      roiCenter = std::vector<int>(dimensionT);    
      for( unsigned int i = 0; i < dimensionT; ++i )
        {
        if( useSkip && (curIndex[i]-start[i]) % skipSize != 0 )
          {
          doCalculation = false;
          break;
          }
        roiCenter[i] = curIndex[i];
        }
      }

    if( doCalculation )
      {
      tube::SubImageGenerator<PixelType,dimensionT> subGenerator;
      subGenerator.SetRoiCenter(roiCenter);
      subGenerator.SetRoiSize(roiSize);
      subGenerator.SetInputVolume(curImage);
      subGenerator.SetInputMask(curPrior);
      subGenerator.Update();
      outImage->SetPixel(curIndex,1);
      }
    ++imageItr;
    }
  
  
  timeCollector.Stop("Pull Sections");


  typedef itk::ImageFileWriter< ImageType  > ImageWriterType;

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( outImage );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    std::cerr << "Exception caught: " << err << std::endl;
    return EXIT_FAILURE;
    }
  timeCollector.Stop("Save data");
  progress = 1.0;
  progressReporter.Report( progress );
  
  timeCollector.Report();

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
