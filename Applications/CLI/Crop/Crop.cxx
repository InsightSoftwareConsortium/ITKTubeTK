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
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"
#include "tubeMessage.h"

// Includes specific to this CLI application
#include "itkCropImageFilter.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHleperFunctions
#include "CropCLP.h"

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
  tube::CLIProgressReporter    progressReporter( "Crop",
                                                 CLPProcessInformation );
  progressReporter.Start();

  typedef pixelT                                        PixelType;
  typedef itk::Image< PixelType,  dimensionT >          ImageType;
  typedef itk::ImageFileReader< ImageType >             ReaderType;
  
  timeCollector.Start("Load data");
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
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

  if( size.size() > 0 || max.size() > 0 || min.size() > 0 )
    {
    if( size.size() > 0 && max.size() > 0 )
      {
      tube::ErrorMessage( "You must specify either --size or --max options.  Not both." );
      return EXIT_FAILURE;
      }

    timeCollector.Start("Crop");

    typedef itk::CropImageFilter< ImageType, ImageType > FilterType;
    typename FilterType::Pointer filter;

    // Progress per iteration
    double progressFraction = 0.8/dimensionT;

    filter = FilterType::New();
    filter->SetInput( curImage );

    if( min.size() == 0 )
      {
      min.resize( dimensionT );
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        min[i] = 0;
        }
      }

    typename ImageType::SizeType imageSize;
    imageSize = curImage->GetLargestPossibleRegion().GetSize();
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      if( min[i] < 0 )
        {
        tube::ErrorMessage( "Min is less than 0." );
        return EXIT_FAILURE;
        }
      if( min[i] >= (int)(imageSize[i]) )
        {
        tube::ErrorMessage( "Min is larger than image size." );
        return EXIT_FAILURE;
        }
      }

    typename ImageType::SizeType outputSize;
    outputSize = curImage->GetLargestPossibleRegion().GetSize();
    if( size.size() > 0 )
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        outputSize[i] = size[i];
        if( outputSize[i] < 1 )
          {
          outputSize[i] = 1;
          }
        }
      }
    else
      {
      for( unsigned int i=0; i<dimensionT; i++ )
        {
        if( min[i] > max[i] )
          {
          int tf = min[i];
          min[i] = max[i];
          max[i] = tf;
          }
        outputSize[i] = max[i]-min[i]+1;
        }
      }
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      if( min[i] + outputSize[i] > imageSize[i] )
        {
        outputSize[i] = imageSize[i] - min[i];
        }
      }

    typename ImageType::SizeType lowerCropSize;
    typename ImageType::SizeType upperCropSize;
    for( unsigned int i=0; i<dimensionT; i++ )
      {
      lowerCropSize[i] = min[i];
      upperCropSize[i] = imageSize[i] - (min[i] + outputSize[i]);
      }

    filter->SetLowerBoundaryCropSize( lowerCropSize );
    filter->SetUpperBoundaryCropSize( upperCropSize );

    tube::CLIFilterWatcher( filter,
                            "Crop",
                            CLPProcessInformation,
                            progressFraction,
                            progress );

    filter->Update();
    curImage = filter->GetOutput();

    timeCollector.Stop("Crop");
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( curImage );
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
