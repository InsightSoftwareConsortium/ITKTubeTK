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

#include <itkTimeProbesCollectorBase.h>
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

// TubeTK includes
#include "tubeConvertShrunkenSeedImageToList.h"

// ITK includes
#include <itkCSVNumericObjectFileWriter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ConvertShrunkenSeedImageToListCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  // CLIProgressReporter is used to communicate progress with the Slicer GUI
  tube::CLIProgressReporter progressReporter( "Shrink Image",
    CLPProcessInformation );
  progressReporter.Start();

  typedef float                                     PixelType;
  typedef itk::Image< PixelType, VDimension >       ImageType;
  typedef itk::ImageFileReader< ImageType >         ImageReaderType;

  typedef itk::Vector< float, VDimension >          PointsPixelType;
  typedef itk::Image< PointsPixelType, VDimension > PointsImageType;
  typedef itk::ImageFileReader< PointsImageType >   PointsImageReaderType;

  typedef tube::ConvertShrunkenSeedImageToList
    < ImageType, PointsImageType > ConvertShrunkenSeedImageToListFilterType;
  typename ConvertShrunkenSeedImageToListFilterType::Pointer filter
    = ConvertShrunkenSeedImageToListFilterType::New();

  double progress = 0.1;
  progressReporter.Report( progress );

  timeCollector.Start( "Load image data" );
  typename ImageReaderType::Pointer inImageReader = ImageReaderType::New();
  inImageReader->SetFileName( inputShrunkenImageFileName.c_str() );
  try
    {
    inImageReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer inImage = inImageReader->GetOutput();
  timeCollector.Stop( "Load image data" );

  progress += 0.1;
  progressReporter.Report( progress );

  timeCollector.Start( "Load scale data" );
  typename ImageReaderType::Pointer inScaleReader = ImageReaderType::New();
  inScaleReader->SetFileName( inputShrunkenScaleImageFileName.c_str() );
  try
    {
    inScaleReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  typename ImageType::Pointer inScale = inScaleReader->GetOutput();
  timeCollector.Stop( "Load scale data" );

  progress += 0.1;
  progressReporter.Report( progress );

  timeCollector.Start( "Load point data" );

  typename PointsImageReaderType::Pointer inPointReader =
    PointsImageReaderType::New();

  inPointReader->SetFileName( inputShrunkenPointsImageFileName.c_str() );
  try
    {
    inPointReader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  typename PointsImageType::Pointer inPoint = inPointReader->GetOutput();
  timeCollector.Stop( "Load point data" );

  progress += 0.1;
  progressReporter.Report( progress );

  timeCollector.Start( "Generate output list" );

  filter->SetInput( inImage );
  filter->SetScaleImage( inScale );
  filter->SetPointsImage( inPoint );
  filter->SetThreshold( shrunkenImageThreshold );

  filter->Update();

  typedef vnl_matrix<PixelType> MatrixType;
  const unsigned int ARows =
    inImage->GetLargestPossibleRegion().GetNumberOfPixels();
  const unsigned int ACols = VDimension + 1;
  MatrixType matrix;
  matrix.set_size( ARows, ACols );
  matrix = filter->GetOutput();

  // write out the vnl_matrix object
  typedef itk::CSVNumericObjectFileWriter<PixelType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFieldDelimiterCharacter( ',' );
  writer->SetFileName( outputListFileName );
  writer->SetInput( &matrix );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject& exp )
    {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  timeCollector.Stop( "Generate output list" );

  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputShrunkenImageFileName.
  return tube::ParseArgsAndCallDoIt( inputShrunkenImageFileName, argc, argv );
}
