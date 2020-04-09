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

#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>
#include <itkBinaryThresholdImageFilter.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/stats.hpp>

#include "SegmentUsingQuantileThresholdCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

#include "../CLI/tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::TimeProbesCollectorBase timeCollector;

  tube::CLIProgressReporter progressReporter(
    "SegmentUsingQuantileThreshold", CLPProcessInformation );
  progressReporter.Start();

  typedef TPixel                                PixelType;
  typedef itk::Image< PixelType, VDimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >     ReaderType;

  timeCollector.Start( "Load data" );
  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputVolume.c_str() );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Load data" );
  double progress = 0.1;
  progressReporter.Report( progress );

  typename ImageType::Pointer image = reader->GetOutput();

  typename ImageType::Pointer maskImage = NULL;
  if( !maskVolume.empty() )
    {
    timeCollector.Start( "Load mask" );
    typename ReaderType::Pointer maskReader = ReaderType::New();
    maskReader->SetFileName( maskVolume.c_str() );
    try
      {
      maskReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading mask: Exception caught: "
                          + std::string( err.GetDescription() ) );
      timeCollector.Report();
      return EXIT_FAILURE;
      }
    timeCollector.Stop( "Load mask" );
    progress = 0.2;
    progressReporter.Report( progress );

    maskImage = maskReader->GetOutput();
    }

  if( thresholdQuantile >= 0 && thresholdQuantile <= 1.0 )
    {
    timeCollector.Start( "Boost accumulate" );

    typedef boost::accumulators::accumulator_set< PixelType,
      boost::accumulators::stats<
        boost::accumulators::tag::p_square_quantile > >
      QuantileAccumulatorType;
    typedef itk::ImageRegionConstIterator< ImageType >
      ImageIteratorType;

    /*
     * Create a and configure a vector of length N of pointers
     * to BOOST accumulators -- Each of the N accumulators will
     * estimate exactly one of the given N desired quantile. If
     * the desired quantile is not within ( 0,1 ), throw an exception.
     */
    QuantileAccumulatorType acc( boost::accumulators::quantile_probability
      = thresholdQuantile );

    /*
     * Use an image iterator to iterate over all pixel/voxel and
     * and then add those values to the accumulators. Adding
     * the values will incrementally compute the quantile estimates.
     */
    ImageIteratorType imIt( image, image->GetLargestPossibleRegion() );
    if( !maskVolume.empty() )
      {
      ImageIteratorType maskIt( maskImage,
        maskImage->GetLargestPossibleRegion() );
      while( !imIt.IsAtEnd() )
        {
        PixelType p = imIt.Get();
        PixelType m = maskIt.Get();
        if( m != 0 )
          {
          acc( p );
          }
        ++imIt;
        ++maskIt;
        }
      }
    else
      {
      while( !imIt.IsAtEnd() )
        {
        PixelType p = imIt.Get();
        acc( p );
        ++imIt;
        }
      }

    PixelType qVal = boost::accumulators::p_square_quantile( acc );

    typedef itk::BinaryThresholdImageFilter< ImageType, ImageType >
      FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput( image );
    filter->SetLowerThreshold( qVal );
    filter->SetOutsideValue( 0 );
    filter->SetInsideValue( 1 );

    filter->Update();

    image = filter->GetOutput();

    timeCollector.Stop( "Boost accumulate" );
    }

  typedef itk::ImageFileWriter< ImageType  >   ImageWriterType;

  timeCollector.Start( "Save data" );
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( image );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string( err.GetDescription() ) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }
  timeCollector.Stop( "Save data" );
  progress = 1.0;
  progressReporter.Report( progress );
  progressReporter.End();

  timeCollector.Report();
  return EXIT_SUCCESS;
}

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  return tube::ParseArgsAndCallDoIt( inputVolume, argc, argv );
}
