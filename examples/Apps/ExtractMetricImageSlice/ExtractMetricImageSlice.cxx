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

#include "../CLI/tubeCLIFilterWatcher.h"
#include "tubeMessage.h"

#include <itkExtractImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "ExtractMetricImageSliceCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  const unsigned int MetricDimension = 6;
  const unsigned int SliceDimension = 2;
  typedef double     FloatType;

  typedef itk::Image< FloatType, MetricDimension > MetricImageType;
  typedef itk::Image< FloatType, SliceDimension >  SliceImageType;

  typedef itk::ImageFileReader< MetricImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputMetricImage );
  try
    {
    reader->UpdateOutputInformation();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading input image information: Exception caught: "
      + std::string( err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  typedef itk::ExtractImageFilter< MetricImageType, SliceImageType >
    ExtractFilterType;
  ExtractFilterType::Pointer extractor = ExtractFilterType::New();
  extractor->SetInput( reader->GetOutput() );
  extractor->SetDirectionCollapseToIdentity();

  MetricImageType::RegionType extractionRegion =
    reader->GetOutput()->GetLargestPossibleRegion();
  MetricImageType::SizeType extractionSize =
    extractionRegion.GetSize();
  MetricImageType::IndexType extractionIndex =
    extractionRegion.GetIndex();
  unsigned int indexCount = 0;
  for( int ii = 0; ii < static_cast< int >( MetricDimension ); ++ii )
    {
    if( ii != sliceDirections[0] && ii != sliceDirections[1] )
      {
      extractionSize[ii] = 0;
      extractionIndex[ii] = indices[indexCount];
      ++indexCount;
      }
    }

  extractionRegion.SetSize( extractionSize );
  extractionRegion.SetIndex( extractionIndex );
  extractor->SetExtractionRegion( extractionRegion );

  typedef itk::ImageFileWriter< SliceImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputMetricSlice );
  writer->SetInput( extractor->GetOutput() );

  tube::CLIFilterWatcher watcher( writer, "Writing Slice" );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing output image: Exception caught: "
      + std::string( err.GetDescription() ) );
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
