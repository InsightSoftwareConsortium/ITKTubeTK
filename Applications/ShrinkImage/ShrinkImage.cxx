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

#include <tubeShrinkWithBlendingImage.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ShrinkImageCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

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

  typedef float                                 PixelType;
  typedef itk::Image< PixelType, VDimension >   InputImageType;
  typedef itk::Image< PixelType, VDimension >   OutputImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< OutputImageType  > ImageWriterType;

  typedef itk::Vector< float, VDimension >         PointPixelType;
  typedef itk::Image< PointPixelType, VDimension > PointImageType;
  typedef itk::ImageFileWriter< PointImageType >   PointImageWriterType;

  typedef tube::ShrinkWithBlendingImage< InputImageType,
    OutputImageType > FilterType;

  timeCollector.Start("Load data");
  typename ImageReaderType::Pointer reader = ImageReaderType::New();
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

  typename InputImageType::Pointer inIm = reader->GetOutput();

  typename InputImageType::SpacingType     inSpacing = inIm->GetSpacing();
  typename InputImageType::PointType       inOrigin = inIm->GetOrigin();
  typename InputImageType::SizeType        inSize =
    inIm->GetLargestPossibleRegion().GetSize();
  typename InputImageType::IndexType       inIndex =
    inIm->GetLargestPossibleRegion().GetIndex();
  typename InputImageType::DirectionType   inDirection =
    inIm->GetDirection();

  typename OutputImageType::SizeType       outSize;
  typename OutputImageType::SpacingType    outSpacing;
  typename OutputImageType::PointType      outOrigin;
  typename OutputImageType::IndexType      outIndex;
  typename OutputImageType::DirectionType  outDirection;
  for( unsigned int i=0; i< VDimension; i++ )
    {
    outSpacing[i] = inSpacing[i];
    outOrigin[i] = inOrigin[i];
    outIndex[i] = inIndex[i];
    outSize[i] = inSize[i];
    }
  outDirection = inDirection;

  typename FilterType::Pointer filter = FilterType::New();
  filter->SetInput( inIm );

  typename FilterType::ShrinkFactorsType shrinkFactors;
  shrinkFactors.Fill( 1 );
  typename FilterType::InputSizeType newInputSize;
  newInputSize.Fill( 1 );

  typename FilterType::InputIndexType overlapVoxels;

  if( overlap.size() > 0 )
    {
    if( overlap.size() != VDimension )
      {
      std::cout << "Error: Size of overlap vector argument does not ";
      std::cout << "match the dimensionality of the image." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      overlapVoxels[ i ] = overlap[ i ];
      }
    filter->SetOverlap( overlapVoxels );
    }

  if( log )
    {
    filter->SetUseLog( true );
    }

  if( mean )
    {
    filter->SetBlendWithMax( false );
    filter->SetBlendWithMean( true );
    filter->SetBlendWithGaussianWeighting( false );
    }
  else if( gaussian )
    {
    filter->SetBlendWithMax( false );
    filter->SetBlendWithMean( false );
    filter->SetBlendWithGaussianWeighting( true );
    }
  else
    {
    filter->SetBlendWithMax( true );
    filter->SetBlendWithMean( false );
    filter->SetBlendWithGaussianWeighting( false );
    }

  if( divideBy.size() > 0 )
    {
    if( divideBy.size() != VDimension )
      {
      std::cout << "Error: Size of divideBy vector argument does not ";
      std::cout << "match the dimensionality of the image." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      shrinkFactors[ i ] = divideBy[ i ];
      }
      filter->SetShrinkFactors( shrinkFactors );
    }

  if( newSize.size() > 0 )
    {
    if( divideBy.size() > 0 )
      {
      std::cout << "Error: Specify only one of divideBy or newSize ";
      std::cout << "arguments." << std::endl;
      return EXIT_FAILURE;
      }
    if( newSize.size() != VDimension )
      {
      std::cout << "Error: Size of newSize vector argument does not ";
      std::cout << "match the dimensionality of the image." << std::endl;
      return EXIT_FAILURE;
      }
    for( unsigned int i = 0; i < VDimension; ++i )
      {
      newInputSize[ i ] = newSize[ i ];
      }
    filter->SetNewSize(newInputSize);
    }

  tube::CLIFilterWatcher watcher( filter, "Shrink Filter",
    CLPProcessInformation, progress, 0.8, true );

  filter->Update();

  timeCollector.Start("Save data");
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetFileName( outputImageFileName.c_str() );
  writer->SetInput( filter->GetOutput() );
  writer->SetUseCompression( true );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Writing volume: Exception caught: "
      + std::string(err.GetDescription()) );
    timeCollector.Report();
    return EXIT_FAILURE;
    }

  if( !mipPointImageFileName.empty() && !mean && !gaussian )
    {
    typename PointImageWriterType::Pointer mipPointImageWriter =
      PointImageWriterType::New();
    mipPointImageWriter->SetFileName( mipPointImageFileName );
    mipPointImageWriter->SetUseCompression( true );
    mipPointImageWriter->SetInput( filter->GetPointImage() );
    try
      {
      mipPointImageWriter->Update();
      }
    catch( itk::ExceptionObject & e )
      {
      std::cout << "Exception caught during image index write:" << std::endl
        << e << std::endl;
      return EXIT_FAILURE;
      }
    }

  timeCollector.Stop("Save data");
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
  //   you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );
}
