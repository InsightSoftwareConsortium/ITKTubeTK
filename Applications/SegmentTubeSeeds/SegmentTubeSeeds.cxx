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

#include "itktubeRidgeSeedFilterIO.h"
#include "itktubeRidgeSeedFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "SegmentTubeSeedsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TImage >
void WriteImageInSequence( const typename TImage::Pointer & img,
  const std::string & base, const std::string & ext, int num )
{
  typedef itk::ImageFileWriter< TImage >     ImageWriterType;

  typename ImageWriterType::Pointer rsImageWriter = ImageWriterType::New();
  std::string fname = base;
  char c[80];
  std::sprintf( c, ext.c_str(), num );
  fname += std::string( c );
  rsImageWriter->SetUseCompression( true );
  rsImageWriter->SetFileName( fname.c_str() );
  rsImageWriter->SetInput( img );
  rsImageWriter->Update();
}

template< class TImage >
void WriteImage( const typename TImage::Pointer & img,
  const std::string & str )
{
  typedef itk::ImageFileWriter< TImage >     ImageWriterType;

  typename ImageWriterType::Pointer rsImageWriter = ImageWriterType::New();
  rsImageWriter->SetUseCompression( true );
  rsImageWriter->SetFileName( str.c_str() );
  rsImageWriter->SetInput( img );
  rsImageWriter->Update();
}

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.:w
  //
  itk::TimeProbesCollectorBase timeCollector;

  typedef TPixel                                      InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >    InputImageType;
  typedef itk::Image< unsigned short, VDimension >    LabelMapImageType;
  typedef itk::Image< float, VDimension >             RidgeSeedImageType;

  typedef itk::ImageFileReader< RidgeSeedImageType >  ImageReaderType;
  typedef itk::ImageFileReader< LabelMapImageType >   LabelMapReaderType;

  typedef itk::tube::RidgeSeedFilter< RidgeSeedImageType,
    LabelMapImageType >   RidgeSeedGeneratorType;
  typename RidgeSeedGeneratorType::Pointer rsGenerator =
    RidgeSeedGeneratorType::New();

  typedef itk::tube::RidgeSeedFilterIO< RidgeSeedImageType,
    LabelMapImageType >   RidgeSeedFilterIOType;

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  reader = ImageReaderType::New();
  reader->SetFileName( inputImage.c_str() );
  reader->Update();
  rsGenerator->SetInput( reader->GetOutput() );

  timeCollector.Stop( "LoadData" );

  if( !labelMap.empty() )
    {
    timeCollector.Start( "LoadLabelMap" );
    typename LabelMapReaderType::Pointer  inMapReader =
      LabelMapReaderType::New();
    inMapReader->SetFileName( labelMap.c_str() );
    inMapReader->Update();
    rsGenerator->SetLabelMap( inMapReader->GetOutput() );
    rsGenerator->SetRidgeId( tubeId );
    rsGenerator->SetBackgroundId( backgroundId );
    rsGenerator->SetUnknownId( unknownId );
    timeCollector.Stop( "LoadLabelMap" );
    }

  if( !loadTubeSeedInfo.empty() )
    {
    timeCollector.Start( "LoadTubeSeed" );

    RidgeSeedFilterIOType rsReader( rsGenerator );
    rsReader.Read( loadTubeSeedInfo.c_str() );

    timeCollector.Stop( "LoadTubeSeed" );

    if( !saveTubeSeedInfo.empty() )
      {
      timeCollector.Start( "SaveTubeSeedInfo" );
      RidgeSeedFilterIOType rsWriter( rsGenerator );
      rsWriter.Write( saveTubeSeedInfo.c_str() );
      timeCollector.Stop( "SaveTubeSeedInfo" );
      }
    }
  else
    {
    if( labelMap.empty() )
      {
      std::cerr << "Must specify a labelMap if training to find seeds"
                << std::endl;
      return EXIT_FAILURE;
      }

    rsGenerator->SetScales( ridgeScales );

    timeCollector.Start( "Update" );
    rsGenerator->Update();
    timeCollector.Stop( "Update" );
    rsGenerator->SetLabelMap( NULL );
    }

  rsGenerator->SetSeedTolerance( seedTolerance );

  timeCollector.Start( "SaveTubeSeedImage" );

  rsGenerator->ClassifyImages();
  WriteImage< LabelMapImageType >( rsGenerator->GetOutput(),
    outputSeedImage );
  timeCollector.Stop( "SaveTubeSeedImage" );

  if( outputSeedScaleImage.size() > 0 )
    {
    timeCollector.Start( "SaveTubeSeedScaleImage" );
    WriteImage< RidgeSeedImageType >( rsGenerator->GetOutputSeedScales(),
      outputSeedScaleImage );
    timeCollector.Stop( "SaveTubeSeedScaleImage" );
    }

  if( saveTubeSeedInfo.size() > 0 )
    {
    timeCollector.Start( "SaveTubeSeedInfo" );
    RidgeSeedFilterIOType rsWriter( rsGenerator );
    rsWriter.Write( saveTubeSeedInfo.c_str() );
    timeCollector.Stop( "SaveTubeSeedInfo" );
    }

  if( saveFeatureImages.size() > 0 )
    {
    timeCollector.Start( "SaveFeatureImages" );
    unsigned int numFeatures = rsGenerator->GetNumberOfBasis();
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      WriteImageInSequence< RidgeSeedImageType >(
        rsGenerator->GetBasisImage( i ),
        saveFeatureImages, ".f%02d.mha", i );
      }
    timeCollector.Stop( "SaveFeatureImages" );
    }

  timeCollector.Report();

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputImage, argc, argv );
}
