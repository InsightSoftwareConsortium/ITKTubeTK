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

#include "itktubeMetaRidgeSeed.h"
#include "itktubeRidgeSeedSupervisedLinearBasisGenerator.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "SegmentTubeSeedsCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TImage >
void WriteRidgeSeed( const typename imageT::Pointer & img,
  std::string base, std::string ext, int num )
{
  typedef itk::ImageFileWriter< TImage >     RidgeSeedImageWriterType;

  typename RidgeSeedImageWriterType::Pointer rsImageWriter =
    RidgeSeedImageWriterType::New();
  std::string fname = base;
  char c[80];
  std::sprintf( c, ext.c_str(), num );
  fname += std::string( c );
  rsImageWriter->SetUseCompression( true );
  rsImageWriter->SetFileName( fname.c_str() );
  rsImageWriter->SetInput( img );
  rsImageWriter->Update();
}

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef TPixel                                      InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >    InputImageType;
  typedef itk::Image< unsigned short, VDimension >    MapImageType;
  typedef itk::Image< float, VDimension >             RidgeSeedImageType;

  typedef itk::ImageFileReader< RidgeSeedImageType >  ImageReaderType;
  typedef itk::ImageFileReader< MapImageType >        MapReaderType;
  typedef itk::ImageFileWriter< RidgeSeedImageType >  RidgeSeedImageWriterType;

  typedef itk::tube::RidgeSeedSupervisedLinearBasisGenerator< RidgeSeedImageType,
    MapImageType >  RidgeSeedGeneratorType;
  typename RidgeSeedGeneratorType::Pointer rsGenerator =
    RidgeSeedGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  reader = ImageReaderType::New();
  reader->SetFileName( inputImage.c_str() );
  reader->Update();
  rsGenerator->SetRidgeImage( reader->GetOutput() );

  timeCollector.Stop( "LoadData" );

  if( labelmap.size() > 0 )
    {
    timeCollector.Start( "LoadLabelMap" );
    typename MapReaderType::Pointer  inMapReader = MapReaderType::New();
    inMapReader->SetFileName( labelmap.c_str() );
    inMapReader->Update();
    rsGenerator->SetLabelMap( inMapReader->GetOutput() );
    rsGenerator->SetObjectId( objectId );
    timeCollector.Stop( "LoadLabelMap" );
    }

  if( loadTubeSeedInfo.size() > 0 )
    {
    timeCollector.Start( "LoadTubeSeed" );

    itk::tube::MetaRidgeSeed rsReader( loadTubeSeedInfo.c_str() );
    rsReader.Read();

    rsGenerator->SetBasisValues( rsReader.GetLDAValues() );
    rsGenerator->SetBasisMatrix( rsReader.GetLDAMatrix() );
    rsGenerator->SetWhitenMeans( rsReader.GetWhitenMeans() );
    rsGenerator->SetWhitenStdDevs( rsReader.GetWhitenStdDevs() );
    rsGenerator->SetScales( rsReader.GetRidgeSeedScales() );

    timeCollector.Stop( "LoadTubeSeed" );
    }
  else
    {
    if( labelmap.size() == 0 )
      {
      std::cerr << "Must specify a labelmap if training to find seeds"
                << std::endl;
      return EXIT_FAILURE;
      }

    rsGenerator->SetScales( ridgeScales );

    timeCollector.Start( "Update" );
    rsGenerator->Update();
    timeCollector.Stop( "Update" );
    rsGenerator->SetLabelMap( NULL );
    }

  timeCollector.Start( "SaveTubeSeedImage" );
  rsGenerator->UpdateBasisImages();
  typename RidgeSeedImageWriterType::Pointer rsImageWriter =
    RidgeSeedImageWriterType::New();
  rsImageWriter->SetUseCompression( true );
  rsImageWriter->SetFileName( outputSeedImage.c_str() );
  rsImageWriter->SetInput( rsGenerator->GetBasisImage( 0 ) );
  rsImageWriter->Update();
  timeCollector.Stop( "SaveTubeSeedImage" );

  if( saveTubeSeedInfo.size() > 0 )
    {
    timeCollector.Start( "SaveTubeSeedInfo" );
    itk::tube::MetaRidgeSeed rsWriter(
      rsGenerator->GetScales(),
      rsGenerator->GetBasisValues(),
      rsGenerator->GetBasisMatrix(),
      rsGenerator->GetWhitenMeans(),
      rsGenerator->GetWhitenStdDevs() );
    rsWriter.Write( saveTubeSeedInfo.c_str() );
    timeCollector.Stop( "SaveTubeSeedInfo" );
    }

  if( saveFeatureImages.size() > 0 )
    {
    timeCollector.Start( "SaveFeatureImages" );
    unsigned int numFeatures = rsGenerator->GetNumberOfFeatures();
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      WriteRidgeSeed< RidgeSeedImageType >(
        rsGenerator->GetFeatureImage( i ),
        saveFeatureImages, ".f%02d.mha", i );
      }
    timeCollector.Stop( "SaveFeatureImages" );
    }

  timeCollector.Report();

  return 0;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputImage, argc, argv );
}
