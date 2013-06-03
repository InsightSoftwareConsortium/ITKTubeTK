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

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// The following four should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include <itkTimeProbesCollectorBase.h>

// Includes specific to this CLI application
#include "tubeStringUtilities.h"
#include "itkTubeNJetLDAGenerator.h"
#include "itkTubeMetaNJetLDA.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "SegmentTubeSeedsCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

template < class imageT >
void WriteLDA( const typename imageT::Pointer & img,
  std::string base, std::string ext, int num )
{
  typedef itk::ImageFileWriter< imageT >     LDAImageWriterType;

  typename LDAImageWriterType::Pointer ldaImageWriter =
    LDAImageWriterType::New();
  std::string fname = base;
  char c[80];
  std::sprintf( c, ext.c_str(), num );
  fname += std::string( c );
  ldaImageWriter->SetUseCompression( true );
  ldaImageWriter->SetFileName( fname.c_str() );
  ldaImageWriter->SetInput( img );
  ldaImageWriter->Update();
}

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef pixelT                                   InputPixelType;
  typedef itk::Image< InputPixelType, dimensionT > InputImageType;
  typedef itk::Image< unsigned short, dimensionT > MapImageType;
  typedef itk::Image< float, dimensionT >          LDAImageType;

  typedef itk::ImageFileReader< LDAImageType >     ImageReaderType;
  typedef itk::ImageFileReader< MapImageType >     MapReaderType;
  typedef itk::ImageFileWriter< LDAImageType >     LDAImageWriterType;

  typedef itk::tube::NJetLDAGenerator< LDAImageType, MapImageType >
    LDAGeneratorType;
  typename LDAGeneratorType::Pointer ldaGenerator = LDAGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  reader = ImageReaderType::New();
  reader->SetFileName( inputImage.c_str() );
  reader->Update();
  ldaGenerator->SetFeatureImage( reader->GetOutput() );

  timeCollector.Stop( "LoadData" );

  if( labelmap.size() > 0 )
    {
    timeCollector.Start( "LoadLabelMap" );
    typename MapReaderType::Pointer  inMapReader = MapReaderType::New();
    inMapReader->SetFileName( labelmap.c_str() );
    inMapReader->Update();
    ldaGenerator->SetLabelmap( inMapReader->GetOutput() );
    ldaGenerator->SetObjectId( objectId );
    timeCollector.Stop( "LoadLabelMap" );
    }

  if( loadVesselSeedInfo.size() > 0 )
    {
    timeCollector.Start( "LoadVesselSeed" );

    itk::tube::MetaNJetLDA ldaReader( loadVesselSeedInfo.c_str() );
    ldaReader.Read();

    ldaGenerator->SetLDAValues( ldaReader.GetLDAValues() );
    ldaGenerator->SetLDAMatrix( ldaReader.GetLDAMatrix() );
    ldaGenerator->SetZeroScales( ldaReader.GetZeroScales() );
    ldaGenerator->SetFirstScales( ldaReader.GetFirstScales() );
    ldaGenerator->SetSecondScales( ldaReader.GetSecondScales() );
    ldaGenerator->SetRidgeScales( ldaReader.GetRidgeScales() );

    timeCollector.Stop( "LoadVesselSeed" );
    }
  else
    {
    if( labelmap.size() == 0 )
      {
      std::cerr << "Must specify a labelmap if training to find seeds"
                << std::endl;
      return EXIT_FAILURE;
      }


    int numScales = ridgeScales.size();

    typename LDAGeneratorType::NJetScalesType midScale(1);
    midScale[0] = ridgeScales[(int)(numScales/2)];
    ldaGenerator->SetZeroScales( midScale );
    ldaGenerator->SetFirstScales( midScale );
    ldaGenerator->SetRidgeScales( ridgeScales );

    ldaGenerator->SetForceOrientationInsensitivity( true );

    ldaGenerator->Update();
    timeCollector.Stop( "Update" );
    }

  timeCollector.Start( "SaveVesselSeedImage" );
  ldaGenerator->UpdateLDAImages();
  typename LDAImageWriterType::Pointer ldaImageWriter =
    LDAImageWriterType::New();
  ldaImageWriter->SetUseCompression( true );
  ldaImageWriter->SetFileName( outputSeedImage.c_str() );
  ldaImageWriter->SetInput( ldaGenerator->GetLDAImage( 0 ) );
  ldaImageWriter->Update();
  timeCollector.Stop( "SaveVesselSeedImage" );

  if( saveVesselSeedInfo.size() > 0 )
    {
    timeCollector.Start( "SaveVesselSeedInfo" );
    itk::tube::MetaNJetLDA ldaWriter(
      ldaGenerator->GetZeroScales(),
      ldaGenerator->GetFirstScales(),
      ldaGenerator->GetSecondScales(),
      ldaGenerator->GetRidgeScales(),
      ldaGenerator->GetLDAValues(),
      ldaGenerator->GetLDAMatrix() );
    ldaWriter.Write( saveVesselSeedInfo.c_str() );
    timeCollector.Stop( "SaveVesselSeedInfo" );
    }

  if( saveFeatureImages.size() > 0 )
    {
    timeCollector.Start( "SaveFeatureImages" );
    unsigned int numFeatures = ldaGenerator->GetNumberOfFeatures();
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      WriteLDA< LDAImageType >( ldaGenerator->GetNJetFeatureImage( i ),
        saveFeatureImages, ".f%02d.mha", i );
      }
    timeCollector.Stop( "SaveFeatureImages" );
    }

  timeCollector.Report();

  return 0;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputImage, argc, argv );
}
