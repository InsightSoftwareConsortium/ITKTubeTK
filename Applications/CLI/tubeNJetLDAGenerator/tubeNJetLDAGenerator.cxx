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
#if defined( _MSC_VER )
#pragma warning ( disable : 4786 )
#endif


#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// The following four should be used in every CLI application
#include "tubeMessage.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "itkTimeProbesCollectorBase.h"

// Includes specific to this CLI application
#include "tubeStringUtilities.h"
#include "itkTubeNJetLDAGenerator.h"
#include "itkTubeMetaNJetLDA.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "tubeNJetLDAGeneratorCLP.h"

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
  sprintf( c, ext.c_str(), num );
  fname += std::string( c );
  ldaImageWriter->SetFileName( fname.c_str() );
  ldaImageWriter->SetInput( img );
  ldaImageWriter->Update();
}

template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef pixelT                                   InputPixelType;
  typedef itk::OrientedImage< InputPixelType, dimensionT >
                                                   InputImageType;
  typedef itk::OrientedImage< unsigned short, dimensionT >
                                                   MaskImageType;
  typedef itk::OrientedImage< float, dimensionT >  LDAImageType;

  typedef itk::ImageFileReader< LDAImageType >     ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< LDAImageType >     LDAImageWriterType;

  typedef itk::tube::NJetLDAGenerator< LDAImageType, MaskImageType >
    LDAGeneratorType;
  typename LDAGeneratorType::Pointer ldaGenerator = LDAGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int vNum=0; vNum<inputVolumesList.size(); vNum++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[vNum].c_str() );
    reader->Update();
    if( vNum == 0 )
      {
      ldaGenerator->SetFeatureImage( reader->GetOutput() );
      }
    else
      {
      ldaGenerator->AddFeatureImage( reader->GetOutput() );
      }
    }

  typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
  inMaskReader->SetFileName( labelmap.c_str() );
  inMaskReader->Update();
  ldaGenerator->SetLabelmap( inMaskReader->GetOutput() );

  timeCollector.Stop( "LoadData" );

  if( objectId.size() > 0 )
    {
    ldaGenerator->SetObjectId( objectId[0] );
    if( objectId.size() > 1 )
      {
      for( unsigned int o=1; o<objectId.size(); o++ )
        {
        ldaGenerator->AddObjectId( objectId[o] );
        }
      }
    }

  if( usePCA )
    {
    ldaGenerator->SetPerformPCA( true );
    }

  if( loadLDAInfo.size() > 0 )
    {
    timeCollector.Start( "LoadLDA" );

    itk::tube::MetaNJetLDA ldaReader( loadLDAInfo.c_str() );
    ldaReader.Read();

    ldaGenerator->SetLDAValues( ldaReader.GetLDAValues() );
    ldaGenerator->SetLDAMatrix( ldaReader.GetLDAMatrix() );
    ldaGenerator->SetZeroScales( ldaReader.GetZeroScales() );
    ldaGenerator->SetFirstScales( ldaReader.GetFirstScales() );
    ldaGenerator->SetSecondScales( ldaReader.GetSecondScales() );
    ldaGenerator->SetRidgeScales( ldaReader.GetRidgeScales() );

    timeCollector.Stop( "LoadLDA" );
    }
  else
    {
    timeCollector.Start( "Update" );

    ldaGenerator->SetZeroScales( zeroScales );
    ldaGenerator->SetFirstScales( firstScales );
    ldaGenerator->SetSecondScales( secondScales );
    ldaGenerator->SetRidgeScales( ridgeScales );

    ldaGenerator->SetForceIntensityConsistency( forceSign );
    ldaGenerator->SetForceOrientationInsensitivity( forceSymmetry );

    ldaGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  unsigned int numLDA = ldaGenerator->GetNumberOfLDA();
  if( useNumberOfLDA>0 && useNumberOfLDA < (int)numLDA )
    {
    numLDA = useNumberOfLDA;
    }

  if( outputBase.size() > 0 )
    {
    timeCollector.Start( "SaveLDAImages" );

    ldaGenerator->UpdateLDAImages();

    for( unsigned int i=0; i<numLDA; i++ )
      {
      typename LDAImageWriterType::Pointer ldaImageWriter =
        LDAImageWriterType::New();
      std::string fname = outputBase;
      char c[80];
      sprintf( c, ".lda%02d.mha", i );
      fname += std::string( c );
      ldaImageWriter->SetFileName( fname.c_str() );
      ldaImageWriter->SetInput( ldaGenerator->GetLDAImage( i ) );
      ldaImageWriter->Update();
      }
    timeCollector.Stop( "SaveLDAImages" );
    }

  if( saveLDAInfo.size() > 0 )
    {
    timeCollector.Start( "SaveLDA" );
    itk::tube::MetaNJetLDA ldaWriter(
      ldaGenerator->GetZeroScales(),
      ldaGenerator->GetFirstScales(),
      ldaGenerator->GetSecondScales(),
      ldaGenerator->GetRidgeScales(),
      ldaGenerator->GetLDAValues(),
      ldaGenerator->GetLDAMatrix() );
    ldaWriter.Write( saveLDAInfo.c_str() );
    timeCollector.Stop( "SaveLDA" );
    }

  if( saveFeatureImages.size() > 0 )
    {
    unsigned int numFeatures = ldaGenerator->GetNumberOfFeatures();
    for( unsigned int i=0; i<numFeatures; i++ )
      {
      WriteLDA< LDAImageType >( ldaGenerator->GetNJetFeatureImage( i ),
        saveFeatureImages, ".f%02d.mha", i );
      }
    }

  if( saveKernelImages.size() > 0 )
    {
    for( unsigned int i=0; i<numLDA; i++ )
      {
      WriteLDA< LDAImageType >( ldaGenerator->GetNJetKernelImage( i ),
        saveKernelImages, ".k%02d.mha", i );
      }
    }

  timeCollector.Report();

  return 0;
}

// Main
int main( int argc, char **argv )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolumesList[0], argc, argv );
}
