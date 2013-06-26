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

#include "itktubeMetaNJetLDA.h"
#include "itktubeNJetSupervisedLinearBasisGenerator.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "EnhanceUsingNJetDiscriminantAnalysisCLP.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class imageT >
void WriteBasis( const typename imageT::Pointer & img,
  std::string base, std::string ext, int num )
{
  typedef itk::ImageFileWriter< imageT >     BasisImageWriterType;

  typename BasisImageWriterType::Pointer basisImageWriter =
    BasisImageWriterType::New();
  std::string basename = base;
  char c[4096];
  std::sprintf( c, ext.c_str(), num );
  basename += std::string( c );
  basisImageWriter->SetUseCompression( true );
  basisImageWriter->SetFileName( basename.c_str() );
  basisImageWriter->SetInput( img );
  basisImageWriter->Update();
}

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // The timeCollector is used to perform basic profiling of the components
  //   of your algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef TPixel                                   InputPixelType;
  typedef itk::Image< InputPixelType, TDimension > InputImageType;
  typedef itk::Image< unsigned short, TDimension > MaskImageType;
  typedef itk::Image< float, TDimension >          BasisImageType;

  typedef itk::ImageFileReader< BasisImageType >   ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< BasisImageType >   BasisImageWriterType;

  typedef itk::tube::NJetSupervisedLinearBasisGenerator< BasisImageType,
          MaskImageType >            SupervisedLinearBasisGeneratorType;

  typename SupervisedLinearBasisGeneratorType::Pointer basisGenerator =
    SupervisedLinearBasisGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int vNum = 0; vNum < inputVolumesList.size(); vNum++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[vNum].c_str() );
    reader->Update();
    if( vNum == 0 )
      {
      basisGenerator->SetNJetImage( reader->GetOutput() );
      }
    else
      {
      basisGenerator->AddNJetImage( reader->GetOutput() );
      }
    }

  if( !labelmap.empty() )
    {
    typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
    inMaskReader->SetFileName( labelmap.c_str() );
    inMaskReader->Update();
    basisGenerator->SetLabelmap( inMaskReader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  if( !objectIdList.empty() )
    {
    basisGenerator->SetObjectId( objectIdList[0] );
    for( unsigned int o = 1; o < objectIdList.size(); o++ )
      {
      basisGenerator->AddObjectId( objectIdList[o] );
      }
    }

  if( usePCA )
    {
    basisGenerator->SetPerformPCA( true );
    }

  if( !loadBasisInfo.empty() )
    {
    timeCollector.Start( "LoadBasis" );

    itk::tube::MetaNJetLDA basisReader( loadBasisInfo.c_str() );
    basisReader.Read();

    basisGenerator->SetZeroScales( basisReader.GetZeroScales() );
    basisGenerator->SetFirstScales( basisReader.GetFirstScales() );
    basisGenerator->SetSecondScales( basisReader.GetSecondScales() );
    basisGenerator->SetRidgeScales( basisReader.GetRidgeScales() );
    basisGenerator->SetWhitenMeans( basisReader.GetWhitenMeans() );
    basisGenerator->SetWhitenStdDevs( basisReader.GetWhitenStdDevs() );
    basisGenerator->SetBasisValues( basisReader.GetLDAValues() );
    basisGenerator->SetBasisMatrix( basisReader.GetLDAMatrix() );

    basisGenerator->SetForceIntensityConsistency( !forceSignOff );
    basisGenerator->SetForceOrientationInsensitivity( forceSymmetry );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    basisGenerator->SetZeroScales( zeroScales );
    basisGenerator->SetFirstScales( firstScales );
    basisGenerator->SetSecondScales( secondScales );
    basisGenerator->SetRidgeScales( ridgeScales );
    basisGenerator->SetWhitenMeans( whitenMeans );
    basisGenerator->SetWhitenStdDevs( whitenStdDevs );

    basisGenerator->SetForceIntensityConsistency( !forceSignOff );
    basisGenerator->SetForceOrientationInsensitivity( forceSymmetry );

    basisGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  unsigned int numBasis = basisGenerator->GetNumberOfBasis();
  if( useNumberOfBasis > 0 &&
    useNumberOfBasis < static_cast<int>( numBasis ) )
    {
    numBasis = useNumberOfBasis;
    }

  if( !outputBase.empty() )
    {
    timeCollector.Start( "SaveBasisImages" );

    basisGenerator->UpdateBasisImages();

    for( unsigned int i = 0; i < numBasis; i++ )
      {
      typename BasisImageWriterType::Pointer basisImageWriter =
        BasisImageWriterType::New();
      std::string basename = outputBase;
      char c[4096];
      std::sprintf( c, ".basis%02u.mha", i );
      basename += std::string( c );
      basisImageWriter->SetUseCompression( true );
      basisImageWriter->SetFileName( basename.c_str() );
      basisImageWriter->SetInput( basisGenerator->GetBasisImage( i ) );
      basisImageWriter->Update();
      }
    timeCollector.Stop( "SaveBasisImages" );
    }

  if( !saveBasisInfo.empty() )
    {
    timeCollector.Start( "SaveBasis" );
    itk::tube::MetaNJetLDA basisWriter(
      basisGenerator->GetZeroScales(),
      basisGenerator->GetFirstScales(),
      basisGenerator->GetSecondScales(),
      basisGenerator->GetRidgeScales(),
      basisGenerator->GetBasisValues(),
      basisGenerator->GetBasisMatrix(),
      basisGenerator->GetWhitenMeans(),
      basisGenerator->GetWhitenStdDevs() );
    basisWriter.Write( saveBasisInfo.c_str() );
    timeCollector.Stop( "SaveBasis" );
    }

  if( !saveFeatureImages.empty() )
    {
    unsigned int numFeatures = basisGenerator->GetNumberOfFeatures();
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      WriteBasis< BasisImageType >( basisGenerator->GetFeatureImage( i ),
        saveFeatureImages, ".f%02d.mha", i );
      }
    }

  timeCollector.Report();

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolumesList[0], argc, argv );
}
