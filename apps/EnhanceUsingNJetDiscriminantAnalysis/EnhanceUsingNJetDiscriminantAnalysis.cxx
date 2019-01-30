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

#include "itktubeMetaNJetLDA.h"
#include "itktubeNJetFeatureVectorGenerator.h"
#include "itktubeBasisFeatureVectorGenerator.h"

#include "../CLI/tubeCLIFilterWatcher.h"
#include "../CLI/tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "EnhanceUsingNJetDiscriminantAnalysisCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "../CLI/tubeCLIHelperFunctions.h"

template< class TImage >
void WriteBasis( const typename TImage::Pointer & img,
  std::string base, std::string ext, int num )
{
  typedef itk::ImageFileWriter< TImage >     BasisImageWriterType;

  typename BasisImageWriterType::Pointer basisImageWriter =
    BasisImageWriterType::New();
  std::string basename = base;
  char c[4096];
  std::snprintf( c, 4095, ext.c_str(), num );
  basename += std::string( c );
  basisImageWriter->SetUseCompression( true );
  basisImageWriter->SetFileName( basename.c_str() );
  basisImageWriter->SetInput( img );
  basisImageWriter->Update();
}

template< class TPixel, unsigned int VDimension >
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
  typedef itk::Image< InputPixelType, VDimension > InputImageType;
  typedef itk::Image< unsigned short, VDimension > MaskImageType;
  typedef itk::Image< float, VDimension >          BasisImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< BasisImageType >   BasisImageWriterType;

  typedef itk::tube::NJetFeatureVectorGenerator< InputImageType >
    NJetFeatureVectorGeneratorType;

  typedef itk::tube::FeatureVectorGenerator< InputImageType >
    FeatureVectorGeneratorType;

  typedef itk::tube::BasisFeatureVectorGenerator< InputImageType,
    MaskImageType >
      BasisFeatureVectorGeneratorType;

  typename NJetFeatureVectorGeneratorType::Pointer fvGenerator =
    NJetFeatureVectorGeneratorType::New();

  typename BasisFeatureVectorGeneratorType::Pointer basisGenerator =
    BasisFeatureVectorGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int vNum = 0; vNum < inputVolumesList.size(); vNum++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[vNum].c_str() );
    reader->Update();
    if( vNum == 0 )
      {
      fvGenerator->SetInput( reader->GetOutput() );
      basisGenerator->SetInput( reader->GetOutput() );
      }
    else
      {
      fvGenerator->AddInput( reader->GetOutput() );
      basisGenerator->AddInput( reader->GetOutput() );
      }
    }

  if( !labelmap.empty() )
    {
    typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
    inMaskReader->SetFileName( labelmap.c_str() );
    inMaskReader->Update();
    basisGenerator->SetLabelMap( inMaskReader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  basisGenerator->SetInputFeatureVectorGenerator( static_cast<
   FeatureVectorGeneratorType * >( fvGenerator.GetPointer() ) );

  if( !objectIdList.empty() )
    {
    basisGenerator->SetObjectId( objectIdList[0] );
    for( unsigned int o = 1; o < objectIdList.size(); o++ )
      {
      basisGenerator->AddObjectId( objectIdList[o] );
      }
    }

  if( !loadBasisInfo.empty() )
    {
    timeCollector.Start( "LoadBasis" );

    itk::tube::MetaNJetLDA basisReader( loadBasisInfo.c_str() );
    basisReader.Read();

    fvGenerator->SetZeroScales( basisReader.GetZeroScales() );
    fvGenerator->SetFirstScales( basisReader.GetFirstScales() );
    fvGenerator->SetSecondScales( basisReader.GetSecondScales() );
    fvGenerator->SetRidgeScales( basisReader.GetRidgeScales() );
    basisGenerator->SetNumberOfPCABasisToUseAsFeatures(
      basisReader.GetNumberOfPCABasisToUseAsFeatures() );
    basisGenerator->SetNumberOfLDABasisToUseAsFeatures(
      basisReader.GetNumberOfLDABasisToUseAsFeatures() );
    basisGenerator->SetInputWhitenMeans( basisReader.
      GetInputWhitenMeans() );
    basisGenerator->SetInputWhitenStdDevs( basisReader.
      GetInputWhitenStdDevs() );
    basisGenerator->SetOutputWhitenMeans( basisReader.
      GetOutputWhitenMeans() );
    basisGenerator->SetOutputWhitenStdDevs( basisReader.
      GetOutputWhitenStdDevs() );
    basisGenerator->SetBasisValues( basisReader.GetLDAValues() );
    basisGenerator->SetBasisMatrix( basisReader.GetLDAMatrix() );

    fvGenerator->SetUpdateWhitenStatisticsOnUpdate( false );
    basisGenerator->SetUpdateWhitenStatisticsOnUpdate( false );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    fvGenerator->SetZeroScales( zeroScales );
    fvGenerator->SetFirstScales( firstScales );
    fvGenerator->SetSecondScales( secondScales );
    fvGenerator->SetRidgeScales( ridgeScales );
    fvGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    fvGenerator->Update();

    basisGenerator->SetNumberOfPCABasisToUseAsFeatures(
      useNumberOfPCABasis );
    if( useNumberOfLDABasis == -1 )
      {
      basisGenerator->SetNumberOfLDABasisToUseAsFeatures(
        objectIdList.size() - 1 );
      }
    else
      {
      basisGenerator->SetNumberOfLDABasisToUseAsFeatures(
        useNumberOfLDABasis );
      }

    basisGenerator->SetUpdateWhitenStatisticsOnUpdate( true );
    basisGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  if( !outputBase.empty() )
    {
    timeCollector.Start( "SaveBasisImages" );

    unsigned int numBasis = basisGenerator->GetNumberOfFeatures();
    for( unsigned int i = 0; i < numBasis; i++ )
      {
      typename BasisImageWriterType::Pointer basisImageWriter =
        BasisImageWriterType::New();
      std::string basename = outputBase;
      char c[4096];
      std::snprintf( c, 4095, ".basis%02u.mha", i );
      basename += std::string( c );
      basisImageWriter->SetUseCompression( true );
      basisImageWriter->SetFileName( basename.c_str() );
      basisImageWriter->SetInput( basisGenerator->GetFeatureImage( i ) );
      basisImageWriter->Update();
      }
    timeCollector.Stop( "SaveBasisImages" );
    }

  if( !saveBasisInfo.empty() )
    {
    timeCollector.Start( "SaveBasis" );
    itk::tube::MetaNJetLDA basisWriter(
      fvGenerator->GetZeroScales(),
      fvGenerator->GetFirstScales(),
      fvGenerator->GetSecondScales(),
      fvGenerator->GetRidgeScales(),
      basisGenerator->GetNumberOfPCABasisToUseAsFeatures(),
      basisGenerator->GetNumberOfLDABasisToUseAsFeatures(),
      basisGenerator->GetBasisValues(),
      basisGenerator->GetBasisMatrix(),
      basisGenerator->GetInputWhitenMeans(),
      basisGenerator->GetInputWhitenStdDevs(),
      basisGenerator->GetOutputWhitenMeans(),
      basisGenerator->GetOutputWhitenStdDevs() );
    basisWriter.Write( saveBasisInfo.c_str() );
    timeCollector.Stop( "SaveBasis" );
    }

  if( !saveFeatureImages.empty() )
    {
    unsigned int numFeatures = fvGenerator->GetNumberOfFeatures();
    for( unsigned int i = 0; i < numFeatures; i++ )
      {
      WriteBasis< BasisImageType >( fvGenerator->GetFeatureImage( i ),
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
