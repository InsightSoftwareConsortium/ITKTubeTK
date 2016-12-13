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

#include "itktubeMetaLDA.h"
#include "itktubeBasisFeatureVectorGenerator.h"

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "EnhanceUsingDiscriminantAnalysisCLP.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // The timeCollector is used to perform basic profiling of the components
  //   of this algorithm.
  itk::TimeProbesCollectorBase timeCollector;

  typedef TPixel                                   InputPixelType;
  typedef itk::Image< InputPixelType, VDimension > InputImageType;
  typedef itk::Image< unsigned short, VDimension > MaskImageType;
  typedef itk::Image< float, VDimension >          BasisImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< BasisImageType >   BasisImageWriterType;

  typedef itk::tube::FeatureVectorGenerator< InputImageType >
    FeatureVectorGeneratorType;

  typedef itk::tube::BasisFeatureVectorGenerator< InputImageType,
    MaskImageType >
      BasisFeatureVectorGeneratorType;

  typename FeatureVectorGeneratorType::Pointer fvGenerator =
    FeatureVectorGeneratorType::New();

  typename BasisFeatureVectorGeneratorType::Pointer basisGenerator =
    BasisFeatureVectorGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int i = 0; i < inputVolumesList.size(); i++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[i].c_str() );
    reader->Update();
    if( i == 0 )
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

  basisGenerator->SetInputFeatureVectorGenerator( fvGenerator );

  basisGenerator->SetObjectId( objectId[0] );
  if( objectId.size() > 1 )
    {
    for( unsigned int o = 1; o < objectId.size(); o++ )
      {
      basisGenerator->AddObjectId( objectId[o] );
      }
    }

  if( loadBasisInfo.size() > 0 )
    {
    timeCollector.Start( "LoadBasis" );

    itk::tube::MetaLDA basisReader( loadBasisInfo.c_str() );
    basisReader.Read();

    basisGenerator->SetNumberOfPCABasisToUseAsFeatures(
      basisReader.GetNumberOfPCABasisToUseAsFeatures() );
    basisGenerator->SetNumberOfLDABasisToUseAsFeatures(
      basisReader.GetNumberOfLDABasisToUseAsFeatures() );

    basisGenerator->SetBasisValues( basisReader.GetLDAValues() );
    basisGenerator->SetBasisMatrix( basisReader.GetLDAMatrix() );
    basisGenerator->SetInputWhitenMeans( basisReader.
      GetInputWhitenMeans() );
    basisGenerator->SetInputWhitenStdDevs( basisReader.
      GetInputWhitenStdDevs() );
    basisGenerator->SetOutputWhitenMeans( basisReader.
      GetOutputWhitenMeans() );
    basisGenerator->SetOutputWhitenStdDevs( basisReader.
      GetOutputWhitenStdDevs() );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    basisGenerator->SetNumberOfPCABasisToUseAsFeatures( useNumberOfPCABasis );
    if( useNumberOfLDABasis == -1 )
      {
      basisGenerator->SetNumberOfLDABasisToUseAsFeatures(
        objectId.size() - 1 );
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

  if( outputBase.size() > 0 )
    {
    timeCollector.Start( "SaveBasisImages" );

    unsigned int numBasis = basisGenerator->GetNumberOfFeatures();
    for( unsigned int i = 0; i < numBasis; i++ )
      {
      typename BasisImageWriterType::Pointer basisImageWriter =
        BasisImageWriterType::New();
      std::string fname = outputBase;
      char c[4096];
      std::sprintf( c, ".basis%02u.mha", i );
      fname += std::string( c );
      basisImageWriter->SetUseCompression( true );
      basisImageWriter->SetFileName( fname.c_str() );
      basisImageWriter->SetInput( basisGenerator->GetFeatureImage( i ) );
      basisImageWriter->Update();
      }
    timeCollector.Stop( "SaveBasisImages" );
    }

  if( saveBasisInfo.size() > 0 )
    {
    timeCollector.Start( "SaveBasis" );
    itk::tube::MetaLDA basisWriter(
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
