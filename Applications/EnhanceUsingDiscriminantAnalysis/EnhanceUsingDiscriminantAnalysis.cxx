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

#include "itktubeMetaLDA.h"
#include "itktubeSupervisedLinearBasisGenerator.h"
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

  typedef itk::tube::SupervisedLinearBasisGenerator< InputImageType,
    MaskImageType >                  SupervisedLinearBasisGeneratorType;

  typename SupervisedLinearBasisGeneratorType::Pointer basisGenerator =
    SupervisedLinearBasisGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int i = 0; i < inputVolumesList.size(); i++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[i].c_str() );
    reader->Update();
    if( i == 0 )
      {
      basisGenerator->SetFeatureImage( reader->GetOutput() );
      }
    else
      {
      basisGenerator->AddFeatureImage( reader->GetOutput() );
      }
    }

  if( labelmap.size() > 0 )
    {
    typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
    inMaskReader->SetFileName( labelmap.c_str() );
    inMaskReader->Update();
    basisGenerator->SetLabelMap( inMaskReader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  basisGenerator->SetObjectId( objectId[0] );
  if( objectId.size() > 1 )
    {
    for( unsigned int o = 1; o < objectId.size(); o++ )
      {
      basisGenerator->AddObjectId( objectId[o] );
      }
    }

  if( usePCA )
    {
    basisGenerator->SetPerformPCA( true );
    }

  if( loadBasisInfo.size() > 0 )
    {
    timeCollector.Start( "LoadBasis" );

    itk::tube::MetaLDA basisReader( loadBasisInfo.c_str() );
    basisReader.Read();

    basisGenerator->SetBasisValues( basisReader.GetLDAValues() );
    basisGenerator->SetBasisMatrix( basisReader.GetLDAMatrix() );
    basisGenerator->SetWhitenMeans( basisReader.GetWhitenMeans() );
    basisGenerator->SetWhitenStdDevs( basisReader.GetWhitenStdDevs() );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    basisGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  if( outputBase.size() > 0 )
    {
    timeCollector.Start( "SaveBasisImages" );

    unsigned int numBasis = basisGenerator->GetNumberOfBasis();
    if( useNumberOfBasis>0 && useNumberOfBasis < (int)numBasis )
      {
      numBasis = useNumberOfBasis;
      }

    basisGenerator->UpdateBasisImages();

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
      basisImageWriter->SetInput( basisGenerator->GetBasisImage( i ) );
      basisImageWriter->Update();
      }
    timeCollector.Stop( "SaveBasisImages" );
    }

  if( saveBasisInfo.size() > 0 )
    {
    timeCollector.Start( "SaveBasis" );
    itk::tube::MetaLDA basisWriter( basisGenerator->GetBasisValues(),
      basisGenerator->GetBasisMatrix(),
      basisGenerator->GetWhitenMeans(),
      basisGenerator->GetWhitenStdDevs() );
    basisWriter.Write( saveBasisInfo.c_str() );
    timeCollector.Stop( "SaveBasis" );
    }

  timeCollector.Report();

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char *[] argv )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputVolume.
  return tube::ParseArgsAndCallDoIt( inputVolumesList[0], argc, argv );
}
