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

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTimeProbesCollectorBase.h>

#include "itktubeMetaRidgeSeed.h"
#include "itktubeRidgeSeedFilter.h"
#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include "EnhanceTubesUsingDiscriminantAnalysisCLP.h"


template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h"
//   and forward declaration of int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

template< class TPixel, unsigned int TDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  std::vector< std::string > inputVolumesList;
  tube::StringToVector< std::string >( inputVolumesString,
    inputVolumesList );

  itk::TimeProbesCollectorBase timeCollector;

  typedef TPixel                                    InputPixelType;
  typedef itk::Image< InputPixelType, TDimension >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >    InputImageReaderType;

  typedef unsigned short                            LabelmapPixelType;
  typedef itk::Image< unsigned short, TDimension >  LabelmapType;
  typedef itk::ImageFileReader< LabelmapType >      LabelmapReaderType;

  typedef itk::tube::MetaRidgeSeed                  TubeParametersIOType;

  typedef itk::tube::RidgeSeedFilter< InputImageType, LabelmapType >
                                                    RidgeSeedFilterType;

  typedef typename RidgeSeedFilterType::ProbabilityPixelType
                                                    OutputPixelType;
  typedef typename RidgeSeedFilterType::ProbabilityImageType
                                                    OutputImageType;
  typedef itk::ImageFileWriter< OutputImageType >   OutputImageWriterType;


  typename RidgeSeedFilterType::Pointer tubeFilter =
    RidgeSeedFilterType::New();

  timeCollector.Start( "LoadData" );

  typename InputImageReaderType::Pointer reader;
  for( unsigned int vNum = 0; vNum < inputVolumesList.size(); vNum++ )
    {
    reader = InputImageReaderType::New();
    reader->SetFileName( inputVolumesList[vNum].c_str() );
    reader->Update();
    if( vNum == 0 )
      {
      tubeFilter->SetInput( reader->GetOutput() );
      }
    else
      {
      tubeFilter->AddInput( reader->GetOutput() );
      }
    }

  if( !labelmap.empty() )
    {
    typename LabelmapReaderType::Pointer  inLabelmapReader =
      LabelmapReaderType::New();
    inLabelmapReader->SetFileName( labelmap.c_str() );
    inLabelmapReader->Update();
    tubeFilter->SetLabelMap( inLabelmapReader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  if( !objectIdList.empty() )
    {
    tubeFilter->SetObjectId( objectIdList[0] );
    for( unsigned int o = 1; o < objectIdList.size(); o++ )
      {
      tubeFilter->AddObjectId( objectIdList[o] );
      }
    }

  if( !loadDiscriminantInfo.empty() )
    {
    timeCollector.Start( "LoadBasis" );

    TubeParametersIOType tubeParametersIO( loadDiscriminantInfo.c_str() );
    tubeParametersIO.Read();

    tubeFilter->SetScales( tubeParametersIO.GetRidgeSeedScales() );
    tubeFilter->SetWhitenMeans( tubeParametersIO.GetWhitenMeans() );
    tubeFilter->SetWhitenStdDevs( tubeParametersIO.GetWhitenStdDevs() );

    tubeFilter->SetBasisValues( tubeParametersIO.GetLDAValues() );
    tubeFilter->SetBasisMatrix( tubeParametersIO.GetLDAMatrix() );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    tubeFilter->SetScales( tubeScales );
    tubeFilter->SetWhitenMeans( whitenMeans );
    tubeFilter->SetWhitenStdDevs( whitenStdDevs );

    tubeFilter->Update();

    timeCollector.Stop( "Update" );
    }

  tubeFilter->ClassifyImages();

  if( !saveDiscriminantInfo.empty() )
    {
    timeCollector.Start( "SaveBasis" );
    TubeParametersIOType tubeParametersIO(
      tubeFilter->GetScales(),
      tubeFilter->GetBasisValues(),
      tubeFilter->GetBasisMatrix(),
      tubeFilter->GetWhitenMeans(),
      tubeFilter->GetWhitenStdDevs() );
    tubeParametersIO.Write( saveDiscriminantInfo.c_str() );
    timeCollector.Stop( "SaveBasis" );
    }

  typename OutputImageWriterType::Pointer outputWriter =
    OutputImageWriterType::New();
  outputWriter->SetFileName( outputVolume.c_str() );
  outputWriter->SetUseCompression( true );
  outputWriter->SetInput( tubeFilter->
    GetClassProbabilityForInput( 0 ) );
  outputWriter->Update();

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
