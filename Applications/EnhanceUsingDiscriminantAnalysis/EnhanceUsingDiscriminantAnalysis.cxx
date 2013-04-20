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
#include "itkTubeSupervisedLinearBasisGenerator.h"
#include "itkTubeMetaLDA.h"

// Must do a forward declaraction of DoIt before including
// tubeCLIHelperFunctions
template< class pixelT, unsigned int dimensionT >
int DoIt( int argc, char * argv[] );

// Must include CLP before including tubeCLIHelperFunctions
#include "EnhanceUsingDiscriminantAnalysisCLP.h"

// Includes tube::ParseArgsAndCallDoIt function
#include "tubeCLIHelperFunctions.h"

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
  typedef itk::Image< InputPixelType, dimensionT > InputImageType;
  typedef itk::Image< unsigned short, dimensionT > MaskImageType;
  typedef itk::Image< float, dimensionT >          BasisImageType;

  typedef itk::ImageFileReader< InputImageType >   ImageReaderType;
  typedef itk::ImageFileReader< MaskImageType >    MaskReaderType;
  typedef itk::ImageFileWriter< BasisImageType >     BasisImageWriterType;

  typedef itk::tube::SupervisedLinearBasisGenerator< InputImageType, MaskImageType >
    SupervisedLinearBasisGeneratorType;
  typename SupervisedLinearBasisGeneratorType::Pointer ldaGenerator = SupervisedLinearBasisGeneratorType::New();

  timeCollector.Start( "LoadData" );

  typename ImageReaderType::Pointer reader;
  for( unsigned int i=0; i<inputVolumesList.size(); i++ )
    {
    reader = ImageReaderType::New();
    reader->SetFileName( inputVolumesList[i].c_str() );
    reader->Update();
    if( i == 0 )
      {
      ldaGenerator->SetFeatureImage( reader->GetOutput() );
      }
    else
      {
      ldaGenerator->AddFeatureImage( reader->GetOutput() );
      }
    }

  if( labelmap.size() > 0 )
    {
    typename MaskReaderType::Pointer  inMaskReader = MaskReaderType::New();
    inMaskReader->SetFileName( labelmap.c_str() );
    inMaskReader->Update();
    ldaGenerator->SetLabelmap( inMaskReader->GetOutput() );
    }

  timeCollector.Stop( "LoadData" );

  ldaGenerator->SetObjectId( objectId[0] );
  if( objectId.size() > 1 )
    {
    for( unsigned int o=1; o<objectId.size(); o++ )
      {
      ldaGenerator->AddObjectId( objectId[o] );
      }
    }

  if( usePCA )
    {
    ldaGenerator->SetPerformPCA( true );
    }

  if( loadBasisInfo.size() > 0 )
    {
    timeCollector.Start( "LoadBasis" );

    itk::tube::MetaLDA ldaReader( loadBasisInfo.c_str() );
    ldaReader.Read();

    ldaGenerator->SetBasisValues( ldaReader.GetLDAValues() );
    ldaGenerator->SetBasisMatrix( ldaReader.GetLDAMatrix() );
    ldaGenerator->SetWhitenMeans( ldaReader.GetWhitenMeans() );
    ldaGenerator->SetWhitenStdDevs( ldaReader.GetWhitenStdDevs() );

    timeCollector.Stop( "LoadBasis" );
    }
  else
    {
    timeCollector.Start( "Update" );

    ldaGenerator->Update();

    timeCollector.Stop( "Update" );
    }

  if( outputBase.size() > 0 )
    {
    timeCollector.Start( "SaveBasisImages" );

    unsigned int numBasis = ldaGenerator->GetNumberOfBasis();
    if( useNumberOfBasis>0 && useNumberOfBasis < (int)numBasis )
      {
      numBasis = useNumberOfBasis;
      }
    std::cout << "number of lda = " << numBasis << std::endl;

    ldaGenerator->UpdateBasisImages();

    for( unsigned int i=0; i<numBasis; i++ )
      {
      typename BasisImageWriterType::Pointer ldaImageWriter =
        BasisImageWriterType::New();
      std::string fname = outputBase;
      char c[80];
      sprintf( c, ".lda%02d.mha", i );
      fname += std::string( c );
      ldaImageWriter->SetUseCompression( true );
      ldaImageWriter->SetFileName( fname.c_str() );
      ldaImageWriter->SetInput( ldaGenerator->GetBasisImage( i ) );
      ldaImageWriter->Update();
      }
    timeCollector.Stop( "SaveBasisImages" );
    }

  if( saveBasisInfo.size() > 0 )
    {
    timeCollector.Start( "SaveBasis" );
    itk::tube::MetaLDA ldaWriter( ldaGenerator->GetBasisValues(),
      ldaGenerator->GetBasisMatrix(),
      ldaGenerator->GetWhitenMeans(),
      ldaGenerator->GetWhitenStdDevs() );
    ldaWriter.Write( saveBasisInfo.c_str() );
    timeCollector.Stop( "SaveBasis" );
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
