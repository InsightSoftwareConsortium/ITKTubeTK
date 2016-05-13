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

#include <ios>

#include "tubeCLIFilterWatcher.h"
#include "tubeCLIProgressReporter.h"
#include "tubeMessage.h"
#include "tubeStringUtilities.h"

#include <itkTimeProbesCollectorBase.h>
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>

#include <metaUtils.h>

// Must include CLP before including tubeCLIHelperFunctions
#include "ConvertCSVToImagesCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of
// int DoIt( ... ).
#include "tubeCLIHelperFunctions.h"

// Your code should be within the DoIt function...
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef float                                     InputPixelType;
  typedef itk::Image< InputPixelType, VDimension >  InputImageType;
  typedef itk::ImageFileReader< InputImageType >    ReaderType;

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( inputImageFileName );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
                        + std::string(err.GetDescription()) );
    return EXIT_FAILURE;
    }
  typename InputImageType::Pointer maskImage = reader->GetOutput();

  std::ifstream inCSVFile( inputCSVFileName.c_str() );
  std::string header;
  std::getline( inCSVFile, header );

  std::vector< std::string > imageFileNameList;
  tube::StringToVector< std::string >( header, imageFileNameList );

  unsigned int numImages = imageFileNameList.size();

  std::vector< typename InputImageType::Pointer > imageList;
  imageList.resize( numImages );
  for( unsigned int i = 0; i < numImages; ++i )
    {
    imageList[i] = InputImageType::New();
    imageList[i]->CopyInformation( maskImage );
    imageList[i]->SetRegions( maskImage->GetLargestPossibleRegion() );
    imageList[i]->Allocate();
    imageList[i]->FillBuffer( 0 );
    }

  typedef typename itk::ImageRegionIterator< InputImageType > ImageIterType;
  std::vector< ImageIterType * > imageIter;
  imageIter.resize( numImages );
  for( unsigned int i = 0; i < numImages; ++i )
    {
    imageIter[i] = new ImageIterType( imageList[i],
      maskImage->GetLargestPossibleRegion() );
    }

  ImageIterType maskIter( maskImage,
    maskImage->GetLargestPossibleRegion() );

  if( stride < 1 )
    {
    stride = 1;
    }

  while( !maskIter.IsAtEnd() )
    {
    if( maskIter.Get() != 0 )
      {
      std::string valueString;
      std::getline( inCSVFile, valueString );

      std::vector< float > valueList;
      tube::StringToVector< float >( valueString, valueList );

      for( unsigned int i=0; i<numImages; ++i )
        {
        imageIter[i]->Set( valueList[i] );
        }
      }
    for( int s=0; s<stride && !maskIter.IsAtEnd(); ++s )
      {
      for( unsigned int i=0; i<numImages; ++i )
        {
        ++(*imageIter[i]);
        }
      ++maskIter;
      }
    }

  typedef itk::ImageFileWriter< InputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  for( unsigned int i=0; i<imageIter.size(); ++i )
    {
    char outName[4096];
    sprintf( outName, "%s.%s.%03d.mha", outputImageBaseFileName.c_str(),
     imageFileNameList[i].c_str(), i );
    writer->SetFileName( outName );
    writer->SetInput( imageList[i] );
    writer->Update();
    }
  for( unsigned int i=0; i<imageIter.size(); ++i )
    {
    delete imageIter[i];
    }
  imageIter.clear();

  inCSVFile.close();

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  //   you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );
}
