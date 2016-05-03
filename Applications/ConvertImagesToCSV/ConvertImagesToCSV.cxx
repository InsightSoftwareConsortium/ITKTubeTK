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
#include "ConvertImagesToCSVCLP.h"

// Must do a forward declaration of DoIt before including
// tubeCLIHelperFunctions
template< class TPixel, unsigned int VDimension >
int DoIt( int argc, char * argv[] );

// Must follow include of "...CLP.h" and forward declaration of int DoIt( ... ).
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

  unsigned int numImages = 0;
  std::vector< typename InputImageType::Pointer > imageList;
  std::vector< std::string > imageFileNameList;
  tube::StringToVector< std::string >( inputImageFileNameList,
    imageFileNameList );
  std::ofstream outFile( outputCSVFileName.c_str() );
  for( unsigned int i = 0; i < imageFileNameList.size(); ++i )
    {
    reader = ReaderType::New();
    reader->SetFileName( imageFileNameList[i] );
    char filePath[4096];
    std::string fileName = imageFileNameList[i];
    if( MET_GetFilePath( imageFileNameList[i].c_str(), filePath ) )
      {
      fileName = &( imageFileNameList[i][ strlen( filePath ) ] );
      }
    outFile << fileName << ", ";
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
    imageList.push_back( reader->GetOutput() );
    ++numImages;
    }
  outFile << "Class" << std::endl;

  typedef itk::ImageRegionIterator< InputImageType > IteratorType;

  std::vector< IteratorType * > iterList;
  for( unsigned int i = 0; i < numImages; ++i )
    {
    iterList.push_back( new IteratorType( imageList[i],
        imageList[i]->GetLargestPossibleRegion() ) );
    }

  if( stride < 1 )
    {
    stride = 1;
    }

  IteratorType maskIter( maskImage, maskImage->GetLargestPossibleRegion() );
  while( !maskIter.IsAtEnd() )
    {
    if( maskIter.Get() != 0 )
      {
      for( unsigned int i=0; i<numImages; ++i )
        {
        outFile << iterList[i]->Get() << ", ";
        }
      outFile << maskIter.Get() << std::endl;
      }
    for( int s=0; s<stride && !maskIter.IsAtEnd(); ++s )
      {
      for( unsigned int i=0; i<numImages; ++i )
        {
        ++(*iterList[i]);
        }
      ++maskIter;
      }
    }

  for( unsigned int i=0; i<iterList.size(); ++i )
    {
    delete iterList[i];
    }
  iterList.clear();

  outFile.close();

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
