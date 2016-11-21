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

#include "tubeMessage.h"
#include "tubeMacro.h"
#include "tubeStringUtilities.h"

// TubeTK includes
#include "tubeConvertImagesToCSV.h"

// ITK includes
#include "itkCSVNumericObjectFileWriter.h"
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

  typedef float                                     MaskPixelType;
  typedef itk::Image< MaskPixelType, VDimension >   InputMaskType;
  typedef itk::ImageFileReader< InputMaskType >     MaskReaderType;

  typedef itk::tube::ConvertImagesToCSVFilter< InputMaskType, InputImageType >
    ConvertImagesToCSVFilterType;
  typename ConvertImagesToCSVFilterType::Pointer filter
          = ConvertImagesToCSVFilterType::New();

  typename MaskReaderType::Pointer readerMask = MaskReaderType::New();

  readerMask->SetFileName( inputImageFileName );

  try
    {
    readerMask->Update();
    }
  catch( itk::ExceptionObject & err )
    {
    tube::ErrorMessage( "Reading volume: Exception caught: "
      + std::string(err.GetDescription()));
    return EXIT_FAILURE;
    }

  typename InputMaskType::Pointer inputMask = readerMask->GetOutput();

  unsigned int numImages = 0;
  std::vector< std::string > imageFileNameList;
  tube::StringToVector< std::string >( inputImageFileNameList,
    imageFileNameList );

  if( stride < 1 )
    {
    stride = 1;
    }
  typename ReaderType::Pointer reader;
  std::vector<std::string> fileName;
  for( unsigned int i = 0; i < imageFileNameList.size(); ++i )
    {
    reader = ReaderType::New();
    reader->SetFileName( imageFileNameList[i] );
    char filePath[4096];
    fileName.push_back( imageFileNameList[i] );
    if( MET_GetFilePath( imageFileNameList[i].c_str(), filePath ) )
      {
      fileName[i] = &( imageFileNameList[i][strlen(filePath)] );
      }
    try
      {
      reader->Update();
      }
    catch ( itk::ExceptionObject & err )
      {
      tube::ErrorMessage( "Reading volume: Exception caught: "
        + std::string(err.GetDescription()) );
      return EXIT_FAILURE;
      }
    filter->AddImage( reader->GetOutput() );
    ++numImages;
    }

  typedef vnl_matrix<InputPixelType> MatrixType;
  const unsigned int ARows =
    inputMask->GetLargestPossibleRegion().GetNumberOfPixels() / stride;
  const unsigned int ACols = imageFileNameList.size() + 1;
  MatrixType matrix;
  matrix.set_size( ARows, ACols );

  filter->SetInputMask( inputMask );
  filter->SetStride( stride );
  filter->SetNumImages( numImages );

  filter->Update();

  matrix = filter->GetOutput()->Get();
  unsigned int numberRows = filter->GetNumberRows();
  MatrixType submatrix = matrix.extract( numberRows, ACols );

  // write out the vnl_matrix object
  typedef itk::CSVNumericObjectFileWriter<InputPixelType> WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetFieldDelimiterCharacter( ',' );
  writer->SetFileName( outputCSVFileName );
  writer->SetInput( &submatrix );

  fileName.push_back( "Class" );
  writer->SetColumnHeaders( fileName );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject& exp )
    {
    std::cerr << "Exception caught!" << std::endl;
    std::cerr << exp << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

// Main
int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  // You may need to update this line if, in the project's .xml CLI file,
  // you change the variable name for the inputImageFileName.
  return tube::ParseArgsAndCallDoIt( inputImageFileName, argc, argv );

}
