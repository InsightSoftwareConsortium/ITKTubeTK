/*=========================================================================

Library:   TubeTK

Copyright 2012 Kitware Inc. 28 Corporate Drive,
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

#include "itkLabelMapToAcousticImpedanceImageFilter.h"

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

template< class TLookupTable >
int ReadLookupTableFromCSV( const char * filename, TLookupTable & lookupTable );

int itkLabelMapToAcousticImpedanceImageFilterTest( int argc, char * argv[] )
{
  // Argument parsing.
  if( argc < 4 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " lookupTable.csv inputLabelMap.mha outputAcousticImpedance.mha" << std::endl;
    return EXIT_FAILURE;
    }
  const char * lookupTableFileName = argv[1];
  const char * labelMap = argv[2];
  const char * acousticImpedance = argv[3];


  // Types.
  enum { Dimension = 2 };

  typedef unsigned char                              LabelMapPixelType;
  typedef itk::Image< LabelMapPixelType, Dimension > LabelMapType;

  typedef float
    AcousticImpedancePixelType;
  typedef itk::Image< AcousticImpedancePixelType, Dimension >
    AcousticImpedanceImageType;

  typedef std::vector< float > LookupTableType;


  // Reader.
  typedef itk::ImageFileReader< LabelMapType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( labelMap );


  // Filter.
  typedef itk::LabelMapToAcousticImpedanceImageFilter< LabelMapType,
    AcousticImpedanceImageType, LookupTableType >
      LabelMapToAcousticImpedanceImageFilterType;

  typedef LabelMapToAcousticImpedanceImageFilterType::FunctorType FunctorType;
  FunctorType::LookupTableType lookupTable;

  if( ReadLookupTableFromCSV< LookupTableType >( lookupTableFileName, lookupTable )
    == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  FunctorType functor;
  functor.SetLookupTable( &lookupTable );

  LabelMapToAcousticImpedanceImageFilterType::Pointer filter =
    LabelMapToAcousticImpedanceImageFilterType::New();
  filter->SetFunctor( functor );
  filter->SetInput( reader->GetOutput() );


  // Writer.
  typedef itk::ImageFileWriter< AcousticImpedanceImageType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( filter->GetOutput() );
  writer->SetFileName( acousticImpedance );
  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & e )
    {
    std::cerr << "Error during pipeline Update: " << e << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}

template< class TLookupTable >
int ReadLookupTableFromCSV( const char * filename, TLookupTable & lookupTable )
{
  std::ifstream inputStream( filename );
  if( !inputStream.is_open() )
    {
    std::cerr << "Could not open input file: " << filename << std::endl;
    return EXIT_FAILURE;
    }

  // Ignore the header
  std::string str;
  std::getline( inputStream, str );
  if( !inputStream.good() )
    {
    inputStream.close();
    return EXIT_FAILURE;
    }

  unsigned int label = 0;
  char tissueType[256];
  float acousticImpedance;
  inputStream >> label;
  while( inputStream.good() )
    {
    lookupTable.resize( label + 1 );
    inputStream.get(); // ','
    inputStream.getline( tissueType, 256, ',' );
    inputStream >> acousticImpedance;
    std::getline( inputStream, str );
    lookupTable[label] = acousticImpedance;
    inputStream >> label;
    }

  inputStream.close();
  return EXIT_SUCCESS;
}
