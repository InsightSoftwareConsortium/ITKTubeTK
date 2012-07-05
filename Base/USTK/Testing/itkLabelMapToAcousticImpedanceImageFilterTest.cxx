/*=========================================================================

Library:   TubeTK

Copyright 2012 Kitware Inc. 28 Corporate Drive,
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
#include <itkLabelMapToAcousticImpedanceImageFilter.h>
#include <fstream>


template< class TLookupTable >
int ReadLookupTableFromCSV( const char * filename, TLookupTable & lookupTable );

int itkLabelMapToAcousticImpedanceImageFilterTest( int argc, char * argv [] )
{
  if( argc < 2 )
    {
    std::cerr << "Missing arguments." << std::endl;
    std::cerr << "Usage: " << std::endl;
    std::cerr << argv[0] << " lookupTable.csv " << std::endl;
    return EXIT_FAILURE;
    }
  const char * lookupTableFileName = argv[1];

  static const unsigned int Dimension = 2;

  typedef unsigned char                              LabelMapPixelType;
  typedef itk::Image< LabelMapPixelType, Dimension > LabelMapType;

  typedef float
    AcousticImpedancePixelType;
  typedef itk::Image< AcousticImpedancePixelType, Dimension >
    AcousticImpedanceImageType;

  typedef std::vector< float > LookupTableType;

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

  size_t label = 0;
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
