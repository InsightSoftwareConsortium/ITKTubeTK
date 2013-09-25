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

#include "itktubeMetaRidgeSeed.h"

int itktubeMetaRidgeSeedTest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cout << "Usage: testname <tempfilename>" << std::endl;
    return EXIT_FAILURE;
    }

  std::vector< double > scales( 2 );
  scales[0] = 1;
  scales[1] = 2;

  vnl_vector< double > v( 10 );
  vnl_matrix< double > m( 10, 10 );
  std::vector< double > wm( 10 );
  std::vector< double > ws( 10 );

  for( unsigned int i = 0; i < 10; i++ )
    {
    v[i] = i;
    wm[i] = i;
    ws[i] = i;
    for( unsigned int j = 0; j < 10; j++ )
      {
      m( i, j ) = i*10 + j;
      }
    }

  itk::tube::MetaRidgeSeed mrs1;
  mrs1.SetLDAValues( v );
  mrs1.SetLDAMatrix( m );
  mrs1.SetWhitenMeans( wm );
  mrs1.SetWhitenStdDevs( ws );
  mrs1.SetRidgeSeedScales( scales );
  mrs1.SetPDFFileName( "test.pdf" );
  if( mrs1.GetLDAValues() != v || mrs1.GetLDAMatrix() != m
    || mrs1.GetRidgeSeedScales() != scales )
    {
    std::cout << "LDA values do not match after set"
      << std::endl;
    return EXIT_FAILURE;
    }


  itk::tube::MetaRidgeSeed mrs2( mrs1 );
  if( mrs2.GetLDAValues() != mrs1.GetLDAValues()
    || mrs2.GetLDAMatrix() != mrs1.GetLDAMatrix()
    || mrs2.GetWhitenMeans() != mrs1.GetWhitenMeans()
    || mrs2.GetWhitenStdDevs() != mrs1.GetWhitenStdDevs()
    || mrs2.GetRidgeSeedScales() != mrs1.GetRidgeSeedScales()
    || mrs2.GetPDFFileName() != mrs1.GetPDFFileName() )
    {
    std::cout << "LDA values do not match after copy constructor"
      << std::endl;
    mrs1.PrintInfo();
    for( unsigned int i = 0; i < mrs2.GetRidgeSeedScales().size(); i++ )
      {
      std::cout << i << " sR : " << mrs2.GetRidgeSeedScales()[i]
        << std::endl;
      }
    std::cout << " PDFFileName : " << mrs2.GetPDFFileName() << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::MetaRidgeSeed mrs3( scales, v, m, wm, ws, "test.pdf" );
  if( mrs3.GetLDAValues() != mrs1.GetLDAValues()
    || mrs3.GetLDAMatrix() != mrs1.GetLDAMatrix()
    || mrs3.GetWhitenMeans() != mrs1.GetWhitenMeans()
    || mrs3.GetWhitenStdDevs() != mrs1.GetWhitenStdDevs()
    || mrs3.GetRidgeSeedScales() != scales
    || strcmp( mrs3.GetPDFFileName().c_str(), "test.pdf" ) )
    {
    std::cout << "LDA values do not match after explicit constructor"
      << std::endl;
    return EXIT_FAILURE;
    }

  mrs1.Write( argv[1] );

  itk::tube::MetaRidgeSeed mrs4( argv[1] );
  if( mrs4.GetLDAValues() != mrs1.GetLDAValues()
    || mrs4.GetLDAMatrix() != mrs1.GetLDAMatrix()
    || mrs4.GetWhitenMeans() != mrs1.GetWhitenMeans()
    || mrs4.GetWhitenStdDevs() != mrs1.GetWhitenStdDevs()
    || mrs4.GetRidgeSeedScales() != scales
    || mrs4.GetPDFFileName() != mrs1.GetPDFFileName() )
    {
    std::cout << "LDA values do not match after write/read constructor"
      << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::MetaRidgeSeed mrs5;
  mrs5.InitializeEssential( scales, v, m, wm, ws, "test.pdf" );
  if( mrs5.GetLDAValues() != mrs1.GetLDAValues()
    || mrs5.GetLDAMatrix() != mrs1.GetLDAMatrix()
    || mrs5.GetWhitenMeans() != mrs1.GetWhitenMeans()
    || mrs5.GetWhitenStdDevs() != mrs1.GetWhitenStdDevs()
    || mrs5.GetRidgeSeedScales() != scales
    || strcmp( mrs5.GetPDFFileName().c_str(), "test.pdf" ) )
    {
    std::cout << "LDA values don't match after InitializeEssential"
      << std::endl;
    return EXIT_FAILURE;
    }

  mrs5.Clear();
  if( mrs5.GetLDAValues().size() != 0 )
    {
    std::cout << "LDA size not 0 after clear." << std::endl;
    return EXIT_FAILURE;
    }

  mrs5.Read( argv[1] );
  if( mrs5.GetLDAValues() != mrs1.GetLDAValues()
    || mrs5.GetLDAMatrix() != mrs1.GetLDAMatrix()
    || mrs5.GetWhitenMeans() != mrs1.GetWhitenMeans()
    || mrs5.GetWhitenStdDevs() != mrs1.GetWhitenStdDevs()
    || mrs5.GetRidgeSeedScales() != scales
    || strcmp( mrs5.GetPDFFileName().c_str(), "test.pdf" ) )
    {
    std::cout << "LDA values do not match after read" << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
