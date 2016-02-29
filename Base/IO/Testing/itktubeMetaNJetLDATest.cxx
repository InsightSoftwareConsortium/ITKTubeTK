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

#include "itktubeMetaNJetLDA.h"

int itktubeMetaNJetLDATest( int argc, char * argv[] )
{
  if( argc != 2 )
    {
    std::cout << "Usage: testname <tempfilename>" << std::endl;
    return EXIT_FAILURE;
    }

  std::vector< double > scales0( 2 );
  scales0[0] = 1;
  scales0[1] = 2;
  std::vector< double > scales1( 1 );
  scales1[0] = 3;
  std::vector< double > scales2( 1 );
  scales2[0] = 4;
  std::vector< double > scalesR( 1 );
  scalesR[0] = 5;

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
      m( i, j ) = i * 10 + j;
      }
    }

  itk::tube::MetaNJetLDA mlda1;
  mlda1.SetLDAValues( v );
  mlda1.SetLDAMatrix( m );
  mlda1.SetInputWhitenMeans( wm );
  mlda1.SetInputWhitenStdDevs( ws );
  mlda1.SetOutputWhitenMeans( ws );
  mlda1.SetOutputWhitenStdDevs( wm );
  mlda1.SetZeroScales( scales0 );
  mlda1.SetFirstScales( scales1 );
  mlda1.SetSecondScales( scales2 );
  mlda1.SetRidgeScales( scalesR );
  if( mlda1.GetLDAValues() != v
      || mlda1.GetLDAMatrix() != m
      || mlda1.GetZeroScales() != scales0
      || mlda1.GetFirstScales() != scales1
      || mlda1.GetSecondScales() != scales2
      || mlda1.GetRidgeScales() != scalesR )
    {
    std::cout << "LDA values do not match after set." << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::MetaNJetLDA mlda2( mlda1 );
  if( mlda2.GetLDAValues() != mlda1.GetLDAValues()
      || mlda2.GetLDAMatrix() != mlda1.GetLDAMatrix()
      || mlda2.GetInputWhitenMeans() != mlda1.GetInputWhitenMeans()
      || mlda2.GetInputWhitenStdDevs() != mlda1.GetInputWhitenStdDevs()
      || mlda2.GetOutputWhitenMeans() != mlda1.GetOutputWhitenMeans()
      || mlda2.GetOutputWhitenStdDevs() != mlda1.GetOutputWhitenStdDevs()
      || mlda2.GetZeroScales() != mlda1.GetZeroScales()
      || mlda2.GetFirstScales() != mlda1.GetFirstScales()
      || mlda2.GetSecondScales() != mlda1.GetSecondScales()
      || mlda2.GetRidgeScales() != mlda1.GetRidgeScales() )
    {
    std::cout << "LDA values do not match after copy constructor."
      << std::endl;
    mlda1.PrintInfo();
    for( unsigned int i = 0; i < mlda2.GetZeroScales().size(); i++ )
      {
      std::cout << i << " s0 : " << mlda2.GetZeroScales()[i] << std::endl;
      }
    for( unsigned int i = 0; i < mlda2.GetFirstScales().size(); i++ )
      {
      std::cout << i << " s1 : " << mlda2.GetFirstScales()[i] << std::endl;
      }
    for( unsigned int i = 0; i < mlda2.GetSecondScales().size(); i++ )
      {
      std::cout << i << " s2 : " << mlda2.GetSecondScales()[i] << std::endl;
      }
    for( unsigned int i = 0; i < mlda2.GetRidgeScales().size(); i++ )
      {
      std::cout << i << " sR : " << mlda2.GetRidgeScales()[i] << std::endl;
      }
    return EXIT_FAILURE;
    }

  itk::tube::MetaNJetLDA mlda3( scales0, scales1, scales2, scalesR, 4, 6, v,
    m, wm, ws, ws, wm );
  if( mlda3.GetLDAValues() != mlda1.GetLDAValues()
      || mlda3.GetLDAMatrix() != mlda1.GetLDAMatrix()
      || mlda3.GetInputWhitenMeans() != mlda1.GetInputWhitenMeans()
      || mlda3.GetInputWhitenStdDevs() != mlda1.GetInputWhitenStdDevs()
      || mlda3.GetOutputWhitenMeans() != mlda1.GetOutputWhitenMeans()
      || mlda3.GetOutputWhitenStdDevs() != mlda1.GetOutputWhitenStdDevs()
      || mlda3.GetZeroScales() != scales0
      || mlda3.GetFirstScales() != scales1
      || mlda3.GetSecondScales() != scales2
      || mlda3.GetRidgeScales() != scalesR )
    {
    std::cout << "LDA values do not match after explicit constructor."
              << std::endl;
    return EXIT_FAILURE;
    }

  mlda1.Write( argv[1] );

  itk::tube::MetaNJetLDA mlda4( argv[1] );
  if( mlda4.GetLDAValues() != mlda1.GetLDAValues()
      || mlda4.GetLDAMatrix() != mlda1.GetLDAMatrix()
      || mlda4.GetInputWhitenMeans() != mlda1.GetInputWhitenMeans()
      || mlda4.GetInputWhitenStdDevs() != mlda1.GetInputWhitenStdDevs()
      || mlda4.GetOutputWhitenMeans() != mlda1.GetOutputWhitenMeans()
      || mlda4.GetOutputWhitenStdDevs() != mlda1.GetOutputWhitenStdDevs()
      || mlda4.GetZeroScales() != scales0
      || mlda4.GetFirstScales() != scales1
      || mlda4.GetSecondScales() != scales2
      || mlda4.GetRidgeScales() != scalesR )
    {
    std::cout << "LDA values do not match after write/read constructor."
              << std::endl;
    return EXIT_FAILURE;
    }

  itk::tube::MetaNJetLDA mlda5;
  mlda5.InitializeEssential( scales0, scales1, scales2, scalesR, 4, 6, v, m,
    wm, ws, ws, wm );
  if( mlda5.GetLDAValues() != mlda1.GetLDAValues()
      || mlda5.GetLDAMatrix() != mlda1.GetLDAMatrix()
      || mlda5.GetInputWhitenMeans() != mlda1.GetInputWhitenMeans()
      || mlda5.GetInputWhitenStdDevs() != mlda1.GetInputWhitenStdDevs()
      || mlda5.GetOutputWhitenMeans() != mlda1.GetOutputWhitenMeans()
      || mlda5.GetOutputWhitenStdDevs() != mlda1.GetOutputWhitenStdDevs()
      || mlda5.GetZeroScales() != scales0
      || mlda5.GetFirstScales() != scales1
      || mlda5.GetSecondScales() != scales2
      || mlda5.GetRidgeScales() != scalesR )
    {
    std::cout << "LDA values do not match after initialize essential."
              << std::endl;
    return EXIT_FAILURE;
    }

  mlda5.Clear();
  if( mlda5.GetLDAValues().size() != 0 )
    {
    std::cout << "LDA size not zero after clear." << std::endl;
    return EXIT_FAILURE;
    }

  mlda5.Read( argv[1] );
  if( mlda5.GetLDAValues() != mlda1.GetLDAValues()
      || mlda5.GetLDAMatrix() != mlda1.GetLDAMatrix()
      || mlda5.GetInputWhitenMeans() != mlda1.GetInputWhitenMeans()
      || mlda5.GetInputWhitenStdDevs() != mlda1.GetInputWhitenStdDevs()
      || mlda5.GetOutputWhitenMeans() != mlda1.GetOutputWhitenMeans()
      || mlda5.GetOutputWhitenStdDevs() != mlda1.GetOutputWhitenStdDevs()
      || mlda1.GetZeroScales() != scales0
      || mlda1.GetFirstScales() != scales1
      || mlda1.GetSecondScales() != scales2
      || mlda1.GetRidgeScales() != scalesR )
    {
    std::cout << "LDA values do not match after read." << std::endl;
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
