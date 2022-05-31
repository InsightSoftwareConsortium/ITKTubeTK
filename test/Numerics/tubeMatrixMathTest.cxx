/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "tubeMatrixMath.h"

#include <itkMersenneTwisterRandomVariateGenerator.h>

template< int VDimension >
int Test( void )
{
  double epsilon = 0.00001;

  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer rndGen
    = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  rndGen->Initialize( 1 );

  int returnStatus = EXIT_SUCCESS;

  for( unsigned int count=0; count<1000; count++ )
    {
    vnl_vector<float> v1( VDimension );
    for( unsigned int d=0; d<VDimension; d++ )
      {
      v1[d] = rndGen->GetNormalVariate( 0.0, 1.0 );
      }
    vnl_vector<float> v2( VDimension );
    if( VDimension == 3 )
      {
      v2 = tube::ComputeOrthogonalVector( v1 );
      if( std::fabs( dot_product( v1, v2 ) ) > epsilon )
        {
        std::cout << count << " : ";
        std::cout << "FAILURE: ComputeOrthogonalVector: DotProduct = "
                  << v1 << " .* " << v2 << " = "
                  << dot_product( v1, v2 ) << std::endl;
        returnStatus = EXIT_FAILURE;
        }

      vnl_vector<float> v3 = tube::ComputeCrossVector( v1, v2 );
      if( std::fabs( dot_product( v1, v3 ) ) > epsilon ||
          std::fabs( dot_product( v2, v3 ) ) > epsilon )
        {
        std::cout << count << " : ";
        std::cout << "FAILURE: ComputeCrossVector: DotProduct = "
          << dot_product( v1, v3 ) << " and "
          << dot_product( v2, v3 ) << std::endl;
        returnStatus = EXIT_FAILURE;
        }
      }
    else
      {
      for( unsigned int d=0; d<VDimension; d++ )
        {
        v2[d] = rndGen->GetNormalVariate( 0.0, 1.0 );
        }
      }

    v2 = v2.normalize();
    vnl_vector<float> v4 = tube::ComputeLineStep( v1, 0.5, v2 );
    if( std::fabs( tube::ComputeEuclideanDistanceVector( v1, v4 ) - 0.5 )
        > epsilon )
      {
      std::cout << count << " : ";
      std::cout << "FAILURE: ComputeLineStep = "
        << v1 << " + 0.5 * " << v2 << " != " << v4 << std::endl;
      std::cout << "FAILURE: ComputeEuclidenDistanceVector = "
        << tube::ComputeEuclideanDistanceVector( v1, v4 ) << std::endl;
      returnStatus = EXIT_FAILURE;
      }

    vnl_matrix<float> m1( VDimension, VDimension );
    for( unsigned int r=0; r<VDimension; r++ )
      {
      for( unsigned int c=r; c<VDimension; c++ )
        {
        m1( r,c ) = rndGen->GetNormalVariate( 0.0, 1.0 );
        m1( c,r ) = m1( r,c );
        }
      }

    vnl_matrix<float> eVects( VDimension, VDimension );
    vnl_vector<float> eVals( VDimension );
    tube::ComputeEigen( m1, eVects, eVals, true );
    for( unsigned int d=0; d<VDimension; d++ )
      {
      v1 = m1 * eVects.get_column( d );
      if( std::fabs( v1.magnitude() - std::fabs( eVals[d] ) ) > epsilon )
        {
        std::cout << count << " : ";
        std::cout << "FAILURE: ComputeEigen : "
          << " v1 * M1 = " << v1
          << " and v1 norm = " << v1.magnitude()
          << " != " << eVals[d]
          << std::endl;
        returnStatus = EXIT_FAILURE;
        }
      }
    }

  return returnStatus;
}

int tubeMatrixMathTest( int itkNotUsed( argc ), char * itkNotUsed( argv )[] )
{
  if( Test<2>() == EXIT_FAILURE ||
      Test<3>() == EXIT_FAILURE ||
      Test<4>() == EXIT_FAILURE )
    {
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
