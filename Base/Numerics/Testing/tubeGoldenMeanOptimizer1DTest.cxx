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

#include "tubeGoldenMeanOptimizer1D.h"

#include <vnl/vnl_math.h>

class MyOGMFunc : public tube::UserFunction< double, double >
{
private:
  double cVal;

public:
  MyOGMFunc( void )
    {
    cVal = 0;
    }

  const double & Value( const double & x )
    {
    cVal = std::sin( x );
    return cVal;
    }

}; // End class MyOGMFunc

class MyOGMFunc2 : public tube::UserFunction< double, double >
{
private:
  double cVal;

public:
  MyOGMFunc2( void )
    {
    cVal = 0;
    }

  const double & Value( const double & x )
    {
    cVal = std::cos( x/100 );
    std::cout << x << " : " << cVal << std::endl;
    return cVal;
    }

}; // End class MyOGMFunc2

int tubeGoldenMeanOptimizer1DTest( int argc, char * argv[] )
{
  if( argc < 9 )
    {
    std::cout << "Requires 9 arguments" << std::endl;
    return EXIT_FAILURE;
    }

  double factor = std::atof( argv[1] );
  int funcType = std::atoi( argv[2] );
  double xMin = std::atof( argv[3] );
  double xMax = std::atof( argv[4] );
  int minimizing = std::atoi( argv[5] );
  double x = std::atof( argv[6] );
  double idealXPiScaling = std::atof( argv[7] );
  double idealV = std::atof( argv[8] );

  MyOGMFunc * myFunc = new MyOGMFunc();
  tube::GoldenMeanOptimizer1D * opt = new tube::GoldenMeanOptimizer1D( myFunc );
  MyOGMFunc2 * myFunc2 = NULL;
  if( funcType == 2 )
    {
    myFunc2 = new MyOGMFunc2();
    opt->Use( myFunc2 );
    }

  double epsilon = 0.000001 * factor;

  int returnStatus = EXIT_SUCCESS;

  opt->SetXMin( xMin );
  if( opt->GetXMin() != xMin )
    {
    std::cout << "xMin should be " << xMin << " and not "
              << opt->GetXMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetXMax( xMax );
  if( opt->GetXMax() != xMax )
    {
    std::cout << "xMax should be " << xMax << " and not "
              << opt->GetXMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetXStep( 0.1*factor );
  if( opt->GetXStep() != 0.1*factor )
    {
    std::cout << "xStep should be " << 0.1*factor
              << " and not " << opt->GetXStep() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetTolerance( 0.000001*factor );
  if( opt->GetTolerance() != 0.000001*factor )
    {
    std::cout << "tolerance should be " << 0.000001*factor
              << " and not " << opt->GetTolerance() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetMaxIterations( 300 );
  if( opt->GetMaxIterations() != 300 )
    {
    std::cout << "maxIterations should be 300 and not "
      << opt->GetMaxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( minimizing )
    {
    opt->SetSearchForMin( true );
    if( !opt->GetSearchForMin() )
      {
      std::cout << "searchForMin should be true and not false." << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    }
  else
    {
    opt->SetSearchForMin( false );
    if( opt->GetSearchForMin() )
      {
      std::cout << "searchForMin should be false and not true." << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    }

  double xVal = 0;
  if( !opt->Extreme( &x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double idealX = vnl_math::pi * idealXPiScaling;
  if( vnl_math_abs( idealX - x ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  x=" << x
      << " diff=" << vnl_math_abs( idealX - x )
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( vnl_math_abs( idealV - xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  delete opt;
  if( funcType == 2 )
    {
    delete myFunc2;
    }
  delete myFunc;

  return returnStatus;
}
