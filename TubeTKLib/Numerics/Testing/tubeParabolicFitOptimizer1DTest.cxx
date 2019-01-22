/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#include "tubeParabolicFitOptimizer1D.h"

#include <vnl/vnl_math.h>

class MyOPFunc : public tube::UserFunction< double, double >
{
private:
  double cVal;

public:
  MyOPFunc( void )
    {
    cVal = 0;
    }
  const double & Value( const double & x )
    {
    cVal = std::sin( x );
    return cVal;
    }

}; // End class MyOPFunc

class MyOPFuncD : public tube::UserFunction< double, double >
{
private:
  double cDeriv;

public:
  MyOPFuncD( void )
    {
    cDeriv = 0;
    }
  const double & Value( const double & x )
    {
    cDeriv = std::cos( x );
    return cDeriv;
    }

}; // End class MyOPFuncD

int tubeParabolicFitOptimizer1DTest( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  double epsilon = 0.000001;

  MyOPFunc * myFunc = new MyOPFunc();
  tube::ParabolicFitOptimizer1D * opt = new tube::ParabolicFitOptimizer1D( myFunc );

  MyOPFunc * myFunc2 = new MyOPFunc();
  opt->Use( myFunc2 );

  delete myFunc;

  int returnStatus = EXIT_SUCCESS;

  opt->SetXMin( -3.5 );
  if( opt->GetXMin() != -3.5 )
    {
    std::cout << "xMin should be -1 and not " << opt->GetXMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetXMax( 3.5 );
  if( opt->GetXMax() != 3.5 )
    {
    std::cout << "xMax should be 1 and not " << opt->GetXMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetXStep( 0.1 );
  if( opt->GetXStep() != 0.1 )
    {
    std::cout << "xStep should be 0.1 and not " << opt->GetXStep()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetTolerance( 0.000001 );
  if( opt->GetTolerance() != 0.000001 )
    {
    std::cout << "tolerance should be 0.001 and not " << opt->GetTolerance()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetMaxIterations( 300 );
  if( opt->GetMaxIterations() != 300 )
    {
    std::cout << "maxIterations should be 100 and not "
      << opt->GetMaxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->SetSearchForMin( true );
  if( !opt->GetSearchForMin() )
    {
    std::cout << "searchForMin should be false and not true." << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double x = 0.01;
  double xVal = 0;
  if( !opt->Extreme( &x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double idealX = -vnl_math::pi / 2;
  if( std::fabs( idealX - x ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  x=" << x
      << " diff=" << std::fabs( idealX - x ) << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( std::fabs( -1 - xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  delete opt;
  delete myFunc2;

  return returnStatus;
}
