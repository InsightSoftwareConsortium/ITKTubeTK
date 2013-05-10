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

#include <cstdlib>
#include <iostream>

#include <vnl/vnl_math.h>
#include <vcl_cmath.h>

#include "tubeMacro.h"

#include "tubeOptParabolicFit1D.h"
#include "tubeOptimizer1D.h"
#include "tubeUserFunc.h"

class MyOPFunc : public tube::UserFunc< double, double >
{
private:
  double cVal;

public:
  MyOPFunc( void )
    {
    cVal = 0;
    }
  const double & value( const double & x )
    {
    cVal = vcl_sin(x);
    return cVal;
    }

}; // End class MyOPFunc

class MyOPFuncD : public tube::UserFunc< double, double >
{
private:
  double cDeriv;

public:
  MyOPFuncD( void )
    {
    cDeriv = 0;
    }
  const double & value( const double & x )
    {
    cDeriv = vcl_cos(x);
    return cDeriv;
    }

}; // End class MyOPFuncD

int tubeOptParabolicFitTest( int tubeNotUsed(argc), char *tubeNotUsed(argv)[] )
{
  double epsilon = 0.000001;

  MyOPFunc * myFunc = new MyOPFunc();
  tube::OptParabolicFit1D * opt = new tube::OptParabolicFit1D( myFunc );

  MyOPFunc * myFunc2 = new MyOPFunc();
  opt->use( myFunc2 );

  delete myFunc;

  int returnStatus = EXIT_SUCCESS;

  opt->xMin( -3.5 );
  if( opt->xMin() != -3.5 )
    {
    std::cout << "xMin should be -1 and not " << opt->xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->xMax( 3.5 );
  if( opt->xMax() != 3.5 )
    {
    std::cout << "xMax should be 1 and not " << opt->xMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->xStep( 0.1 );
  if( opt->xStep() != 0.1 )
    {
    std::cout << "xStep should be 0.1 and not " << opt->xStep()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->tolerance( 0.000001 );
  if( opt->tolerance() != 0.000001 )
    {
    std::cout << "tolerance should be 0.001 and not " << opt->tolerance()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->maxIterations( 300 );
  if( opt->maxIterations() != 300 )
    {
    std::cout << "maxIterations should be 100 and not "
      << opt->maxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->searchForMin( true );
  if( !opt->searchForMin() )
    {
    std::cout << "searchForMin should be false and not true." << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double x = 0.01;
  double xVal = 0;
  if( !opt->extreme( &x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double idealX = -vnl_math::pi / 2;
  if( vnl_math_abs( idealX - x ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  x=" << x
      << " diff=" << vnl_math_abs( idealX - x ) << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( vnl_math_abs( -1 - xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  delete opt;
  delete myFunc2;

  return returnStatus;
}
