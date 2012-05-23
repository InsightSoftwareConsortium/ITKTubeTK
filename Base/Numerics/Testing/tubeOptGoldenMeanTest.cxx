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

#include "tubeOptGoldenMean1D.h"
#include "tubeOptimizer1D.h"
#include "tubeUserFunc.h"

class MyOGMFunc:
  public tube::UserFunc< double, double >
  {
  private:
    double cVal;
  public:
    MyOGMFunc( )
      {
      cVal = 0;
      };
    const double & value( const double & x )
      {
      cVal = vcl_sin(x);
      std::cout << x << " : " << cVal << std::endl;
      return cVal;
      };
  };

class MyOGMFunc2:
  public tube::UserFunc< double, double >
  {
  private:
    double cVal;
  public:
    MyOGMFunc2( )
      {
      cVal = 0;
      };
    const double & value( const double & x )
      {
      cVal = vcl_cos( x/100 );
      std::cout << x << " : " << cVal << std::endl;
      return cVal;
      };
  };

int tubeOptGoldenMeanTest( int argc, char *argv[] )
{
  int testNum = 0;
  if( argc > 1 )
    {
    testNum = atoi( argv[ 1 ] );
    }

  MyOGMFunc * myFunc = new MyOGMFunc();
  tube::OptGoldenMean1D * opt = new tube::OptGoldenMean1D( myFunc );

  double factor = 1;
  MyOGMFunc2 * myFunc2 = NULL;
  if( testNum != 0 )
    {
    factor = 100;
    myFunc2 = new MyOGMFunc2();
    opt->use( myFunc2 );

    delete myFunc;
    }

  double epsilon = 0.000001 * factor;

  int returnStatus = EXIT_SUCCESS;

  opt->xMin( -3*factor );
  if( opt->xMin() != -3*factor )
    {
    std::cout << "xMin should be -3.5 and not " << opt->xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->xMax( 3.5*factor );
  if( opt->xMax() != 3.5*factor )
    {
    std::cout << "xMax should be 3.5 and not " << opt->xMax() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->xStep( 0.1*factor );
  if( opt->xStep() != 0.1*factor )
    {
    std::cout << "xStep should be 0.1 and not " << opt->xStep()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->tolerance( 0.000001*factor );
  if( opt->tolerance() != 0.000001*factor )
    {
    std::cout << "tolerance should be 0.000001 and not " << opt->tolerance()
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
  if( testNum == 1 )
    {
    opt->searchForMin( false );
    if( opt->searchForMin() )
      {
      std::cout << "searchForMin should be false and not true."
        << std::endl;
      returnStatus = EXIT_FAILURE;
      }
    }

  double x = 0.01 * factor;
  if( testNum != 0 )
    {
    x = opt->xMax();
    }
  double xVal = 0;
  if( !opt->extreme( &x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double idealX = - vnl_math::pi / 2;
  if( testNum != 0 )
    {
    idealX = vnl_math::pi * factor;
    }
  if( vnl_math_abs( idealX - x ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  x=" << x
      << " diff=" << vnl_math_abs( idealX - x )
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  double idealV = -1;
  if( vnl_math_abs( idealV - xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  delete opt;
  if( testNum != 0 )
    {
    delete myFunc2;
    }

  return returnStatus;
}
