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

#include "itkMacro.h"

#include "../itkOptBrent1D.h"
#include "../itkOptimizer1D.h"
#include "../UserFunc.h"

class MyFunc:
  public itk::UserFunc< double, double > 
  {
  public:
    MyFunc( )
      {
      };
    double value( double x )
      {
      return vcl_sin(x);
      };
  };

int itkOptBrentTest( int itkNotUsed(argc), char **itkNotUsed(argv) )
{
  double epsilon = 0.000001;

  MyFunc * myFunc = new MyFunc();
  itk::OptBrent1D * opt = new itk::OptBrent1D( myFunc );

  opt->smallDouble( epsilon );

  MyFunc * myFunc2 = new MyFunc();
  opt->use( myFunc2 );

  delete myFunc;

  int returnStatus = EXIT_SUCCESS;

  opt->xMin( -1 );
  if( opt->xMin() != -1 )
    {
    std::cout << "xMin should be -1 and not " << opt->xMin() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->xMax( 1 );
  if( opt->xMax() != 1 )
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

  opt->tolerance( 0.001 );
  if( opt->tolerance() != 0.001 )
    {
    std::cout << "tolerance should be 0.001 and not " << opt->tolerance() 
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->maxIterations( 100 );
  if( opt->maxIterations() != 100 )
    {
    std::cout << "maxIterations should be 100 and not " 
      << opt->maxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt->searchForMin( false );
  if( opt->searchForMin() )
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

  if( vnl_math_abs( x ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  x=" << x 
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( vnl_math_abs( 1-xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal 
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  delete opt;
  delete myFunc2;

  return returnStatus;
}

