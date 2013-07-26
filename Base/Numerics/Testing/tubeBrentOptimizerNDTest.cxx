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

#include "tubeBrentOptimizer1D.h"
#include "tubeOptimizer1D.h"
#include "tubeOptimizerND.h"
#include "tubeUserFunction.h"

#include <itkMacro.h>

#include <vcl_cmath.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>

#include <iostream>

class MyNDFunc : public tube::UserFunction< vnl_vector< double >, double >
{
public:
  MyNDFunc( void )
    {
    cVal = 0;
    }
  const double & value( const vnl_vector<double> & x )
    {
    cVal = vcl_sin(x(0)) + vcl_sin(x(1));
    return cVal;
    }

private:
  double cVal;

}; // End class MyNDFunc

class MyNDFuncD : public tube::UserFunction< vnl_vector< double >,
                                             vnl_vector< double > >
{
public:
  MyNDFuncD( void )
    {
    cDx.set_size(2);
    }
  const vnl_vector<double> & value( const vnl_vector<double> & x )
    {
    cDx[0] = vcl_cos(x(0));
    cDx[1] = vcl_cos(x(1));
    std::cout << "deriv = " << cDx[0] << ", " << cDx[1] << std::endl;
    return cDx;
    }

private:
  vnl_vector<double> cDx;

}; // End class MyNDFuncD

int tubeBrentOptimizerNDTest( int itkNotUsed( argc ), char * itkNotUsed( argv )[] )
{
  double epsilon = 0.000001;

  MyNDFunc myFunc;
  MyNDFuncD myFuncD;
  tube::BrentOptimizer1D opt1D;

  tube::OptimizerND opt( 2, &myFunc, &myFuncD, &opt1D );

  typedef vnl_vector< double > VectorType;

  int returnStatus = EXIT_SUCCESS;

  VectorType xMin(2);
  xMin.fill( -3.5 );
  opt.xMin( xMin );
  if( opt.xMin()[0] != -3.5 )
    {
    std::cout << "xMin should be -3.5 and not " << opt.xMin()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType xMax(2);
  xMax.fill( 3.5 );
  opt.xMax( xMax );
  if( opt.xMax()[0] != 3.5 )
    {
    std::cout << "xMax should be 3.5 and not " << opt.xMax()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType xStep(2);
  xStep.fill( 0.1 );
  opt.xStep( xStep );
  if( opt.xStep()[0] != 0.1 )
    {
    std::cout << "xStep should be 0.1 and not " << opt.xStep()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.tolerance( 0.000001 );
  if( opt.tolerance() != 0.000001 )
    {
    std::cout << "tolerance should be 0.001 and not " << opt.tolerance()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.maxIterations( 300 );
  if( opt.maxIterations() != 300 )
    {
    std::cout << "maxIterations should be 100 and not "
      << opt.maxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.searchForMin( true );
  if( !opt.searchForMin() )
    {
    std::cout << "searchForMin should be false and not true." << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType x(2);
  x[0] = 0.01;
  x[1] = 1.01;
  double xVal = 0;
  if( !opt.extreme( x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType idealX(2);
  idealX[0] = - vnl_math::pi / 2;
  idealX[1] = - vnl_math::pi / 2;
  double diff0 = idealX[0] - x[0];
  double diff1 = idealX[1] - x[1];
  double diff = diff0*diff0 + diff1*diff1;
  if( diff > epsilon )
    {
    std::cout << "Optimization not within tolerance!" << std::endl;
    std::cout << "  x=" << x << std::endl;
    std::cout << "  idealX=" << idealX << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  if( vnl_math_abs( -2 - xVal ) > epsilon )
    {
    std::cout << "Optimization not within tolerance!  xVal=" << xVal
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  return returnStatus;
}
