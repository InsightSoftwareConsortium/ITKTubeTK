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

#include "tubeBrentOptimizer1D.h"
#include "tubeOptimizerND.h"

class MyNDFunc : public tube::UserFunction< vnl_vector< double >, double >
{
public:
  MyNDFunc( void )
    {
    cVal = 0;
    }
  const double & Value( const vnl_vector<double> & x )
    {
    cVal = std::sin( x( 0 ) ) + std::sin( x( 1 ) );
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
    cDx.set_size( 2 );
    }
  const vnl_vector<double> & Value( const vnl_vector<double> & x )
    {
    cDx[0] = std::cos( x( 0 ) );
    cDx[1] = std::cos( x( 1 ) );
    std::cout << "deriv = " << cDx[0] << ", " << cDx[1] << std::endl;
    return cDx;
    }

private:
  vnl_vector<double> cDx;

}; // End class MyNDFuncD

int tubeBrentOptimizerNDTest( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  double epsilon = 0.000001;

  MyNDFunc myFunc;
  MyNDFuncD myFuncD;
  tube::BrentOptimizer1D opt1D;

  tube::OptimizerND opt( 2, &myFunc, &myFuncD, &opt1D );

  typedef vnl_vector< double > VectorType;

  int returnStatus = EXIT_SUCCESS;

  VectorType xMin( 2 );
  xMin.fill( -3.5 );
  opt.SetXMin( xMin );
  if( opt.GetXMin()[0] != -3.5 )
    {
    std::cout << "xMin should be -3.5 and not " << opt.GetXMin()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType xMax( 2 );
  xMax.fill( 3.5 );
  opt.SetXMax( xMax );
  if( opt.GetXMax()[0] != 3.5 )
    {
    std::cout << "xMax should be 3.5 and not " << opt.GetXMax()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType xStep( 2 );
  xStep.fill( 0.1 );
  opt.SetXStep( xStep );
  if( opt.GetXStep()[0] != 0.1 )
    {
    std::cout << "xStep should be 0.1 and not " << opt.GetXStep()[0]
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.SetTolerance( 0.000001 );
  if( opt.GetTolerance() != 0.000001 )
    {
    std::cout << "tolerance should be 0.001 and not " << opt.GetTolerance()
      << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.SetMaxIterations( 300 );
  if( opt.GetMaxIterations() != 300 )
    {
    std::cout << "maxIterations should be 100 and not "
      << opt.GetMaxIterations() << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  opt.SetSearchForMin( true );
  if( !opt.GetSearchForMin() )
    {
    std::cout << "searchForMin should be false and not true." << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType x( 2 );
  x[0] = 0.01;
  x[1] = 1.01;
  double xVal = 0;
  if( !opt.Extreme( x, &xVal ) )
    {
    std::cout << "Optimization failed!" << std::endl;
    returnStatus = EXIT_FAILURE;
    }

  VectorType idealX( 2 );
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
