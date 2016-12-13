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

#include "tubeOptimizerND.h"

namespace tube
{

class OptimizerNDValueFunction : public UserFunction< double, double >
{
public:

  typedef OptimizerNDValueFunction        Self;
  typedef UserFunction< double, double >  Superclass;
  typedef Self *                          Pointer;
  typedef const Self *                    ConstPointer;

  typedef double                    InputType;
  typedef double                    OutputType;

  tubeTypeMacro( OptimizerNDValueFunction );

  OptimizerNDValueFunction( OptimizerND::Pointer optimizer )
    {
    m_Optimizer = optimizer;
    m_Value = 0.0;
    }

  const OutputType & Value( const InputType & input )
    {
    m_Value = m_Optimizer->FuncVal( input );
    return m_Value;
    }

private:

  OptimizerNDValueFunction( const Self & self );
  void operator=( const Self & self );

  OptimizerND::Pointer  m_Optimizer;
  OutputType            m_Value;

}; // End class OptimizerNDValueFunction

class OptimizerNDDerivativeFunction : public UserFunction< double, double >
{
public:

  typedef OptimizerNDDerivativeFunction    Self;
  typedef UserFunction< double, double >   Superclass;
  typedef Self *                           Pointer;
  typedef const Self *                     ConstPointer;

  typedef double                           InputType;
  typedef double                           OutputType;

  tubeTypeMacro( OptimizerNDDerivativeFunction );

  OptimizerNDDerivativeFunction( OptimizerND::Pointer optimizer )
    {
    m_Optimizer = optimizer;
    m_Derivative = 0.0;
    }

  const OutputType & Value( const InputType & input )
    {
    m_Derivative = m_Optimizer->FuncDeriv( input );
    return m_Derivative;
    }

private:

  OptimizerNDDerivativeFunction( const Self & self );
  void operator=( const Self & self );

  OptimizerND::Pointer  m_Optimizer;
  OutputType            m_Derivative;

}; // End class OptimizerNDDerivativeFunction


OptimizerND
::OptimizerND( void )
 : m_FuncValND( NULL ), m_FuncDerivND( NULL )
{
  m_Dimension = 0;

  m_SearchForMin = false;
  m_Tolerance = 0.0001;

  m_MaxIterations = 300;
  m_MaxLineSearches = 10;

  m_Optimizer1D = NULL;

  m_Optimizer1DVal = new OptimizerNDValueFunction( this );
  m_Optimizer1DDeriv = new OptimizerNDDerivativeFunction( this );
}


OptimizerND
::OptimizerND( unsigned int dimension,
  ValueFunctionType::Pointer funcValND,
  DerivativeFunctionType::Pointer funcDerivND,
  Optimizer1D::Pointer optimizer1D )
{
  m_Dimension = 0;

  m_SearchForMin = false;
  m_Tolerance = 0.0001;

  m_MaxIterations = 300;
  m_MaxLineSearches = 10;

  m_Optimizer1D = NULL;
  m_FuncValND = NULL;
  m_FuncDerivND = NULL;

  m_Optimizer1DVal = new OptimizerNDValueFunction( this );
  m_Optimizer1DDeriv = new OptimizerNDDerivativeFunction( this );

  this->Use( dimension, funcValND, funcDerivND, optimizer1D );
}


OptimizerND
::~OptimizerND( void )
{
  delete m_Optimizer1DVal;
  delete m_Optimizer1DDeriv;
}


void
OptimizerND
::Use( unsigned int dimension,
  ValueFunctionType::Pointer funcValND,
  DerivativeFunctionType::Pointer funcDerivND,
  Optimizer1D::Pointer optimizer1D )
{
  m_Dimension = dimension;

  m_XMin.set_size( m_Dimension );
  m_XMin.fill( 0.0 );
  m_XMax.set_size( m_Dimension );
  m_XMax.fill( 1.0 );
  m_XStep.set_size( m_Dimension );
  m_XStep.fill( 0.01 );
  m_X0.set_size( m_Dimension );
  m_X0.fill( 0.0 );
  m_X0Dir.set_size( m_Dimension );
  m_X0Dir.fill( 1.0/std::sqrt( ( float )m_Dimension ) );
  m_X0Temp.set_size( m_Dimension );
  m_X0Temp.fill( 0.0 );

  m_FuncValND = funcValND;
  m_FuncDerivND = funcDerivND;

  m_Optimizer1D = optimizer1D;
  if( m_Optimizer1D != NULL )
    {
    m_Optimizer1D->Use( m_Optimizer1DVal, m_Optimizer1DDeriv );
    m_Optimizer1D->SetSearchForMin( m_SearchForMin );
    m_Optimizer1D->SetTolerance( m_Tolerance );
    m_Optimizer1D->SetMaxIterations( m_MaxIterations );
    }
}


void
OptimizerND
::SetTolerance( double tolerance )
{
  if( m_Optimizer1D != NULL )
    {
    m_Optimizer1D->SetTolerance( tolerance );
    }
  m_Tolerance = tolerance;
}


void
OptimizerND
::SetMaxIterations( unsigned int maxIterations )
{
  if( m_Optimizer1D != NULL )
    {
    m_Optimizer1D->SetMaxIterations( maxIterations );
    }
  m_MaxIterations = maxIterations;
}


void
OptimizerND
::SetSearchForMin( bool searchForMin )
{
  if( m_Optimizer1D != NULL )
    {
    m_Optimizer1D->SetSearchForMin( searchForMin );
    }
  m_SearchForMin = searchForMin;
}


double
OptimizerND
::FuncVal( double a )
{
  m_X0Temp = ComputeLineStep( m_X0, a, m_X0Dir );

  return m_FuncValND->Value( m_X0Temp );
}


double
OptimizerND
::FuncDeriv( double a )
{
  m_X0Temp = ComputeLineStep( m_X0, a, m_X0Dir );

  return dot_product( m_FuncDerivND->Value( m_X0Temp ), m_X0Dir );
}


bool
OptimizerND
::Extreme( VectorType & x, double * xVal )
{
  m_X0 = x;
  double a = 1;

  double xmin;
  double xmax;

  unsigned int count = 0;
  while( vnl_math_abs( a ) > m_Tolerance )
    {
    m_X0Dir = m_FuncDerivND->Value( m_X0 );
    if( m_X0Dir.magnitude() < m_Tolerance * m_Tolerance )
      {
      a = 0;
      break;
      }

    m_X0Dir.normalize();

    if( m_X0Dir( 0 ) != 0 )
      {
      xmin = ( m_XMin( 0 ) - m_X0( 0 ) ) / m_X0Dir( 0 );
      xmax = ( m_XMax( 0 ) - m_X0( 0 ) ) / m_X0Dir( 0 );
      }
    else
      {
      xmin = -9999999999;
      xmax = 9999999999;
      }
    if( xmin > xmax )
      {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      }
    for( unsigned int i=1; i<m_Dimension; i++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( m_X0Dir( i ) != 0 )
        {
        tmin = ( m_XMin( i ) - m_X0( i ) )/m_X0Dir( i );
        tmax = ( m_XMax( i ) - m_X0( i ) )/m_X0Dir( i );
        }
      if( tmin > tmax )
        {
        double tmp = tmin;
        tmin = tmax;
        tmax = tmp;
        }
      if( tmin > xmin )
        {
        xmin = tmin;
        }
      if( tmax < xmax )
        {
        xmax = tmax;
        }
      }

    if( xmin == xmax )
      {
      a = 0;
      m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
      break;
      }

    double xstep = vnl_math_abs( dot_product( m_XStep, m_X0Dir ) );

    a = 0;
    m_Optimizer1D->SetXMin( xmin );
    m_Optimizer1D->SetXMax( xmax );
    m_Optimizer1D->SetXStep( xstep );

    m_Optimizer1D->Extreme( &a, xVal );

    m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
    if( count++ > m_MaxLineSearches )
      {
      break;
      }
    }

  x = m_X0;

  if( vnl_math_abs( a ) > m_Tolerance )
    {
    std::cout << "Scanned " << count
      << " directions without convergence - aborting" << std::endl;
    std::cout << "  Last step size = " << a << std::endl;
    std::cout << "  Tolerance = " << m_Tolerance << std::endl;
    return false;
    }

  return true;
}

bool
OptimizerND
::Extreme( VectorType & x, double * xVal, unsigned int n, MatrixType & directions )
{
  m_X0 = x;
  double a;
  double xmin;
  double xmax;
  double xstep;
  for( unsigned int i=0; i<m_MaxLineSearches; i++ )
    {
    for( unsigned int j=0; j<x.size(); j++ )
      {
      m_X0Dir( j ) = directions.get( j, i%n );
      }

    m_X0Dir.normalize();

    if( m_X0Dir( 0 ) != 0 )
      {
      xmin = ( m_XMin( 0 ) - m_X0( 0 ) ) / m_X0Dir( 0 );
      xmax = ( m_XMax( 0 ) - m_X0( 0 ) ) / m_X0Dir( 0 );
      }
    else
      {
      xmin = -9999999999;
      xmax = 9999999999;
      }
    if( xmin > xmax )
      {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      }
    for( unsigned int k=1; k<m_Dimension; k++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( m_X0Dir( k ) != 0 )
        {
        tmin = ( m_XMin( k ) - m_X0( k ) )/m_X0Dir( k );
        tmax = ( m_XMax( k ) - m_X0( k ) )/m_X0Dir( k );
        }
      if( tmin > tmax )
        {
        double tmp = tmin;
        tmin = tmax;
        tmax = tmp;
        }
      if( tmin > xmin )
        {
        xmin = tmin;
        }
      if( tmax < xmax )
        {
        xmax = tmax;
        }
      }

    if( xmin == xmax )
      {
      a = 0;
      m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
      continue;
      }

    xstep = vnl_math_abs( dot_product( m_XStep, m_X0Dir ) );

    a = 0;

    m_Optimizer1D->SetXMin( xmin );
    m_Optimizer1D->SetXMax( xmax );
    m_Optimizer1D->SetXStep( xstep );

    m_Optimizer1D->Extreme( &a, xVal );
    m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
    }

  x = m_X0;

  return true;
}


void
OptimizerND
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Dimension:        " << m_Dimension << std::endl;
  os << indent << "XMin:             " << m_XMin << std::endl;
  os << indent << "XMax:             " << m_XMax << std::endl;
  os << indent << "XStep:            " << m_XStep << std::endl;
  os << indent << "X0:               " << m_X0 << std::endl;
  os << indent << "X0Dir:            " << m_X0Dir << std::endl;
  os << indent << "X0Temp:           " << m_X0Temp << std::endl;
  os << indent << "SearchForMin:     " << m_SearchForMin << std::endl;
  os << indent << "Tolerance:        " << m_Tolerance << std::endl;
  os << indent << "MaxIterations:    " << m_MaxIterations << std::endl;
  os << indent << "MaxLineSearches:  " << m_MaxLineSearches << std::endl;
  os << indent << "Optimizer1DVal:   " << m_Optimizer1DVal << std::endl;
  os << indent << "Optimizer1DDeriv: " << m_Optimizer1DDeriv << std::endl;
  os << indent << "Optimizer1D:      " << m_Optimizer1D << std::endl;
  os << indent << "FuncValND:        " << m_FuncValND << std::endl;
  os << indent << "FuncDerivND:      " << m_FuncDerivND << std::endl;
}

} // End namespace tube
