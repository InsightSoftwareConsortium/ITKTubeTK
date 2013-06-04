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

#include "tubeOptimizerND.h"

#include <iostream>

#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

namespace tube
{

class OptValFuncND : public UserFunction<double, double>
{
public:
  OptValFuncND( OptimizerND * newOpt )
    {
    m_Opt = newOpt;
    m_Val = 0;
    }

  const double & value( const double & x )
    {
    m_Val = m_Opt->funcVal( x );
    return m_Val;
    }

private:
  OptimizerND * m_Opt;
  double        m_Val;

}; // End class OptValFuncND

class OptDerivFuncND : public UserFunction<double, double>
{
public:
  OptDerivFuncND( OptimizerND * newOpt )
    {
    m_Opt = newOpt;
    m_Deriv = 0;
    }

  const double & value( const double & x )
    {
    m_Deriv = m_Opt->funcDeriv( x );
    return m_Deriv;
    }

private:
  OptimizerND * m_Opt;
  double m_Deriv;

}; // End class OptDerivFuncND


OptimizerND::OptimizerND( void )
{
  m_NDims = 0;

  m_SearchForMin = false;
  m_Tolerance = 0.0001;

  m_MaxIterations = 300;
  m_MaxLineSearches = 10;

  m_Opt1D = NULL;
  m_FuncValND = NULL;
  m_FuncDerivND = NULL;

  m_Opt1DVal = new OptValFuncND( this );
  m_Opt1DDeriv = new OptDerivFuncND( this );
}


OptimizerND::OptimizerND( int newNDims,
  UserFunction< vnl_vector<double>, double > * newFuncValND,
  UserFunction< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
  Optimizer1D * newOpt1D )
{
  m_NDims = 0;

  m_SearchForMin = false;
  m_Tolerance = 0.0001;

  m_MaxIterations = 300;
  m_MaxLineSearches = 10;

  m_Opt1D = NULL;
  m_FuncValND = NULL;
  m_FuncDerivND = NULL;

  m_Opt1DVal = new OptValFuncND( this );
  m_Opt1DDeriv = new OptDerivFuncND( this );

  this->use( newNDims, newFuncValND, newFuncDerivND, newOpt1D );
}

OptimizerND::~OptimizerND( void )
{
  delete m_Opt1DVal;
  delete m_Opt1DDeriv;
}

void OptimizerND::use( int newNDims,
  UserFunction< vnl_vector<double>, double > * newFuncValND,
  UserFunction< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
  Optimizer1D * newOpt1D )
{
  m_NDims = newNDims;

  m_XMin.set_size( m_NDims );
  m_XMin.fill( 0.0 );
  m_XMax.set_size( m_NDims );
  m_XMax.fill( 1.0 );
  m_XStep.set_size( m_NDims );
  m_XStep.fill( 0.01 );
  m_X0.set_size( m_NDims );
  m_X0.fill( 0.0 );
  m_X0Dir.set_size( m_NDims );
  m_X0Dir.fill( 1.0/vcl_sqrt((float)m_NDims) );
  m_X0Temp.set_size( m_NDims );
  m_X0Temp.fill( 0.0 );

  m_FuncValND = newFuncValND;
  m_FuncDerivND = newFuncDerivND;

  m_Opt1D = newOpt1D;
  if( m_Opt1D != NULL )
    {
    m_Opt1D->use( m_Opt1DVal, m_Opt1DDeriv );
    m_Opt1D->searchForMin( m_SearchForMin );
    m_Opt1D->tolerance( m_Tolerance );
    m_Opt1D->maxIterations( m_MaxIterations );
    }
}


vnl_vector<double> & OptimizerND::xMin( void )
{
  return m_XMin;
}

void OptimizerND::xMin( vnl_vector<double> & newXMinn )
{
  m_XMin = newXMinn;
}

vnl_vector<double> & OptimizerND::xMax( void )
{
  return m_XMax;
}

void OptimizerND::xMax( vnl_vector<double> & newXMaxx )
{
  m_XMax = newXMaxx;
}

vnl_vector<double> & OptimizerND::xStep( void )
{
  return m_XStep;
}

void OptimizerND::xStep( vnl_vector<double> & newXStepp )
{
  m_XStep = newXStepp;
}

double OptimizerND::tolerance( void )
{
  return m_Tolerance;
}


void OptimizerND::tolerance( double newTolerance )
{
  if( m_Opt1D != NULL )
    {
    m_Opt1D->tolerance( newTolerance );
    }
  m_Tolerance = newTolerance;
}

unsigned int OptimizerND::maxIterations( void )
{
  return m_MaxIterations;
}

void OptimizerND::maxIterations( unsigned int newMaxIterations )
{
  if( m_Opt1D != NULL )
    {
    m_Opt1D->maxIterations( newMaxIterations );
    }
  m_MaxIterations = newMaxIterations;
}

unsigned int OptimizerND::maxLineSearches( void )
{
  return m_MaxLineSearches;
}

void OptimizerND::maxLineSearches( unsigned int newMaxLineSearches )
{
  m_MaxLineSearches = newMaxLineSearches;
}


bool OptimizerND::searchForMin( void )
{
  return m_SearchForMin;
}


void OptimizerND::searchForMin( bool newSearchForMin )
{
  if( m_Opt1D != NULL )
    {
    m_Opt1D->searchForMin( newSearchForMin );
    }
  m_SearchForMin = newSearchForMin;
}


double OptimizerND::funcVal(double a )
{
  m_X0Temp = ComputeLineStep( m_X0, a, m_X0Dir );

  return m_FuncValND->value( m_X0Temp );
}


double OptimizerND::funcDeriv( double a )
{
  m_X0Temp = ComputeLineStep( m_X0, a, m_X0Dir );

  return dot_product( m_FuncDerivND->value( m_X0Temp ), m_X0Dir );
}


bool OptimizerND::extreme( vnl_vector<double> & x, double * xVal )
{
  m_X0 = x;
  double a = 1;

  double xmin;
  double xmax;

  unsigned int count = 0;
  while( vnl_math_abs(a) > m_Tolerance )
    {
    m_X0Dir = m_FuncDerivND->value( m_X0 );
    if( m_X0Dir.magnitude() < m_Tolerance * m_Tolerance )
      {
      a = 0;
      break;
      }

    m_X0Dir.normalize();

    if( m_X0Dir(0) != 0 )
      {
      xmin = (m_XMin(0) - m_X0(0)) / m_X0Dir(0);
      xmax = (m_XMax(0) - m_X0(0)) / m_X0Dir(0);
      }
    else
      {
      xmin = 0;
      xmax = 0;
      }
    if( xmin > xmax )
      {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      }
    for( unsigned int i=1; i<m_NDims; i++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( m_X0Dir(i) != 0 )
        {
        tmin = (m_XMin(i) - m_X0(i))/m_X0Dir(i);
        tmax = (m_XMax(i) - m_X0(i))/m_X0Dir(i);
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

    //std::cout << "  x = " << m_X0 << std::endl;
    //std::cout << "  dir = " << m_X0Dir << std::endl;
    //std::cout << "  xmin = " << xmin << std::endl;
    //std::cout << "  xmax = " << xmax << std::endl;

    if( xmin == xmax )
      {
      a = 0;
      m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
      break;
      }

    double xstep = vnl_math_abs( dot_product(m_XStep, m_X0Dir) );

    a = 0;
    m_Opt1D->xMin( xmin );
    m_Opt1D->xMax( xmax );
    m_Opt1D->xStep( xstep );

    m_Opt1D->extreme( &a, xVal );

    m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
    if( count++ > m_MaxLineSearches )
      {
      break;
      }
    }

  x = m_X0;

  if( vnl_math_abs(a) > m_Tolerance )
    {
    std::cout << "Scanned " << count
      << " directions without convergence - aborting" << std::endl;
    return false;
    }

  return true;
}


bool OptimizerND::extreme( vnl_vector<double> & x, double * xVal,
  unsigned int n, MatrixType &dirs )
{
  m_X0 = x;
  double a;
  double xmin, xmax, xstep;
  for(unsigned int i=0; i<m_MaxLineSearches; i++ )
    {
    for(unsigned int j=0; j<x.size(); j++ )
      {
      m_X0Dir( j ) = dirs.get( j, i%n );
      }

    m_X0Dir.normalize();

    if( m_X0Dir(0) != 0 )
      {
      xmin = (m_XMin(0) - m_X0(0)) / m_X0Dir(0);
      xmax = (m_XMax(0) - m_X0(0)) / m_X0Dir(0);
      }
    else
      {
      xmin = 0;
      xmax = 0;
      }
    if( xmin > xmax )
      {
      double tmp = xmin;
      xmin = xmax;
      xmax = tmp;
      }
    for( unsigned int k=1; k<m_NDims; k++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( m_X0Dir(0) != 0 )
        {
        tmin = (m_XMin(k) - m_X0(k))/m_X0Dir(k);
        tmax = (m_XMax(k) - m_X0(k))/m_X0Dir(k);
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

    xstep = vnl_math_abs( dot_product( m_XStep, m_X0Dir) );

    a = 0;

    m_Opt1D->xMin( xmin );
    m_Opt1D->xMax( xmax );
    m_Opt1D->xStep( xstep );

    m_Opt1D->extreme( &a, xVal );
    m_X0 = ComputeLineStep( m_X0, a, m_X0Dir );
    }

  x = m_X0;

  return true;
}

} // End namespace tube
