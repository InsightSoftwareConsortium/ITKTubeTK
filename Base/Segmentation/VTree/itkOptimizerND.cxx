/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#include "itkOptimizerND.h"
#include "itkOptimizer1D.h"
#include "itkUserFunc.h"
#include <iostream>
#include <cmath>

namespace itk
{


class OptValFuncND : public UserFunc<double, double> 
{

private:

  OptimizerND * cOpt;
  double cVal;

public:

  OptValFuncND( OptimizerND * newOpt )
  {
    cOpt = newOpt;
    cVal = 0;
  };

  const double & value( const double & x )
  {
    cVal = cOpt->funcVal( x );
    return cVal;
  };

};

class OptDerivFuncND : public UserFunc<double, double> 
{

private:

  OptimizerND * cOpt;
  double cDeriv;

public:

  OptDerivFuncND( OptimizerND * newOpt )
  {
    cOpt = newOpt;
    cDeriv = 0;
  };

  const double & value( const double & x )
  {
    cDeriv = cOpt->funcDeriv( x );
    return cDeriv;
  };

};


OptimizerND::OptimizerND( void )
{
  cNDims = 0;

  cSearchForMin = false;
  cTolerance = 0.0001;

  cMaxIterations = 300;
  cMaxLineSearches = 10;

  cOpt1D = NULL;
  cFuncValND = NULL;
  cFuncDerivND = NULL;

  cOpt1DVal = new OptValFuncND( this );
  cOpt1DDeriv = new OptDerivFuncND( this );
}


OptimizerND::OptimizerND( int newNDims,
  UserFunc< vnl_vector<double>, double > * newFuncValND,
  UserFunc< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
  Optimizer1D * newOpt1D )
{
  cNDims = 0;

  cSearchForMin = false;
  cTolerance = 0.0001;

  cMaxIterations = 300;
  cMaxLineSearches = 10;

  cOpt1D = NULL;
  cFuncValND = NULL;
  cFuncDerivND = NULL;

  cOpt1DVal = new OptValFuncND( this );
  cOpt1DDeriv = new OptDerivFuncND( this );

  this->use( newNDims, newFuncValND, newFuncDerivND, newOpt1D );
}

OptimizerND::~OptimizerND( void )
{
  delete cOpt1DVal;
  delete cOpt1DDeriv;
}

void OptimizerND::use( int newNDims,
  UserFunc< vnl_vector<double>, double > * newFuncValND,
  UserFunc< vnl_vector<double>, vnl_vector<double> > * newFuncDerivND,
  Optimizer1D * newOpt1D )
{
  cNDims = newNDims;

  cXMin.set_size( cNDims );
  cXMin.fill( 0.0 );
  cXMax.set_size( cNDims );
  cXMax.fill( 1.0 );
  cXStep.set_size( cNDims );
  cXStep.fill( 0.01 );
  cX0.set_size( cNDims );
  cX0.fill( 0.0 );
  cX0Dir.set_size( cNDims );
  cX0Dir.fill( 1.0/vcl_sqrt((float)cNDims) );
  cX0Temp.set_size( cNDims );
  cX0Temp.fill( 0.0 );

  cFuncValND = newFuncValND;
  cFuncDerivND = newFuncDerivND;

  cOpt1D = newOpt1D;
  if( cOpt1D != NULL )
    {
    cOpt1D->use( cOpt1DVal, cOpt1DDeriv );
    cOpt1D->searchForMin( cSearchForMin );
    cOpt1D->tolerance( cTolerance );
    cOpt1D->maxIterations( cMaxIterations );
    }
}


vnl_vector<double> & OptimizerND::xMin( void )
{
  return cXMin;
}

void OptimizerND::xMin( vnl_vector<double> & newXMinn )
{
   cXMin = newXMinn;
}

vnl_vector<double> & OptimizerND::xMax( void )
{
  return cXMax;
}

void OptimizerND::xMax( vnl_vector<double> & newXMaxx )
{
  cXMax = newXMaxx;
}

vnl_vector<double> & OptimizerND::xStep( void )
{
  return cXStep;
}

void OptimizerND::xStep( vnl_vector<double> & newXStepp )
{
  cXStep = newXStepp;
}

double OptimizerND::tolerance( void )
{
  return cTolerance;
}


void OptimizerND::tolerance( double newTolerance )
{
  if( cOpt1D != NULL )
    {
    cOpt1D->tolerance( newTolerance );
    }
  cTolerance = newTolerance;
}

unsigned int OptimizerND::maxIterations( void )
{
  return cMaxIterations;
}

void OptimizerND::maxIterations( unsigned int newMaxIterations )
{
  if( cOpt1D != NULL )
    {
    cOpt1D->maxIterations( newMaxIterations );
    }
  cMaxIterations = newMaxIterations;
}

unsigned int OptimizerND::maxLineSearches( void )
{
  return cMaxLineSearches;
}

void OptimizerND::maxLineSearches( unsigned int newMaxLineSearches )
{
  cMaxLineSearches = newMaxLineSearches;
}


bool OptimizerND::searchForMin( void )
{
  return cSearchForMin;
}


void OptimizerND::searchForMin( bool newSearchForMin )
{
  if( cOpt1D != NULL )
    {
    cOpt1D->searchForMin( newSearchForMin );
    }
  cSearchForMin = newSearchForMin;
}


double OptimizerND::funcVal(double a )
{
  cX0Temp = ComputeLineStep( cX0, a, cX0Dir );

  return cFuncValND->value( cX0Temp );
}


double OptimizerND::funcDeriv( double a )
{
  cX0Temp = ComputeLineStep( cX0, a, cX0Dir );

  return dot_product( cFuncDerivND->value( cX0Temp ), cX0Dir );
}


bool OptimizerND::extreme( vnl_vector<double> & x, double * xVal )
{
  cX0 = x;
  double a = 1;
  double xmin, xmax, xstep;
  unsigned int count = 0;
  while( fabs(a) > cTolerance ) 
    {
    cX0Dir = cFuncDerivND->value( cX0 );
    if( cX0Dir.magnitude() < cTolerance*cTolerance ) 
      {
      a = 0;
      break;
      }

    cX0Dir.normalize();

    if( cX0Dir(0) != 0 )
      {
      xmin = (cXMin(0) - cX0(0)) / cX0Dir(0);
      xmax = (cXMax(0) - cX0(0)) / cX0Dir(0);
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
    for( unsigned int i=1; i<cNDims; i++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( cX0Dir(i) != 0 )
        {
        tmin = (cXMin(i) - cX0(i))/cX0Dir(i);
        tmax = (cXMax(i) - cX0(i))/cX0Dir(i);
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

    //std::cout << "  x = " << cX0 << std::endl;
    //std::cout << "  dir = " << cX0Dir << std::endl;
    //std::cout << "  xmin = " << xmin << std::endl;
    //std::cout << "  xmax = " << xmax << std::endl;
    
    if( xmin == xmax )
      {
      a = 0;
      cX0 = ComputeLineStep( cX0, a, cX0Dir );
      break;
      }

    xstep = fabs( dot_product(cXStep, cX0Dir) );

    a = 0;
    cOpt1D->xMin( xmin );
    cOpt1D->xMax( xmax );
    cOpt1D->xStep( xstep );

    cOpt1D->extreme( &a, xVal );

    cX0 = ComputeLineStep( cX0, a, cX0Dir );
    if( count++ > cMaxLineSearches )
      {
      break;
      }
    }

  x = cX0;

  if( fabs(a) > cTolerance )
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
  cX0 = x;
  double a;
  double xmin, xmax, xstep;
  for(unsigned int i=0; i<cMaxLineSearches; i++ )
    {
    for(unsigned int j=0; j<x.size(); j++ )
      {
      cX0Dir( j ) = dirs.get( j, i%n );
      }

    cX0Dir.normalize();

    if( cX0Dir(0) != 0 )
      {
      xmin = (cXMin(0) - cX0(0)) / cX0Dir(0);
      xmax = (cXMax(0) - cX0(0)) / cX0Dir(0);
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
    for( unsigned int k=1; k<cNDims; k++ )
      {
      double tmin = xmin;
      double tmax = xmax;
      if( cX0Dir(0) != 0 )
        {
        tmin = (cXMin(k) - cX0(k))/cX0Dir(k);
        tmax = (cXMax(k) - cX0(k))/cX0Dir(k);
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
      cX0 = ComputeLineStep( cX0, a, cX0Dir );
      continue;
      }

    xstep = fabs( dot_product(cXStep, cX0Dir) );

    a = 0;

    cOpt1D->xMin( xmin );
    cOpt1D->xMax( xmax );
    cOpt1D->xStep( xstep );

    cOpt1D->extreme( &a, xVal );
    cX0 = ComputeLineStep( cX0, a, cX0Dir );
    }

  x = cX0;

  return true;
}

}; // namespace ITK

