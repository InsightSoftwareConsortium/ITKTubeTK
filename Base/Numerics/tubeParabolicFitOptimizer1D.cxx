/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

namespace tube
{

ParabolicFitOptimizer1D
::ParabolicFitOptimizer1D( void )
  : Optimizer1D()
{
}


ParabolicFitOptimizer1D
::ParabolicFitOptimizer1D( ValueFunctionType::Pointer funcVal )
  : Optimizer1D( funcVal, NULL )
{
}


ParabolicFitOptimizer1D
::~ParabolicFitOptimizer1D( void )
{
}


void
ParabolicFitOptimizer1D
::Use( ValueFunctionType::Pointer funcVal )
{
  this->Superclass::Use( funcVal, NULL );
}


double
ParabolicFitOptimizer1D
::m_Center( double x1, double y1, double x2, double y2, double x3, double y3 )
{
  double a = ( y1 - ( ( y2-y3 )*x1 )/( x2-x3 ) - y3
      + ( ( y2-y3 )*x3 )/( x2-x3 ) ) /
    ( x1*x1 - x3*x3 + ( ( x3*x3-x2*x2 )*x1 )/( x2-x3 ) -
     ( ( x3*x3-x2*x2 )*x3 )/( x2-x3 ) );
  double b = ( y2 - a * x2*x2 - y3 + a * x3*x3 ) / ( x2 - x3 );

  return -b/( 2*a );
}


bool
ParabolicFitOptimizer1D
::m_Extreme( double * extX, double * extVal )
{
  double minSign = 1;
  if( !m_SearchForMin )
    {
    minSign = -1;
    }

  double d;
  double v;
  double fv;
  double x;
  double fx;
  double u;
  double fu;
  double w;

  d = -1;
  v = ( *extX );
  fv = minSign * m_FuncVal->Value( v );
  x = v + d * m_XStep;
  u = v;
  fu = fv;

  if( x < m_XMin || x > m_XMax )
    {
    d *= -1;
    x = v + d * m_XStep;
    }

  fx = minSign * m_FuncVal->Value( x );

  if( fx>=fv )
    {
    u = x;
    fu = fx;
    d *= -1;
    x = v + d * m_XStep;
    fx = minSign * m_FuncVal->Value( x );
    }
  w = 1;
  while( fx < fv )
    {
    u = v;
    fu = fv;
    v = x;
    fv = fx;
    x = v + d * m_XStep*w;
    while( ( x < m_XMin || x > m_XMax ) && w > m_Tolerance )
      {
      w /= 2;
      x = v + d * m_XStep*w;
      }

    if( x < m_XMin || x > m_XMax )
      {
      x = v - d * m_XStep * w;
      fx = minSign * m_FuncVal->Value( x );
      *extX = x;
      *extVal = minSign * fx;
      std::cout << " Exiting legal parameter space - aborting" << std::endl;
      std::cout << "   x = " << x << std::endl;
      std::cout << "   Range = " << m_XMin << " - " << m_XMax << std::endl;
      return false;
      }

    fx = minSign * m_FuncVal->Value( x );
    w *= 1.1;
    }

    // Bracket = u, v, x;
  if( x<u )
    {
    double tf = u;
    u = x;
    x = tf;
    tf = fu;
    fu = fx;
    fx = tf;
    }

  d = ( x-u )/2;


  while( d > m_Tolerance )
    {
    w = m_Center( u, fu, v, fv, x, fx );
    if( w>=x || w<=u )
      {
      *extX = v;
      *extVal = minSign*fv;
      std::cout << "w " << w << " ";
      std::cout << ": u " << u << "=" << fu << " ";
      std::cout << ": v " << v << "=" << fv << " ";
      std::cout << ": x " << x << "=" << fx << std::endl;
      std::cout << "d = " << d << std::endl;
      std::cout << " parabola outside bounds - aborting" << std::endl;
      return false;
      }

    double fw = minSign * m_FuncVal->Value( w );

    if( fw<fv )
      {
      if( w>v )
        {
        u = v;
        fu = fv;
        v = w;
        fv = fw;
        d = ( x-u )/2;
        }
      else // w<v
        {
        x = v;
        fx = fv;
        v = w;
        fv = fw;
        d = ( x-u )/2;
        }
      }
    else
      {
      if( w>v )
        {
        x = w;
        fx = fw;
        d = ( x-u )/2;
        }
      else // w<v
        {
        u = w;
        fu = fw;
        d = ( x-u )/2;
        }
      }

    if( x - v < m_Tolerance && d > m_Tolerance )
      {
      x = v;
      fx = fv;
      v = x-( x-u )/3;
      fv = minSign * m_FuncVal->Value( v );
      d = ( x-u )/2;
      while( fv > fx && d > m_Tolerance )
        {
        u = v;
        fu = fv;
        v = x-( x-u )/3;
        fv = minSign * m_FuncVal->Value( v );
        d = ( x-u )/2;
        }
      }
    else if( v-u < m_Tolerance && d > m_Tolerance )
      {
      u = v;
      fu = fv;
      v = u+( x-u )/3;
      fv = minSign * m_FuncVal->Value( v );
      d = ( x-u )/2;
      while( fv > fu && d > m_Tolerance )
        {
        x = v;
        fx = fv;
        v = u+( x-u )/3;
        fv = minSign * m_FuncVal->Value( v );
        d = ( x-u )/2;
        }
      }
    }

  *extX = v;
  *extVal = minSign*fv;

  return true;
}

} // End namespace tube
