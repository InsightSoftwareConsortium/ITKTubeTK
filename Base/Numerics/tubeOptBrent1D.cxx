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
#include "tubeOptBrent1D.h"
#include "vnl/vnl_math.h"
#include <cmath>
#include <iostream>

namespace tube
{

OptBrent1D::OptBrent1D( void )
: Optimizer1D()
{
  cSmall = 1.0e-20;
}

OptBrent1D::OptBrent1D( UserFunc< double, double > * newFuncVal,
                        UserFunc< double, double > * newFuncDeriv )
: Optimizer1D( newFuncVal, newFuncDeriv )
{
  cSmall = 1.0e-20;
}


OptBrent1D::~OptBrent1D( void )
{
}


double OptBrent1D::smallDouble( void )
{
  return cSmall;
}

void OptBrent1D::smallDouble( double newSmall )
{
  cSmall = newSmall;
}

void OptBrent1D::use( UserFunc< double, double > * newFuncVal,
                      UserFunc< double, double > * newFuncDeriv )
{
  Optimizer1D::use( newFuncVal, newFuncDeriv );
}


void OptBrent1D::cMove( double & a, double & b, double & c,
  double d, double e, double f )
{
   a = d;
   b = e;
   c = f;
}

bool OptBrent1D::cExtreme( double *extX, double *extVal )
{
  unsigned int iter;

  int ok1, ok2;
  double a, b, d, d1, d2, du, dv, dw, dx, e=0.0;
  double fu, fv, fw, fx, olde, tol1, tol2, u, u1, u2, v, w, x, xm;

  double maxSign = -1;

  if( cSearchForMin )
    {
    maxSign = 1;
    }

  d = -1;
  v = *extX;
  fv = maxSign*cFuncVal->value(v);
  x = v+d*cXStep;
  if(x<cXMin || x>cXMax)
    {
    d *= -1;
    x = v+d*cXStep;
    }
  fx = maxSign*cFuncVal->value(x);
  if(fx>fv)
    {
    d *= -1;
    x = v+d*cXStep;
    fx = maxSign*cFuncVal->value(x);
    }
  w = 1;

  while(fx < fv)
    {
    v = x;
    fv = fx;
    x = v+d*cXStep*w;
    if(x<cXMin || x>cXMax)
      {
      if(x<cXMin)
        {
        x = cXMin;
        }
      else
        {
        x = cXMax;
        }
      fx = maxSign*cFuncVal->value(x);
      if( fx >= fv )
        {
        *extX = v;
        *extVal = maxSign*fv;
        //std::cout << " limit: x0 = " << x << std::endl;
        //std::cout << " limit: v = " << v << std::endl;
        //std::cout << " limit: vVal = " << maxSign*fv << std::endl;
        return false;
        }
      }
    else
      {
      fx = maxSign*cFuncVal->value(x);
      //std::cout << " x0 = " << x << std::endl;
      //std::cout << " xVal = " << fx << std::endl;
      w *= 1.1;
      }
    }

  u = v-d*cXStep*w;

  a = (u < x ? u : x);
  b = (u > x ? u : x);

  w = x = v;
  fw = fv = fx = maxSign*cFuncVal->value(v);
  dw = dv = dx = maxSign*cFuncDeriv->value(v);

  for(iter = 0; iter < cMaxIterations; iter++)
    {
    xm = 0.5 * (a+b);
    //std::cout << "x = " << x << std::endl;
    //std::cout << "  fx = " << fx << std::endl;
    //std::cout << "  a = " << a << std::endl;
    //std::cout << "  b = " << b << std::endl;
    tol1 = cTolerance * fabs(x) + cSmall;
    tol2 = 2.0 * tol1;
    if(fabs(x-xm) <= (tol2 - 0.5*(b-a)))
      {
      *extX = x;
      *extVal = maxSign*fx;
      return true;
      }
    if(fabs(e) > tol1)
      {
      d1 = 2.0*(b-a);
      d2 = d1;
      if(dw != dx)
        {
        d1 = (w-x)*dx/(dx-dw);
        }

      if(dv != dx)
        {
        d2 = (v-x)*dx/(dx-dv);
        }

      u1 = x + d1;
      u2 = x + d2;
      ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
      ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
      olde = e;
      e = d;
      if(ok1 || ok2)
        {
        if(ok1 && ok2)
          {
          d = (fabs(d1) < fabs(d2) ? d1 : d2);
          }
        else
          {
          if(ok1)
            {
            d = d1;
            }
          else
            {
            d = d2;
            }
          }

        if(fabs(d) <= fabs(0.5 * olde))
          {
          u = x+d;
          if(u-a < tol2 || b-u < tol2)
            {
            d = tol1 * vnl_math_sgn(xm-x);
            }
          }
        else
          {
          d = (double)0.5 * (e = (dx >= 0.0 ? a-x : b-x));
          }
        }
      else
        {
        d = (double)0.5 * (e = (dx >= 0.0 ? a-x : b-x));
        }
      }
    else
      {
      d = (double)0.5 * (e = (dx >= 0.0 ? a-x : b-x));
      }
    if(fabs(d) >= tol1)
      {
      u = x + d;
      fu = maxSign*cFuncVal->value(u);
      }
    else
      {
      u = x + tol1 * vnl_math_sgn(d);
      fu = maxSign*cFuncVal->value(u);
      if(fu > fx)
        {
        *extX = x;
        *extVal = maxSign*fx;
        return true;
        }
      }
    du = maxSign*cFuncDeriv->value(u);
    if(fu <= fx)
      {
      if(u >= x)
        {
        a = x;
        }
      else
        {
        b = x;
        }

      cMove(v,fv,dv, w,fw,dw);
      cMove(w,fw,dw, x,fx,dx);
      cMove(x,fx,dx, u,fu,du);
      }
    else
      {
      if(u < x)
        {
        a = u;
        }
      else
        {
        b = u;
        }

      if(fu <= fw || w == x)
        {
        cMove(v,fv,dv, w,fw,dw);
        cMove(w,fw,dw, u,fu,du);
        }
      else
        {
        if(fu < fv || v == x || v == w)
          {
          cMove(v,fv,dv, u,fu,du);
          }
        }
      }
    }
  *extX = x;
  *extVal = maxSign*fx;

  std::cout << "Reached max iterations." << std::endl;

  return false;
}


}; // namespace tube
