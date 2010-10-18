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
#include "itkOptGoldenMean1D.h"
#include <cmath>
#include <iostream>

namespace itk 
{

OptGoldenMean1D::OptGoldenMean1D( )
: Optimizer1D( )
{
}

OptGoldenMean1D::OptGoldenMean1D( UserFunc<double, double> *newFuncVal )
: Optimizer1D( newFuncVal, NULL )
{
}


OptGoldenMean1D::~OptGoldenMean1D( )
{
}


void OptGoldenMean1D::use( UserFunc<double, double> * newFuncVal )
{
    Optimizer1D::use( newFuncVal, NULL );
}



bool OptGoldenMean1D::cExtreme( double *extX, double *extVal )
{  
  double maxSign = 1;
  if( cSearchForMin )
    {
    maxSign = -1;
    }

  double v;
  double prevV = cFuncVal->value( *extX );
  double xstep = cXStep;
  int dir = 1;
  double v1 = maxSign*cFuncVal->value( *extX+dir*xstep );
  double v2 = maxSign*cFuncVal->value( *extX-dir*xstep );
  
  if( v2>v1 ) 
    {
    dir = -1;
    v = v2;
    }
  else
    {
    v = v1;
    }
  
  unsigned int iter = 0;
  int dirInit = dir;
  while( v<prevV && xstep>2*cTolerance && iter<cMaxIterations )
    {
    dir *= -1;
    v = maxSign*cFuncVal->value( *extX+dir*xstep );
    if( v<prevV && dir == dirInit )
      {
      xstep /= 1.618;
      }
    ++iter;
    }

  prevV = v;

  *extX += dir*xstep;

  if( *extX>cXMax ) 
    {
    *extX = cXMax;
    *extVal = prevV;
    return false;
    }

  if( *extX<cXMin ) 
    {
    *extX = cXMin;
    *extVal = prevV;
    return false;
    }

  dirInit = dir;
  while( xstep>cTolerance && iter<cMaxIterations )
    {
    v = maxSign*cFuncVal->value( *extX+dir*xstep );
    while( v>prevV && iter<cMaxIterations ) 
      {
      dirInit = dir;
      prevV = v;
      *extX += dir*xstep;
      if( *extX>cXMax ) 
        {    
        *extX = cXMax;        
        *extVal = prevV;
        return false;
        }

      if( *extX<cXMin ) 
        {            
        *extX = cXMin;           
        *extVal = prevV;        
        return false;    
        } 
      v = maxSign*cFuncVal->value( *extX+dir*xstep );     
      ++iter;
      }  
    if( dir == dirInit )
      {
      xstep /= 1.618;
      }
    dir *= -1;
    } 

  *extVal = maxSign*prevV;
  return true;
}


}; // namespace itk

