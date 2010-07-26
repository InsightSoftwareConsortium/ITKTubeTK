/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptGoldenMean1D.cxx,v $
  Language:  C++
  Date:      $Date: 2005/09/04 16:48:56 $
  Version:   $Revision: 1.6 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkOptGoldenMean1D.h"
#include <cmath>
#include <iostream>

namespace itk 
{

OptGoldenMean1D::OptGoldenMean1D()
: Optimizer1D()
{
}

OptGoldenMean1D::OptGoldenMean1D(UserFunc<double, double> *newFuncVal)
: Optimizer1D(newFuncVal, NULL)
{
}


OptGoldenMean1D::~OptGoldenMean1D()
{
}


void OptGoldenMean1D::use(UserFunc<double, double> * newFuncVal)
{
    Optimizer1D::use(newFuncVal, NULL);
}



bool OptGoldenMean1D::cExtreme(double *extX, double *extVal)
{  
  double maxSign = 1;
  if(cSearchForMin)
  {
    maxSign = -1;
  }

  double v;
  double prevV = cFuncVal->value(*extX);
  double xstep = cXStep;
  int dir = 1;
  double v1 = maxSign*cFuncVal->value(*extX+dir*xstep);
  double v2 = maxSign*cFuncVal->value(*extX-dir*xstep);
  
  if(v2>v1) 
  {
    dir = -1;
    v = v2;
  }
  else
  {
    v = v1;
  }
   
  while(v<prevV && xstep>2*cTolerance)
  {
    dir *= -1;
    v = maxSign*cFuncVal->value(*extX+dir*xstep);
    if(v<prevV && dir == 1)
    {
      xstep /= 1.618;
    }
  }

  prevV = v;

  *extX += dir*xstep;

  if(*extX>cXMax) 
  {
    *extX = cXMax;
    *extVal = prevV;
    return false;
  }

  if(*extX<cXMin) 
  {
    *extX = cXMin;
    *extVal = prevV;
    return false;
  }

  while(xstep>cTolerance)
  {
    v = maxSign*cFuncVal->value(*extX+dir*xstep);
    while(v>prevV) 
    {
      prevV = v;
      *extX += dir*xstep;
      if(*extX>cXMax) 
      {    
        *extX = cXMax;        
        *extVal = prevV;
        return false;
      }

      if(*extX<cXMin) 
      {            
        *extX = cXMin;           
        *extVal = prevV;        
        return false;    
      } 
      v = maxSign*cFuncVal->value(*extX+dir*xstep);     
    }  
    xstep /= 1.618;
  } 
 *extVal = maxSign*prevV;
  return true;

}


}; // namespace itk

