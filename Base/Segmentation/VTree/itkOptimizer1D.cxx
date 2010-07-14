/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptimizer1D.cxx,v $
  Language:  C++
  Date:      $Date: 2005/09/25 15:27:51 $
  Version:   $Revision: 1.7 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkOptimizer1D.h"
#include <iostream>

#include <itkMacro.h>

namespace itk
{

bool 
Optimizer1D
::cExtreme(double * itkNotUsed(x), double * itkNotUsed(xVal))
{
  return false;
}


Optimizer1D
::Optimizer1D()
{
  cSearchForMin = true;
  cTolerance = 0.0001;
  cMaxIterations = 300;
  cXMin = 0;
  cXMax = 1;
  cXStep = 0.01;
  cDefined = false;
}



Optimizer1D::Optimizer1D(UserFunc<double, double> *newFuncVal, UserFunc<double, double> *newFuncDeriv)
{
  cSearchForMin = true;
  cTolerance = 0.0001;
  cMaxIterations = 300;
  cXMin = 0;
  cXMax = 1;
  cXStep = 0.01;
  cFuncVal = newFuncVal;
  cFuncDeriv = newFuncDeriv;
  cDefined = true;
}


Optimizer1D::~Optimizer1D()
{
}

void Optimizer1D::use(UserFunc<double, double> *newFuncVal, UserFunc<double, double> *newFuncDeriv)
{
  cFuncVal = newFuncVal;
  cFuncDeriv = newFuncDeriv;
  cDefined = true;
}


double Optimizer1D::xMin(void)
{
  return cXMin;
}


void Optimizer1D::xMin(double newXMin)
{
  cXMin = newXMin;
}

double Optimizer1D::xMax(void)
{
  return cXMax;
}


void Optimizer1D::xMax(double newXMax)
{
  cXMax = newXMax;
}

double Optimizer1D::xStep(void)
{
  return cXStep;
}


void Optimizer1D::xStep(double newXStep)
{
  cXStep = newXStep;
}


double Optimizer1D::tolerance(void)
{
  return cTolerance;
}

void Optimizer1D::tolerance(double newTolerance)
{
  cTolerance = newTolerance;
}

unsigned int Optimizer1D::maxIterations(void)
{
    return cMaxIterations;
}

void Optimizer1D::maxIterations(unsigned int newMaxIterations)
{
  cMaxIterations = newMaxIterations;
}


bool Optimizer1D::searchForMin(void)
{
  return cSearchForMin;
}


void Optimizer1D::searchForMin(bool newSearchForMin)
{
  cSearchForMin = newSearchForMin;
}

bool Optimizer1D::extreme(double *x, double *xVal)
{
  if(!cDefined)
  {
    return false;
  }
  return cExtreme(x, xVal);
}


} // end namespace itk

