/*=========================================================================

Library:   TubeTK/VTree3D

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
#include "itkOptimizer1D.h"
#include <iostream>

#include <itkMacro.h>

namespace itk
{

bool 
Optimizer1D
::cExtreme( double * itkNotUsed(x), double * itkNotUsed(xVal) )
{
  return false;
}


Optimizer1D
::Optimizer1D( void )
{
  cSearchForMin = true;
  cTolerance = 0.0001;
  cMaxIterations = 300;
  cXMin = 0;
  cXMax = 1;
  cXStep = 0.01;
  cDefined = false;
}



Optimizer1D::Optimizer1D( UserFunc< double, double > * newFuncVal,
  UserFunc< double, double > * newFuncDeriv )
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


Optimizer1D::~Optimizer1D( void )
{
}

void Optimizer1D::use( UserFunc< double, double > * newFuncVal,
  UserFunc< double, double > * newFuncDeriv )
{
  cFuncVal = newFuncVal;
  cFuncDeriv = newFuncDeriv;
  cDefined = true;
}


double Optimizer1D::xMin( void )
{
  return cXMin;
}


void Optimizer1D::xMin( double newXMin )
{
  cXMin = newXMin;
}

double Optimizer1D::xMax( void )
{
  return cXMax;
}


void Optimizer1D::xMax( double newXMax )
{
  cXMax = newXMax;
}

double Optimizer1D::xStep( void )
{
  return cXStep;
}


void Optimizer1D::xStep( double newXStep )
{
  cXStep = newXStep;
}


double Optimizer1D::tolerance( void )
{
  return cTolerance;
}

void Optimizer1D::tolerance( double newTolerance )
{
  cTolerance = newTolerance;
}

unsigned int Optimizer1D::maxIterations( void )
{
    return cMaxIterations;
}

void Optimizer1D::maxIterations( unsigned int newMaxIterations )
{
  cMaxIterations = newMaxIterations;
}


bool Optimizer1D::searchForMin( void )
{
  return cSearchForMin;
}


void Optimizer1D::searchForMin( bool newSearchForMin )
{
  cSearchForMin = newSearchForMin;
}

bool Optimizer1D::extreme( double * x, double * xVal )
{
  if( !cDefined )
    {
    return false;
    }

  return cExtreme( x, xVal );
}


} // end namespace itk
