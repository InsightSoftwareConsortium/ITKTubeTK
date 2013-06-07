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

#include "tubeOptimizer1D.h"

#include <iostream>

#include "tubeMacro.h"

namespace tube
{

bool
Optimizer1D
::m_Extreme( double * tubeNotUsed(x), double * tubeNotUsed(xVal) )
{
  return false;
}


Optimizer1D
::Optimizer1D( void )
{
  m_SearchForMin = true;
  m_Tolerance = 0.0001;
  m_MaxIterations = 300;
  m_XMin = 0;
  m_XMax = 1;
  m_XStep = 0.01;
  m_Defined = false;
  m_FuncVal = 0;
  m_FuncDeriv = 0;
}


Optimizer1D::Optimizer1D( UserFunction< double, double > * newFuncVal,
  UserFunction< double, double > * newFuncDeriv )
{
  m_SearchForMin = true;
  m_Tolerance = 0.0001;
  m_MaxIterations = 300;
  m_XMin = 0;
  m_XMax = 1;
  m_XStep = 0.01;
  m_FuncVal = newFuncVal;
  m_FuncDeriv = newFuncDeriv;
  m_Defined = true;
}


Optimizer1D::~Optimizer1D( void )
{
}

void Optimizer1D::use( UserFunction< double, double > * newFuncVal,
  UserFunction< double, double > * newFuncDeriv )
{
  m_FuncVal = newFuncVal;
  m_FuncDeriv = newFuncDeriv;
  m_Defined = true;
}


double Optimizer1D::xMin( void )
{
  return m_XMin;
}


void Optimizer1D::xMin( double newXMin )
{
  m_XMin = newXMin;
}

double Optimizer1D::xMax( void )
{
  return m_XMax;
}


void Optimizer1D::xMax( double newXMax )
{
  m_XMax = newXMax;
}

double Optimizer1D::xStep( void )
{
  return m_XStep;
}


void Optimizer1D::xStep( double newXStep )
{
  m_XStep = newXStep;
}


double Optimizer1D::tolerance( void )
{
  return m_Tolerance;
}

void Optimizer1D::tolerance( double newTolerance )
{
  m_Tolerance = newTolerance;
}

unsigned int Optimizer1D::maxIterations( void )
{
  return m_MaxIterations;
}

void Optimizer1D::maxIterations( unsigned int newMaxIterations )
{
  m_MaxIterations = newMaxIterations;
}


bool Optimizer1D::searchForMin( void )
{
  return m_SearchForMin;
}


void Optimizer1D::searchForMin( bool newSearchForMin )
{
  m_SearchForMin = newSearchForMin;
}

bool Optimizer1D::extreme( double * x, double * xVal )
{
  if( !m_Defined )
    {
    return false;
    }

  return m_Extreme( x, xVal );
}

void Optimizer1D::PrintSelf( std::ostream & os ) const
{
  if( m_Defined )
    {
    os << "m_Defined = True" << std::endl;
    }
  else
    {
    os << "m_Defined = False" << std::endl;
    }
  os << "m_XMin = " << m_XMin << std::endl;
  os << "m_XMax = " << m_XMax << std::endl;
  os << "m_XStep = " << m_XStep << std::endl;
  if( m_SearchForMin )
    {
    os << "m_SearchForMin = True" << std::endl;
    }
  else
    {
    os << "m_SearchForMin = False" << std::endl;
    }
  os << "m_Tolerance = " << m_Tolerance << std::endl;
  os << "m_MaxIterations = " << m_MaxIterations << std::endl;
  os << "m_FuncVal = " << m_FuncVal << std::endl;
  os << "m_FuncDeriv = " << m_FuncDeriv << std::endl;
}

} // End namespace tube
