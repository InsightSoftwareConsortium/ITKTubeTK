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

#include "tubeOptimizer1D.h"

namespace tube
{

bool
Optimizer1D
::m_Extreme( double * tubeNotUsed( x ), double * tubeNotUsed( xval ) )
{
  return false;
}


Optimizer1D
::Optimizer1D( void )
  : m_FuncVal( NULL ), m_FuncDeriv( NULL )
{
  m_SearchForMin = true;
  m_Tolerance = 0.0001;
  m_MaxIterations = 300;
  m_XMin = 0;
  m_XMax = 1;
  m_XStep = 0.01;
  m_Defined = false;
}


Optimizer1D
::Optimizer1D( ValueFunctionType::Pointer funcVal,
  DerivativeFunctionType::Pointer funcDeriv )
{
  m_SearchForMin = true;
  m_Tolerance = 0.0001;
  m_MaxIterations = 300;
  m_XMin = 0;
  m_XMax = 1;
  m_XStep = 0.01;
  m_FuncVal = funcVal;
  m_FuncDeriv = funcDeriv;
  m_Defined = true;
}


Optimizer1D
::~Optimizer1D( void )
{
}


void
Optimizer1D
::Use( ValueFunctionType::Pointer funcVal,
  DerivativeFunctionType::Pointer funcDeriv )
{
  m_FuncVal = funcVal;
  m_FuncDeriv = funcDeriv;
  m_Defined = true;
}


bool
Optimizer1D
::Extreme( double * x, double * xVal )
{
  if( !m_Defined )
    {
    return false;
    }

  return m_Extreme( x, xVal );
}


void
Optimizer1D
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Defined:       " << m_Defined << std::endl;
  os << indent << "XMin:          " << m_XMin << std::endl;
  os << indent << "XMax:          " << m_XMax << std::endl;
  os << indent << "XStep:         " << m_XStep << std::endl;
  os << indent << "SearchForMin:  " << m_SearchForMin << std::endl;
  os << indent << "Tolerance:     " << m_Tolerance << std::endl;
  os << indent << "MaxIterations: " << m_MaxIterations << std::endl;
  os << indent << "FuncVal:       " << m_FuncVal << std::endl;
  os << indent << "FuncDeriv:     " << m_FuncDeriv << std::endl;
}

} // End namespace tube
