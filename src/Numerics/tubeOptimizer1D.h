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

#ifndef __tubeOptimizer1D_h
#define __tubeOptimizer1D_h

#include "tubeObject.h"
#include "tubeUserFunction.h"

namespace tube
{

/** Solve for local extremes of 1D functions
 *  Must be derived to specify specific optimization method
 *  ( e.g., BrentOptimizer1D )
 *  \class Optimizer1D
 *  \author Stephen R. Aylward
 *  \rewritten Julien Jomier
 *  \rewritten Stephen R. Aylward
 *  \date 11/22/99
 *  \todo Transform this to ITK optimizer */
class Optimizer1D : public Object
{
public:

  typedef Optimizer1D                     Self;
  typedef Object                          Superclass;
  typedef Self *                          Pointer;
  typedef const Self *                    ConstPointer;

  typedef UserFunction< double, double >  ValueFunctionType;
  typedef UserFunction< double, double >  DerivativeFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( Optimizer1D );

  /** Constructor. */
  Optimizer1D( void );

  /** Constructor.
   * \param funcVal User derivation of UserFunction to define
   * function to be optimized
   * \param funcDeriv User derivation of UserFunction to define
   * derivative of function to be optimized */
  Optimizer1D( ValueFunctionType::Pointer funcVal,
    DerivativeFunctionType::Pointer funcDeriv );

  /** Destructor. */
  virtual ~Optimizer1D( void );

  tubeGetMacro( MaxIterations, unsigned int );

  tubeSetMacro( MaxIterations, unsigned int );

  tubeGetMacro( SearchForMin, bool );

  tubeSetMacro( SearchForMin, bool );

  tubeBooleanMacro( SearchForMin );

  tubeGetMacro( Tolerance, double );

  tubeSetMacro( Tolerance, double );

  tubeGetMacro( XMax, double );

  tubeSetMacro( XMax, double );

  tubeGetMacro( XMin, double );

  tubeSetMacro( XMin, double );

  tubeGetMacro( XStep, double );

  tubeSetMacro( XStep, double );

  bool Extreme( double * x, double * xVal );

  /** Specify new functions to be optimized
  * \param funcVal User derivation of UserFunction to define
  * function to be optimized
  * \param funcDeriv User derivation of UserFunction to define
  * derivative of function to be optimized */
  void Use( ValueFunctionType::Pointer funcVal,
    DerivativeFunctionType::Pointer funcDeriv );

protected:

  virtual bool m_Extreme( double * x, double * xVal );

  bool                             m_Defined;
  double                           m_XMin;
  double                           m_XMax;
  double                           m_XStep;
  bool                             m_SearchForMin;
  double                           m_Tolerance;
  unsigned int                     m_MaxIterations;
  ValueFunctionType::Pointer       m_FuncVal;
  DerivativeFunctionType::Pointer  m_FuncDeriv;

protected:

  /** Print out information about this object. */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

private:

  // Copy constructor not implemented.
  Optimizer1D( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class Optimizer1D

} // End namespace tube

#endif // End !defined( __tubeOptimizer1D_h )
