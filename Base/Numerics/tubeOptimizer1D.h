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

#ifndef __tubeOptimizer1D_h
#define __tubeOptimizer1D_h

#include "tubeUserFunction.h"

#include <ostream>

namespace tube
{

/** Solve for local extremes of 1D functions
 *  Must be derived to specify specific optimization method (e.g., BrentOptimizer1D)
 *  \class Optimizer1D
 *  \author Stephen R. Aylward
 *  \rewritten Julien Jomier
 *  \rewritten Stephen R. Aylward
 *  \date 11/22/99
 *  \todo Transform this to ITK optimizer */
class Optimizer1D
{
public:

  /**
   * Null constructor - insufficient to define class; use "use" function */
  Optimizer1D( void );

  /** Constructor
   * \param newFuncVal User derivation of UserFunction to define
   * function to be optimized
   * \param newFuncDeriv User derivation of UserFunction to define
   * derivative of function to be optimized */
  Optimizer1D( UserFunction< double, double > * newFuncVal,
    UserFunction< double, double > * newFuncDeriv );

  /**  Destructor */
  virtual ~Optimizer1D( void );

  /** Specify new functions to be optimized
  * \param newFuncVal User derivation of UserFunction to define
  * function to be optimized
  * \param newFuncDeriv User derivation of UserFunction to define
  * derivative of function to be optimized */
  void     use( UserFunction< double, double > * newFuncVal,
    UserFunction< double, double > * newFuncDeriv );

  double   xMin( void );
  void     xMin( double newXMin );

  double   xMax( void );
  void     xMax( double newXMax );

  double   xStep( void );
  void     xStep( double newXStep );

  double   tolerance( void );
  void     tolerance( double newTolerance );

  unsigned int     maxIterations( void );
  void             maxIterations( unsigned int newMaxIterations );

  bool     searchForMin( void );
  void     searchForMin( bool newSearchForMin );

  bool     extreme( double * x, double * xVal );

  void     PrintSelf( std::ostream & os ) const;

protected:

  bool         m_Defined;
  double       m_XMin;
  double       m_XMax;
  double       m_XStep;
  bool         m_SearchForMin;
  double       m_Tolerance;
  unsigned int m_MaxIterations;

  UserFunction< double, double > * m_FuncVal;
  UserFunction< double, double > * m_FuncDeriv;

  virtual  bool m_Extreme( double * x, double * xVal );

}; // End class Optimizer1D

} // End namespace tube

#endif // End !defined(__tubeOptimizer1D_h)
