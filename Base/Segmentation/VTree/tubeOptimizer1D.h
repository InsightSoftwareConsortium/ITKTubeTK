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
#ifndef __tubeOptimizer1D_h
#define __tubeOptimizer1D_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "tubeUserFunc.h"

/** Solve for local extremes of 1D functions
 *  Must be derived to specify specific optimization method (e.g., OptBrent1D)
 *  \author Stephen R. Aylward
 *  \rewritten Julien Jomier
 *  \date 11/22/99
 *  \todo Transform this to ITK optimizer
 */
namespace tube
{

class Optimizer1D
{
protected :
  bool         cDefined;
  double       cXMin;
  double       cXMax;
  double       cXStep;
  bool         cSearchForMin;
  double       cTolerance;
  unsigned int cMaxIterations;

  UserFunc< double, double > * cFuncVal;
  UserFunc< double, double > * cFuncDeriv;

  virtual  bool cExtreme( double * x, double * xVal );

public :

  /** Null constructor - insufficient to define class; use "use" function */
  Optimizer1D( void );

  /** Constructor
   *   \param newFuncVal User derivation of UserFunc to define
   *   function to be optimized
   *   \param newFuncDeriv User derivation of UserFunc to define
   *   derivative of function to be optimized */
  Optimizer1D( UserFunc< double, double > * newFuncVal,
    UserFunc< double, double > * newFuncDeriv );

  /**  Destructor */
  virtual ~Optimizer1D();

  /** Specify new functions to be optimized
   *   \param newFuncVal User derivation of UserFunc to define
   *   function to be optimized
   *   \param newFuncDeriv User derivation of UserFunc to define
   *   derivative of function to be optimized */
  void     use( UserFunc< double, double > * newFuncVal,
    UserFunc< double, double > * newFuncDeriv );

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

};


}; // end namespace tube

#endif
