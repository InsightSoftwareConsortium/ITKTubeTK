/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptimizer1D.h,v $
  Language:  C++
  Date:      $Date: 2005/09/04 16:48:56 $
  Version:   $Revision: 1.5 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptimizer1D_h
#define __itkOptimizer1D_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "UserFunc.h"

/** Solve for local extremes of 1D functions
 *  Must be derived to specify specific optimization method (e.g., OptBrent1D)
 *  \author Stephen R. Aylward
 *  \rewritten Julien Jomier
 *  \date 11/22/99
 *  \todo Transform this to ITK optimizer
 */
namespace itk
{

class Optimizer1D
{
protected :
  bool     cDefined;
  double   cXMin;
  double   cXMax;
  double   cXStep;         
  bool     cSearchForMin;
  double   cTolerance;
  unsigned int     cMaxIterations;
  UserFunc<double, double> *cFuncVal;
  UserFunc<double, double> *cFuncDeriv;
  virtual  bool cExtreme(double *x, double *xVal);

public :
       
  /** Null constructor - insufficient to define class - use "use" function */
  Optimizer1D();
  /** Constructor
   *   \param newFuncVal User derivation of UserFunc to define function to be optimized
   *   \param newFuncDeriv User derivation of UserFunc to define derivative of function to be optimized */
  Optimizer1D(UserFunc<double, double> *newFuncVal, UserFunc<double, double> *newFuncDeriv);
  /**  Destructor */
  virtual ~Optimizer1D();  
  /** Specify new functions to be optimized
   *   \param newFuncVal User derivation of UserFunc to define function to be optimized
   *   \param newFuncDeriv User derivation of UserFunc to define derivative of function to be optimized */
  void     use(UserFunc<double, double> *newFuncVal, UserFunc<double, double> *newFuncDeriv);
  double   xMin();
  void     xMin(double newXMin);
  double   xMax();
  void     xMax(double newXMax);
  double   xStep();
  void     xStep(double newXStep);
  double   tolerance();
  void     tolerance(double newTolerance);
  unsigned int     maxIterations(void);
  void     maxIterations(unsigned int newMaxIterations); 
  bool     searchForMin();
  void     searchForMin(bool newSearchForMin);
  bool     extreme(double *x, double *xVal);
};


}; // end namespace itk

#endif

