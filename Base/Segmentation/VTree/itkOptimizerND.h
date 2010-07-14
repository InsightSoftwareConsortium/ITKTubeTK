/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptimizerND.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:25 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptimizerND_h
#define __itkOptimizerND_h

#include "UserFunc.h"
#include "itkOptimizer1D.h"
#include <vnl/vnl_vector.h>
#include <vnl/vnl_matrix.h>
#include "itkMatrixMath.h"

namespace itk 
{


class OptimizerND
{

public :

  /** Typedef for the vector type used */
  typedef vnl_vector<double> VectorType;

  /** Typedef for the matrix type used */
  typedef vnl_matrix<double> MatrixType;

  OptimizerND();
  OptimizerND(int newNDims, UserFunc<VectorType*, double> * newFuncValND, UserFunc<VectorType*, VectorType &> * newFuncDerivND, Optimizer1D *newOpt1D);
  virtual ~OptimizerND();

  void use(int newNDims, UserFunc<VectorType*, double> * newFuncValND, UserFunc<VectorType*, VectorType &> * newFuncDerivND, Optimizer1D *newOpt1D);

  VectorType & xMin(void);
  void    xMin(VectorType & newXMinn);

  VectorType & xMax(void);

  void    xMax(VectorType & newXMaxx);
  VectorType & xStep(void);
  void    xStep(VectorType & newXStepp);

  double  tolerance();
  void    tolerance(double newTolerance);
  unsigned int    maxIterations(void);
  void    maxIterations(unsigned int newMaxIterations);
  bool    searchForMin();
  void    searchForMin(bool newSearchForMin);
  double  funcVal(double x);
  double  funcDeriv(double x);
  bool    extreme(VectorType &x, double *xVal);
  bool    extreme(VectorType &x, double *xVal, unsigned int n, MatrixType &dirs);

protected :
       
  bool            cDefined;
  unsigned int            cNDims;
  
  VectorType *cXMin;
  VectorType *cXMax;
  VectorType *cXStep;
  VectorType *cX0;
  VectorType *cX0Dir;
  VectorType *cX0Temp;
  bool            cSearchForMin;
  double          cTolerance;
  unsigned int            cMaxIterations;
  UserFunc<double, double> * cOpt1DVal;
  UserFunc<double, double> * cOpt1DDeriv;
  Optimizer1D *   cOpt1D;
  UserFunc<VectorType*, double> *      cFuncValND;
  UserFunc<VectorType*, VectorType &> *      cFuncDerivND;
  
};



}; // end namespace itk



#endif

