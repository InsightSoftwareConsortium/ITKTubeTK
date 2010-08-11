/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkSplineApproximation1D.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:26 $
  Version:   $Revision: 1.3 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkSplineApproximation1D_h
#define itkSplineApproximation1D_h

#include "itkSpline1D.h"
#include "itkUserFunc.h"
#include "itkOptimizer1D.h"

namespace itk 
{

class SplineApproximation1D : public Spline1D
{

public :

  typedef vnl_vector<double> VectorType;
  SplineApproximation1D();
  SplineApproximation1D(UserFunc<int, double> *newFunval, Optimizer1D * newOpt1D);
  virtual ~SplineApproximation1D();
  double  dataValue(VectorType y, double x);
  double  dataValueD(VectorType y, double x);
  double  dataValueD2(VectorType y, double x);
  double  dataValueJet(VectorType y, double x, double *d, double *d2);

protected :

  float cSplineApproximation1DMatrixConst;
  float cSplineApproximation1DMatrix[4][4];



};

}; // end namespace itk


#endif

