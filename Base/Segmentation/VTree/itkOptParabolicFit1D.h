/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptParabolicFit1D.h,v $
  Language:  C++
  Date:      $Date: 2004/10/25 14:23:31 $
  Version:   $Revision: 1.4 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptParabolicFit1D_h
#define __itkOptParabolicFit1D_h

#include "itkOptimizer1D.h"
#include "UserFunc.h"

namespace itk
{

class OptParabolicFit1D : public Optimizer1D 
{

public:

  OptParabolicFit1D();
  OptParabolicFit1D(UserFunc<double, double> *newFuncVal);
  ~OptParabolicFit1D();
  
  void use(UserFunc<double, double> *newFuncVal, UserFunc<double,double>* deriv = NULL);

protected:

  double cCenter(double x1, double y1, 
                 double x2, double y2, 
                 double x3, double y3);
  bool cExtreme(double * x, double * xVal);


};

}; // end namespace itk

#endif

