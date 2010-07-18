/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptGoldenMean1D.h,v $
  Language:  C++
  Date:      $Date: 2005/09/04 16:48:56 $
  Version:   $Revision: 1.5 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptGoldenMean1D_h
#define __itkOptGoldenMean1D_h

#include "itkOptimizer1D.h"
#include "UserFunc.h"

namespace itk
{

class OptGoldenMean1D : public Optimizer1D 
{

public:

  OptGoldenMean1D();
  OptGoldenMean1D(UserFunc<double, double> *newFuncVal);
  ~OptGoldenMean1D();

  void use(UserFunc<double, double> *newFuncVal);

protected:

  bool cExtreme(double * x, double * xVal);


};

}; // end namespace itk

#endif

