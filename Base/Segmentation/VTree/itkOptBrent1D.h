/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkOptBrent1D.h,v $
  Language:  C++
  Date:      $Date: 2005/09/04 16:48:56 $
  Version:   $Revision: 1.6 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkOptBrent1D_h
#define __itkOptBrent1D_h

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkOptimizer1D.h"
#include "UserFunc.h"

namespace itk
{

class OptBrent1D : public Optimizer1D 
{

public:

  OptBrent1D();
  OptBrent1D(UserFunc<double, double> *newFuncVal);
  ~OptBrent1D();

  void use(UserFunc<double, double> *newFuncVal);

  double   smallDouble();
  void     smallDouble(double newSmall);


protected:

  void cMove(double & a, double & b, double & c, double d, double e, double f);
  bool cExtreme(double * x, double * xVal);

  double cSmall;


};

}; // end namespace itk

#endif

