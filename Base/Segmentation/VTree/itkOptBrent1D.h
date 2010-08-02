/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

  OptBrent1D( void );
  OptBrent1D( UserFunc< double, double > * newFuncVal,
              UserFunc< double, double > * newFuncDeriv );

  ~OptBrent1D( void );

  void use( UserFunc< double, double > * newFuncVal,
            UserFunc< double, double > * newFuncDeriv );

  double   smallDouble( void );
  void     smallDouble( double newSmall );


protected:

  void cMove( double & a, double & b, double & c,
              double d, double e, double f );

  bool cExtreme( double * x, double * xVal );

  double cSmall;

};

}; // end namespace itk

#endif

