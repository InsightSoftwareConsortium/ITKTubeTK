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

#ifndef __tubeBrentOptimizer1D_h
#define __tubeBrentOptimizer1D_h

#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

namespace tube
{

class BrentOptimizer1D : public Optimizer1D
{

public:

  BrentOptimizer1D( void );
  BrentOptimizer1D( UserFunction< double, double > * newFuncVal,
              UserFunction< double, double > * newFuncDeriv );

  ~BrentOptimizer1D( void );

  void use( UserFunction< double, double > * newFuncVal,
            UserFunction< double, double > * newFuncDeriv );

  double   smallDouble( void );
  void     smallDouble( double newSmall );


protected:

  void m_Move( double & a, double & b, double & c,
              double d, double e, double f );

  bool m_Extreme( double * x, double * xVal );

  double m_Small;

}; // End class BrentOptimizer1D

} // End namespace tube

#endif // End !defined(__tubeBrentOptimizer1D_h)
