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
#ifndef __tubeOptBrent1D_h
#define __tubeOptBrent1D_h

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#endif

#include "tubeOptimizer1D.h"
#include "tubeUserFunc.h"

namespace tube
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

  void m_Move( double & a, double & b, double & c,
              double d, double e, double f );

  bool m_Extreme( double * x, double * xVal );

  double m_Small;

}; // End class OptBrent1D

} // End namespace tube

#endif // End !defined(__tubeOptBrent1D_h)
