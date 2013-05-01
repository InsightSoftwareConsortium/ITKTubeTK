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
#ifndef __tubeOptParabolicFit1D_h
#define __tubeOptParabolicFit1D_h

#include "tubeOptimizer1D.h"
#include "tubeUserFunc.h"

namespace tube
{

class OptParabolicFit1D : public Optimizer1D
{

public:

  OptParabolicFit1D( void );
  OptParabolicFit1D( UserFunc< double, double > * newFuncVal );
  ~OptParabolicFit1D( void );

  void use( UserFunc< double, double > * newFuncVal,
            UserFunc< double, double > * deriv = NULL );

protected:

  double m_Center( double x1, double y1,
                  double x2, double y2,
                  double x3, double y3 );
  bool m_Extreme( double * x, double * xVal );


};

} // end namespace tube

#endif
