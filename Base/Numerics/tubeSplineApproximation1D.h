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

#ifndef __tubeSplineApproximation1D_h
#define __tubeSplineApproximation1D_h

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>

#include "tubeOptimizer1D.h"
#include "tubeSpline1D.h"
#include "tubeUserFunc.h"

namespace tube
{

class SplineApproximation1D : public Spline1D
{
public:

  typedef vnl_vector<double> VectorType;

  SplineApproximation1D( void );
  SplineApproximation1D( UserFunc<int, double> *newFunval,
    Optimizer1D * newOpt1D );

  virtual ~SplineApproximation1D( void );

  double  dataValue( const VectorType & y, double x );

  double  dataValueD( const VectorType & y, double x );

  double  dataValueD2( const VectorType & y, double x );

  double  dataValueJet( const VectorType & y,
    double x, double *d, double *d2 );

protected:

  typedef vnl_matrix_fixed<double, 4, 4> MatrixType;

  double      m_SplineApproximation1DMatrixConst;

  MatrixType  m_SplineApproximation1DMatrix;

}; // End class SplineApproximation1D

} // End namespace tube

#endif // End !defined(__tubeSplineApproximation1D_h)
