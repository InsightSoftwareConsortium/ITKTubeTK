/*=========================================================================

Library:   TubeTK/VTree

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
  SplineApproximation1D(UserFunc<int, double> *newFunval,
    Optimizer1D * newOpt1D);

  virtual ~SplineApproximation1D();

  double  dataValue(const VectorType & y, double x);

  double  dataValueD(const VectorType & y, double x);

  double  dataValueD2(const VectorType & y, double x);

  double  dataValueJet(const VectorType & y,
    double x, double *d, double *d2);

protected :

  float cSplineApproximation1DMatrixConst;
  float cSplineApproximation1DMatrix[4][4];

};

}; // end namespace itk

#endif

