/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef tubeSplineApproximation1D_h
#define tubeSplineApproximation1D_h

#include "tubeOptimizer1D.h"
#include "tubeSpline1D.h"
#include "tubeUserFunction.h"

#include <vnl/vnl_matrix_fixed.h>
#include <vnl/vnl_vector.h>

namespace tube
{

class SplineApproximation1D : public Spline1D
{
public:
  using Self = SplineApproximation1D;
  using Superclass = Spline1D;
  using Pointer = Self *;
  using ConstPointer = const Self *;

  using MatrixType = vnl_matrix_fixed<double, 4, 4>;
  using VectorType = Superclass::VectorType;
  using ValueFunctionType = Superclass::ValueFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro(SplineApproximation1D);

  /** Constructor. */
  SplineApproximation1D(void);

  /** Constructor. */
  SplineApproximation1D(ValueFunctionType::Pointer funcVal, Optimizer1D::Pointer optimizer1D);

  /** Destructor. */
  ~SplineApproximation1D(void) override;

  double
  DataValue(const VectorType & y, double x) override;

  double
  DataValueD(const VectorType & y, double x) override;

  double
  DataValueD2(const VectorType & y, double x) override;

  double
  DataValueJet(const VectorType & y, double x, double * d, double * d2) override;

protected:
  /** Print out information about this object. */
  void
  PrintSelf(std::ostream & os, Indent indent) const override;

  double     m_SplineApproximation1DMatrixConst;
  MatrixType m_SplineApproximation1DMatrix;

private:
  // Copy constructor not implemented.
  SplineApproximation1D(const Self & self);

  // Copy assignment operator not implemented.
  void
  operator=(const Self & self);

}; // End class SplineApproximation1D

} // End namespace tube

#endif // End !defined( __tubeSplineApproximation1D_h )
