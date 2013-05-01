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
#ifndef __tubeSplineND_h
#define __tubeSplineND_h

#include "itkImage.h"
#include "itkVectorContainer.h"

#include "tubeSpline1D.h"
#include "tubeOptimizer1D.h"
#include "tubeOptimizerND.h"

namespace tube
{

/*! Multidimensional spline class
 *  Evaluates user supplied function to determine values at integers
 *  and will then interpolate (based on derivation used) non-integer
 *  values.  Includes methods for determining value, first derivative,
 *  hessian, and local extrema.   Relies on a derivation of Spline1D to
 *  specify how interpolation is performed.
 *  \author Stephen R. Aylward
 *  \date 11/21/99
 */

class SplineND
{
public:

  /** Typedef for the vector type used */
  typedef vnl_vector< double > VectorType;

  /** Typedef for the vector type used */
  typedef vnl_vector< int > IntVectorType;

  /** Typedef for the matrix type used */
  typedef vnl_matrix< double > MatrixType;

  /** Typedef for the multidimensional structure (ie an image) */
  typedef itk::Image< double, 4 >   ImageType;

  /** Default constructor: insufficient for using class */
  SplineND( void );

  /** Constructor produces usable instance of this class
  *  \param newNDims dimensionality of the space being interpolated
  *  \param newFuncVal function evaluated at integer points to determine
  *         control point values
  *  \param newSpline1D derivation of Spline1D (e.g., SplApprox1D) used
  *         to marginally (i.e., per dimension) interpolate values
  *  \param newOptND an instance (can be NULL constructed) of OptimizerND
  *         that is used to find local extrema
  *  \warning xMin and xMax must be set!
  */
  SplineND(unsigned int newNDims,
    UserFunc< IntVectorType, double > * newFuncVal,
    Spline1D * newSpline1D,
    Optimizer1D * newOptND);

  /** Destructor
  *  This class is not usually derived
  */
  virtual ~SplineND( void );

  /** Specify a new spline function
  *  \param newNDims dimensionality of the space being interpolated
  *  \param newFuncVal function evaluated at integer points to determine
  *         control point values
  *  \param newSpline1D derivation of Spline1D (e.g., SplApprox1D) used
  *         to marginally (i.e., per dimension) interpolate values
  *  \param newOptND an instance (can be NULL constructed) of OptimizerND
  *         that is used to find local extrema
  *  \warning xMin and xMax must be set!
  */
  void use(unsigned int newNDims,
    UserFunc<IntVectorType, double > * newFuncVal,
    Spline1D * newSpline1D,
    Optimizer1D * newOptND);

  /** Returns the characteristics of spline evaluations near data bounds
  *   (xMin and xMax). If true, values beyond edges (xMin and xMax) are
  *   set to zero.  If false, values beyond edges are faded to 0 as a
  *   function of distance from edge, squeared.
  */
  bool clipEdge( void );

  /** Sets the characteristics of spline evaluations near data bounds
  *   (xMin and xMax). If true, values beyond edges (xMin and xMax) are
  *   set to zero.  If false, values beyond edges are faded to 0 as a
  *   function of distance from edge, squeared.
  */
  void clipEdge(bool newClip);

  /** Returns the number of dimensions of the problem's domain.  */
  int nDims( void ) { return m_NDims; }

  /** User specification of lower bound */
  void xMin(const IntVectorType & newXMin);

  /** Returns the control points' (integer value locations) lower bound */
  const IntVectorType & xMin( void );

  /** User Specification of upper bound */
  void xMax(const IntVectorType & newXMax);

  /** Sets control points' (integer value locations) upper bound */
  const IntVectorType & xMax( void );

  /** Tracks the validity of internally maintained intermediate
   * calculations and data. Returns true if a new spline instance has been
   * created (e.g., use has been called)
   */
  bool newData( void );

  /** User sets to true to force recalcuation of internal data.
   * For example, use to flag that UserFunc has changed externally
   */
  void newData(bool newNewData);

  /** Returns spline interpolated value at x.
   * Calculates the values at control (integer) points by calling the
   * UserFunc and returns the interpolated value between those points.
   * Type of interpolation is dependent on which spline derivation is used
   * (e.g., SplApprox1D).   Intermediate calculations and control point
   * evaluations are stored to speed subsequent calls.
   */
  const double & value(const VectorType & x);

  /** Returns spline interpolated first derivative at x projected onto dx.
   *  Calculates the values at control (integer) points by calling the
   *  UserFunc and returns the interpolated first derivative between those
   *  points.  Type of interpolation is dependent on which spline derivation
   *  is used (e.g., SplApprox1D).   Intermediate calculations and control
   *  point evaluations are stored to speed subsequent calls.
   */
  double valueD(const VectorType & x, IntVectorType & dx);

  /** Returns spline interpolated first derivative at x.
   * Calculates the values at control (integer) points by calling the
   * UserFunc and returns the interpolated first derivative between those
   * points.  Type of interpolation is dependent on which spline derivation
   * is used (e.g., SplApprox1D).  Intermediate calculations and control
   * point evaluations are stored to speed subsequent calls.
   */
  VectorType & valueD(const VectorType & x);

  /** Returns spline interpolated Hessian at x.
   * Calculates the values at control (integer) points by calling the
   * UserFunc and returns the interpolated Hessian between those points.
   * Type of interpolation is dependent on which spline derivation is used
   * (e.g., SplApprox1D).  Intermediate calculations and control point
   * evaluations are stored to speed subsequent calls
   */
  MatrixType & hessian(const VectorType & x);

  /** Returns spline interpolated derivative jet (value, 1st deriv,
   * Hessian) at x. Calculates the values at control (integer) points by
   * calling the UserFunc and returns the interpolated derivative jet
   * between those points.   Type of interpolation is dependent on which
   * spline derivation is used (e.g., SplApprox1D).   Intermediate
   * calculations and control point evaluations are stored to speed
   * subsequent calls.
   */
  double valueJet(const VectorType & x, VectorType & d, MatrixType & h);

  /** Returns spline interpolated 1st derivatives and 2nd derivatives at x
   * Calculates the values at control (integer) points by calling the
   * UserFunc and returns the interpolated 1st derivatives and 2nd
   * derivatives between those points.  Type of interpolation is dependent
   * on which spline derivation is used (e.g., SplApprox1D). Intermediate
   * calculations and control point evaluations are stored to speed
   * subsequent calls
   */
  double valueVDD2(const VectorType & x, VectorType & d, VectorType & d2);

  OptimizerND * optimizerND( void );

  /** Calculates the local extreme using the supplied instance of a
   * derivation of OptimizerND.  Function returns true on successful local
   * extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location
   *         of extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool extreme(VectorType & extX, double * extVal);

  /** Calculates the local extreme in the direction dir using the supplied
   * instance of a derivation of OptimizerND. Function returns true on
   * successful local extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   *  \param dir Direction to search for local extreme
   */
  bool extreme(VectorType & extX, double * extVal, VectorType &dir);

  /** Calculates the local extreme in the basis space dirs using the
   * supplied instance of a derivation of OptimizerND. Function returns
   * true on successful local extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location
   *         of extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   *  \param n number of vectors in dirs to use to define the basis space
   *  \param dirs TNT::Vectors that define the basis space to search for
   *         local extreme
   */
  bool extreme(VectorType & extX, double * extVal, unsigned int n,
    MatrixType &dirs);

  /** Calculates the local extreme using an approximation to the conjugate
   * gradient descent method. Function returns true on successful local
   * extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool extremeConjGrad(VectorType & extX, double * extVal);

protected:

  typedef itk::VectorContainer< unsigned int, ImageType::Pointer >
    VectorImageType;

  /** Used to enable/disable cout of intermediate calculations */
  bool           m_Debug;

  unsigned int   m_NDims;
  bool           m_Clip;
  IntVectorType  m_XMin;
  IntVectorType  m_XMax;
  bool           m_NewData;
  IntVectorType  m_Xi;
  double         m_Val;
  VectorType     m_D;
  MatrixType     m_H;

  ImageType::Pointer  m_Data;
  ImageType::Pointer  m_DataWS;
  VectorType          m_Data1D;

  VectorImageType::Pointer m_DataWSX;
  VectorImageType::Pointer m_DataWSXX;

  UserFunc< IntVectorType, double >   * m_FuncVal;

  UserFunc< VectorType, double >      * m_OptNDVal;
  UserFunc< VectorType, VectorType >  * m_OptNDDeriv;

  OptimizerND       * m_OptND;
  Spline1D          * m_Spline1D;

  void  m_GetData(const VectorType &x);

private:

  /** Prevent copying and assignment */
  SplineND(const SplineND &);
  SplineND& operator=(const SplineND &);
};

} //end namespace tube

#endif
