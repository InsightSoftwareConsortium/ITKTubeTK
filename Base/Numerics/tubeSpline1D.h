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
#ifndef __tubeSpline1D_h
#define __tubeSpline1D_h

#include "tubeUserFunc.h"
#include "tubeOptimizer1D.h"

namespace tube
{

/** Example of how to use Spline1D.   Program is run to verify methods
*  \example TestSpline1D/testSpline1D.cpp
*/

/** 1D Spline abstract base class
*
*  Provides a consistent interface to setting the common
*      parameters of spline functions.   Also, hides the
*      archiving of previous calculations and control point
*      evaluations to speed subsequent spline evaluations.
*      Example derivation in SplApprox1D class
*  \class Spline1D
*  \author Stephen R. Aylward
*  \date   11/19/99 */
class  Spline1D
{
public:

  typedef vnl_vector<double> VectorType;

  Spline1D();

  /** Construct a viable class instance
   * \param newFuncVal an instance of a derivation of the UserFunc class
   *        used to specify values at control points (i.e., integer values)
   * \param newOpt1D a (possibly NULL constructed) instance of a derivation
   *        of the Optimizer1D class (e.g., OptBrent1D). Use to find local
   *        maxs and mins.
   * \warning Must set xMin and xMax
   */
  Spline1D( UserFunc<int, double> * newFuncVal, Optimizer1D * newOpt1D);

  //! Virtual destructor
  virtual ~Spline1D();

  /** Supply a new spline definition
   * \param newFuncVal an instance of a derivation of the UserFunc class
   *        used to specify values at control points (i.e., integer values)
   * \param newOpt1D a (possibly NULL constructed) instance of a derivation
   *        of the Optimizer1D class (e.g., OptBrent1D). Use to find local
   *        maxs and mins.
   * \warning Must set xMin and xMax
   */
  void    use( UserFunc<int, double> * newFuncVal, Optimizer1D * newOpt1D);

  /** Returns the characteristics of spline evaluations near data bounds
   * (xMin and xMax). If true, values beyond edges (xMin and xMax) are set
   * to zero. If false, values beyond edges are faded to 0 as a function of
   * distance from edge, squeared.
   */
  bool    clipEdge();

  /** User specification of characteristics */
  void    clipEdge(bool newClip);

  /** Returns the control points' (integer value locations) lower bound */
  int     xMin();
  /** User specification of lower bound */
  void    xMin(int newXMin);

  /** Sets control points' (integer value locations) upper bound */
  int     xMax();
  /** User Specification of upper bound */
  void    xMax(int newXMax);

  /** Tracks the validity of internally maintained intermediate
   * calculations and data.  Returns true if a new spline instance has been
   * created (e.g., use has been called)
   */
  bool    newData();

  /** User sets to true to force recalcuation of internal data
   *  For example, use to flag that UserFunc has changed externally
   */
  void    newData(bool newNewData);

  /** Value virtual function, defined by spline types, e.g., SplApprox1D
   *  \param y provides the data local to x; for SplApprox1D, must have
   *         quantity()>=4 where y(0) is the value at x=-1, y(1) is the
   *         value at x=0, y(2) for x=1, and y(3) for x=2
   *  \param x must be between 0 and 1
   *  \warning This function is not normally called directly by the class
   *           user.  The value of y is normally set and x is remapped
   *           internally when the function value is called
   */
  virtual double  dataValue(const VectorType & y, double x) = 0;

  /** First derivative virtual function, defined by spline types
   * (e.g., SplApprox1D)
   * \param y provides the data local to x; for SplApprox1D, must have
   *        quantity()>=4 where y(0) is the value at x=-1, y(1) is the
   *        value at x=0, y(2) for x=1, and y(3) for x=2
   * \param x must be between 0 and 1
   * \warning This function is not normally called directly by the class
   *          user,   The value of y is normally set and x is remapped
   *          internally when the function valueD is called
   */
  virtual double  dataValueD(const VectorType & y, double x) = 0;

  /** Second derivative virtual function, defined by specific spline types
   * (e.g., SplApprox1D)
   * \param y provides the data local to x; for SplApprox1D, must have
   *        quantity()>=4 where y(0) is the value at x=-1, y(1) is the
   *        value at x=0, y(2) for x=1, and y(3) for x=2
   * \param x must be between 0 and 1
   * \warning This function is not normally called directly by the class
   *          user,   The value of y is normally set and x is remapped
   *          internally when the function valueD2 is called.
   */
  virtual double  dataValueD2(const VectorType & y, double x) = 0;

  /** Jet virtual function - returns value and sets first and second
   * derivative, defined by specific spline types (e.g., SplApprox1D)
   * \param y provides the data local to x; for SplApprox1D, for
   *        SplApprox1D, must have quantity()>=4 where y(0) is
   *        the value at x=-1, y(1) is the value at x=0, y(2) for x=1,
   *        and y(3) for x=2
   * \param x must be between 0 and 1
   * \param d returns first derivative value at x
   * \param d2 returns second derivative value at x
   * \warning This function is not normally called directly by the class
   *          user.  The value of y is normally set and x is remapped
   *          internally when the function valueJet is called
   */
  virtual double  dataValueJet(const VectorType & y, double x,
    double * d, double * d2) = 0;

  /** Returns spline interpolated value at x
   * Calculates the values at control (integer) points by calling the
   * UserFunc and returns the interpolated value between those points.
   * Type of interpolation is dependent on which spline derivation is
   * used (e.g., SplApprox1D).   Intermediate calculations and control
   * point evaluations are stored to speed subsequent calls
   */
  double  value(double x);

  /** Returns spline interpolated first derivative at x
   *  Calculates the values at control (integer) points by calling the
   *  UserFunc and returns the interpolated first derivative between those
   *  points.   Type of interpolation is dependent on which spline
   *  derivation is used (e.g., SplApprox1D).   Intermediate calculations
   *  and control point evaluations are stored to speed subsequent calls
   */
  double  valueD(double x);

  /** Returns spline interpolated second derivative at x
   *  Calculates the values at control (integer) points by calling the
   *  UserFunc and returns the interpolated second derivative between those
   *  points.   Type of interpolation is dependent on which spline
   *  derivation is used (e.g., SplApprox1D).   Intermediate calculations
   *  and control point evaluations are stored to speed subsequent calls
   */
  double  valueD2(double x);

  /** Returns spline interpolated derivative jet at x
   *  Calculates the values at control (integer) points by calling the
   *  UserFunc and returns the interpolated value, first, and second
   *  derivatives between those points.  Type of interpolation is dependent
   *  on which spline derivation is used (e.g., SplApprox1D).  Intermediate
   *  calculations and control point evaluations are stored to speed
   *  subsequent calls
   */
  double  valueJet(double x, double * d, double * d2);

  /** Returns spline interpolated curvature at x
   *  Calculates the values at control (integer) points by calling the
   *  UserFunc and returns the interpolated curvature value between those
   *  points. Type of interpolation is dependent on which spline derivation
   *  is used (e.g., SplApprox1D).
   *  Intermediate calculations and control point evaluations are stored to
   *  speed subsequent calls
   */
  double  curv(double x);

  /** Calculates the local extreme using the supplied instance of a
   *  derivation of Optimizer1D.  Function returns true on successful local
   *  extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool    extreme(double *extX, double *extVal);

  void PrintSelf( std::ostream & os ) const;

protected:

  bool                       m_Defined;
  UserFunc<int, double>    * m_FuncVal;
  bool                       m_Clip;
  int                        m_XMin;
  int                        m_XMax;

  bool                       m_NewData;
  VectorType                 m_Data;

  UserFunc<double, double> * m_Opt1DVal;
  UserFunc<double, double> * m_Opt1DDeriv;
  Optimizer1D              * m_Opt1D;

  void                       m_GetData(double x);

};

} // end namespace itk

#endif /* __itkRadiusExtractor_h */
