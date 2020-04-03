/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
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

#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

#include <vnl/vnl_vector.h>

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
class Spline1D : public Object
{
public:

  typedef Spline1D                        Self;
  typedef Object                          Superclass;
  typedef Self *                          Pointer;
  typedef const Self *                    ConstPointer;

  typedef vnl_vector< double >            VectorType;
  typedef UserFunction< int, double >     ValueFunctionType;
  typedef UserFunction< double, double >  OptimizerValueFunctionType;
  typedef UserFunction< double, double >  OptimizerDerivativeFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( Spline1D );

  /** Destructor. */
  Spline1D( void );

  /** Constructor.
   * \param funcVal an instance of a derivation of the UserFunction class
   *        used to specify values at control points
   *        ( i.e., integer values )
   * \param optimizer1D a ( possibly NULL constructed ) instance of a
   *        derivation of the Optimizer1D class ( e.g., BrentOptimizer1D ).
   *        Use to find local maximums and minimums.
   * \warning Must set xMin and xMax
   */
  Spline1D( ValueFunctionType::Pointer funcVal,
    Optimizer1D::Pointer optimizer1D );

  /** Destructor. */
  virtual ~Spline1D( void );

  /** Returns the characteristics of spline evaluations near data bounds
   * ( xMin and xMax ). If true, values beyond edges ( xMin and xMax ) are set
   * to zero. If false, values beyond edges are faded to 0 as a function of
   * squared distance from edge.
   */
  tubeGetMacro( Clip, bool );

  /** User specification of characteristics */
  tubeSetMacro( Clip, bool );

  tubeBooleanMacro( Clip );

  /** Returns the control points' ( integer value locations ) lower bound */
  tubeGetMacro( XMin, int );

  /** User specification of lower bound */
  virtual void SetXMin( int xMin );

  /** Sets control points' ( integer value locations ) upper bound */
  tubeGetMacro( XMax, int );

  /** User Specification of upper bound */
  virtual void SetXMax( int xMax );

  /** Tracks the validity of internally maintained intermediate
   * calculations and data.  Returns true if a new spline instance has been
   * created ( e.g., use has been called )
   */
  tubeGetMacro( NewData, bool );

  /** User sets to true to force recalculation of internal data
   *  For example, use to flag that UserFunction has changed externally
   */
  tubeSetMacro( NewData, bool );

  tubeBooleanMacro( NewData );

  /** Returns spline interpolated curvature at x
   *  Calculates the values at control ( integer ) points by calling the
   *  UserFunction and returns the interpolated curvature value between those
   *  points. Type of interpolation is dependent on which spline derivation
   *  is used ( e.g., SplApprox1D ).
   *  Intermediate calculations and control point evaluations are stored to
   *  speed subsequent calls
   */
  double Curv( double x );

  /** Value virtual function, defined by spline types, e.g., SplApprox1D
   *  \param y provides the data local to x; for SplApprox1D, must have
   *         quantity()>=4 where y( 0 ) is the value at x=-1, y( 1 ) is the
   *         value at x=0, y( 2 ) for x=1, and y( 3 ) for x=2
   *  \param x must be between 0 and 1
   *  \warning This function is not normally called directly by the class
   *           user.  The value of y is normally set and x is remapped
   *           internally when the function value is called
   */
  virtual double DataValue( const VectorType & y, double x ) = 0;

  /** First derivative virtual function, defined by spline types
   * ( e.g., SplApprox1D )
   * \param y provides the data local to x; for SplApprox1D, must have
   *        quantity()>=4 where y( 0 ) is the value at x=-1, y( 1 ) is the
   *        value at x=0, y( 2 ) for x=1, and y( 3 ) for x=2
   * \param x must be between 0 and 1
   * \warning This function is not normally called directly by the class
   *          user,   The value of y is normally set and x is remapped
   *          internally when the function ValueD is called
   */
  virtual double DataValueD( const VectorType & y, double x ) = 0;

  /** Second derivative virtual function, defined by specific spline types
   * ( e.g., SplApprox1D )
   * \param y provides the data local to x; for SplApprox1D, must have
   *        quantity()>=4 where y( 0 ) is the value at x=-1, y( 1 ) is the
   *        value at x=0, y( 2 ) for x=1, and y( 3 ) for x=2
   * \param x must be between 0 and 1
   * \warning This function is not normally called directly by the class
   *          user,   The value of y is normally set and x is remapped
   *          internally when the function ValueD2 is called.
   */
  virtual double DataValueD2( const VectorType & y, double x ) = 0;

  /** Jet virtual function - returns value and sets first and second
   * derivative, defined by specific spline types ( e.g., SplApprox1D )
   * \param y provides the data local to x; for SplApprox1D, for
   *        SplApprox1D, must have quantity()>=4 where y( 0 ) is
   *        the value at x=-1, y( 1 ) is the value at x=0, y( 2 ) for x=1,
   *        and y( 3 ) for x=2
   * \param x must be between 0 and 1
   * \param d returns first derivative value at x
   * \param d2 returns second derivative value at x
   * \warning This function is not normally called directly by the class
   *          user.  The value of y is normally set and x is remapped
   *          internally when the function valueJet is called
   */
  virtual double DataValueJet( const VectorType & y,
    double x,
    double * d,
    double * d2 ) = 0;

  /** Calculates the local extreme using the supplied instance of a
   *  derivation of Optimizer1D.  Function returns true on successful local
   *  extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool Extreme( double * extX, double * extVal );

  /** Supply a new spline definition
   * \param funcVal an instance of a derivation of the UserFunction class
   *        used to specify values at control points
   *        ( i.e., integer values )
   * \param optimizer1D a ( possibly NULL constructed ) instance of a
   *        derivation of the Optimizer1D class ( e.g., BrentOptimizer1D ).
   *        Use to find local maximums and minimums.
   * \warning Must set xMin and xMax
   */
  void Use( ValueFunctionType::Pointer funcVal,
    Optimizer1D::Pointer optimizer1D );

  /** Returns spline interpolated value at x
   * Calculates the values at control ( integer ) points by calling the
   * UserFunction and returns the interpolated value between those points.
   * Type of interpolation is dependent on which spline derivation is
   * used ( e.g., SplApprox1D ).   Intermediate calculations and control
   * point evaluations are stored to speed subsequent calls
   */
  double Value( double x );

  /** Returns spline interpolated first derivative at x
   *  Calculates the values at control ( integer ) points by calling the
   *  UserFunction and returns the interpolated first derivative between those
   *  points.   Type of interpolation is dependent on which spline
   *  derivation is used ( e.g., SplApprox1D ).   Intermediate calculations
   *  and control point evaluations are stored to speed subsequent calls
   */
  double ValueD( double x );

  /** Returns spline interpolated second derivative at x
   *  Calculates the values at control ( integer ) points by calling the
   *  UserFunction and returns the interpolated second derivative between those
   *  points.   Type of interpolation is dependent on which spline
   *  derivation is used ( e.g., SplApprox1D ).   Intermediate calculations
   *  and control point evaluations are stored to speed subsequent calls
   */
  double ValueD2( double x );

  /** Returns spline interpolated derivative jet at x
   *  Calculates the values at control ( integer ) points by calling the
   *  UserFunction and returns the interpolated value, first, and second
   *  derivatives between those points.  Type of interpolation is dependent
   *  on which spline derivation is used ( e.g., SplApprox1D ).  Intermediate
   *  calculations and control point evaluations are stored to speed
   *  subsequent calls
   */
  double ValueJet( double x, double * d, double * d2 );

protected:

  /** Print out information about this object. */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void m_GetData( double x );

  bool                                      m_Defined;
  ValueFunctionType::Pointer                m_FuncVal;
  bool                                      m_Clip;
  int                                       m_XMin;
  int                                       m_XMax;
  bool                                      m_NewData;
  VectorType                                m_Data;
  OptimizerValueFunctionType::Pointer       m_Optimizer1DVal;
  OptimizerDerivativeFunctionType::Pointer  m_Optimizer1DDeriv;
  Optimizer1D::Pointer                      m_Optimizer1D;

private:

  // Copy constructor not implemented.
  Spline1D( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class Spline1D

} // End namespace itk

#endif // End !defined( __tubeSpline1D_h )
