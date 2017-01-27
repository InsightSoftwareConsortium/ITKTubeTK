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

#ifndef __tubeSplineND_h
#define __tubeSplineND_h

#include "tubeOptimizer1D.h"
#include "tubeOptimizerND.h"
#include "tubeSpline1D.h"
#include "tubeUserFunction.h"

#include <itkImage.h>
#include <itkVectorContainer.h>

namespace tube
{

/** Multidimensional spline class
 *  Evaluates user supplied function to determine values at integers
 *  and will then interpolate ( based on derivation used ) non-integer
 *  values.  Includes methods for determining value, first derivative,
 *  Hessian, and local extrema.   Relies on a derivation of Spline1D to
 *  specify how interpolation is performed.
 *  \author Stephen R. Aylward
 *  \date 11/21/99
 */

class SplineND : public Object
{
public:

  typedef SplineND                                Self;
  typedef Object                                  Superclass;
  typedef Self *                                  Pointer;
  typedef const Self *                            ConstPointer;

  typedef itk::Image< double, 4 >                 ImageType;
  typedef vnl_vector< int >                       IntVectorType;
  typedef vnl_matrix< double >                    MatrixType;
  typedef vnl_vector< double >                    VectorType;
  typedef UserFunction< IntVectorType, double >   ValueFunctionType;

  typedef UserFunction< VectorType, double >
    OptimizerValueFunctionType;
  typedef UserFunction< VectorType, VectorType >
    OptimizerDerivativeFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( SplineND );

  /** Constructor. */
  SplineND( void );

  /** Constructor.
   *  \param dimension dimensionality of the space being interpolated
   *  \param funcVal function evaluated at integer points to determine
   *         control point values
   *  \param spline1D derivation of Spline1D ( e.g., SplApprox1D ) used
   *         to marginally ( i.e., per dimension ) interpolate values
   *  \param optimizer1D an instance ( can be NULL constructed ) of OptimizerND
   *         that is used to find local extrema
   *  \warning xMin and xMax must be set!
   */
  SplineND( unsigned int dimension,
    ValueFunctionType::Pointer funcVal,
    Spline1D::Pointer spline1D,
    Optimizer1D::Pointer optimizer1D );

  /** Destructor. */
  virtual ~SplineND( void );

  /** Returns the characteristics of spline evaluations near data bounds
   *   ( xMin and xMax ). If true, values beyond edges ( xMin and xMax ) are
   *   set to zero.  If false, values beyond edges are faded to 0 as a
   *   function of squared distance from edge.
   */
  tubeGetMacro( Clip, bool );

  /** Sets the characteristics of spline evaluations near data bounds
   *   ( xMin and xMax ). If true, values beyond edges ( xMin and xMax ) are
   *   set to zero.  If false, values beyond edges are faded to 0 as a
   *   function of squared distance from edge.
   */
  tubeSetMacro( Clip, bool );

  tubeBooleanMacro( Clip );

  /** Returns the dimension of the problem's domain.  */
  tubeGetMacro( Dimension, unsigned int );

  OptimizerND::Pointer GetOptimizerND( void );

  /** User specification of lower bound */
  virtual void SetXMin( IntVectorType xMin );

  /** Returns the control points' ( integer value locations ) lower bound */
  tubeGetMacro( XMin, IntVectorType );

  /** User Specification of upper bound */
  virtual void SetXMax( IntVectorType xMax );

  /** Sets control points' ( integer value locations ) upper bound */
  tubeGetMacro( XMax, IntVectorType );

  /** Tracks the validity of internally maintained intermediate
   * calculations and data. Returns true if a new spline instance has been
   * created ( e.g., use has been called )
   */
  tubeGetMacro( NewData, bool );

  /** User sets to true to force recalculation of internal data.
   * For example, use to flag that UserFunction has changed externally
   */
  tubeSetMacro( NewData, bool );

  tubeBooleanMacro( NewData );

  /** Calculates the local extreme using the supplied instance of a
   * derivation of OptimizerND.  Function returns true on successful local
   * extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location
   *         of extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool Extreme( VectorType & extX, double * extVal );

  /** Calculates the local extreme in the specified direction using the supplied
   * instance of a derivation of OptimizerND. Function returns true on
   * successful local extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   *  \param dir Direction to search for local extreme
   */
  bool Extreme( VectorType & extX, double * extVal, VectorType & dir );

  /** Calculates the local extreme in the basis space directions using the
   * supplied instance of a derivation of OptimizerND. Function returns
   * true on successful local extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location
   *         of extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   *  \param n number of vectors in directions to use to define the basis
   *         space
   *  \param dirs TNT::Vectors that define the basis space to search for
   *         local extreme
   */
  bool Extreme( VectorType & extX, double * extVal, unsigned int n,
    MatrixType & dirs );

  /** Calculates the local extreme using an approximation to the conjugate
   * gradient descent method. Function returns true on successful local
   * extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of
   *         extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   */
  bool ExtremeConjGrad( VectorType & extX, double * extVal );

  /** Returns spline interpolated Hessian at x.
   * Calculates the values at control ( integer ) points by calling the
   * UserFunction and returns the interpolated Hessian between those points.
   * Type of interpolation is dependent on which spline derivation is used
   * ( e.g., SplApprox1D ).  Intermediate calculations and control point
   * evaluations are stored to speed subsequent calls
   */
  MatrixType & Hessian( const VectorType & x );

  /** Specify a new spline function
   *  \param dimension dimensionality of the space being interpolated
   *  \param funcVal function evaluated at integer points to determine
   *         control point values
   *  \param spline1D derivation of Spline1D ( e.g., SplApprox1D ) used
   *         to marginally ( i.e., per dimension ) interpolate values
   *  \param optimizer1D an instance ( can be NULL constructed ) of OptimizerND
   *         that is used to find local extrema
   *  \warning xMin and xMax must be set!
   */
  void Use( unsigned int dimension,
    ValueFunctionType::Pointer funcVal,
    Spline1D::Pointer spline1D,
    Optimizer1D::Pointer optimizer1D );

  /** Returns spline interpolated value at x.
   * Calculates the values at control ( integer ) points by calling the
   * UserFunction and returns the interpolated value between those points.
   * Type of interpolation is dependent on which spline derivation is used
   * ( e.g., SplApprox1D ).   Intermediate calculations and control point
   * evaluations are stored to speed subsequent calls.
   */
  double Value( const VectorType & x );

  /** Returns spline interpolated first derivative at x projected onto dx.
   *  Calculates the values at control ( integer ) points by calling the
   *  UserFunction and returns the interpolated first derivative between those
   *  points.  Type of interpolation is dependent on which spline derivation
   *  is used ( e.g., SplApprox1D ).   Intermediate calculations and control
   *  point evaluations are stored to speed subsequent calls.
   */
  double ValueD( const VectorType & x, IntVectorType & dx );

  /** Returns spline interpolated first derivative at x.
   * Calculates the values at control ( integer ) points by calling the
   * UserFunction and returns the interpolated first derivative between those
   * points.  Type of interpolation is dependent on which spline derivation
   * is used ( e.g., SplApprox1D ).  Intermediate calculations and control
   * point evaluations are stored to speed subsequent calls.
   */
  VectorType & ValueD( const VectorType & x );

  /** Returns spline interpolated derivative jet ( value, 1st derivative,
   * Hessian ) at x. Calculates the values at control ( integer ) points by
   * calling the UserFunction and returns the interpolated derivative jet
   * between those points.   Type of interpolation is dependent on which
   * spline derivation is used ( e.g., SplApprox1D ).   Intermediate
   * calculations and control point evaluations are stored to speed
   * subsequent calls.
   */
  double ValueJet( const VectorType & x, VectorType & d, MatrixType & h );

  /** Returns spline interpolated 1st derivatives and 2nd derivatives at x
   * Calculates the values at control ( integer ) points by calling the
   * UserFunction and returns the interpolated 1st derivatives and 2nd
   * derivatives between those points.  Type of interpolation is dependent
   * on which spline derivation is used ( e.g., SplApprox1D ). Intermediate
   * calculations and control point evaluations are stored to speed
   * subsequent calls
   */
  double ValueVDD2( const VectorType & x, VectorType & d, VectorType & d2 );

protected:

  typedef itk::VectorContainer< unsigned int, ImageType::Pointer >
    VectorImageType;

  /** Print out information about this object. */
  void PrintSelf( std::ostream & os, Indent indent ) const;

  void m_GetData( const VectorType & x );

  unsigned int                              m_Dimension;
  bool                                      m_Clip;
  IntVectorType                             m_XMin;
  IntVectorType                             m_XMax;
  bool                                      m_NewData;
  IntVectorType                             m_Xi;
  double                                    m_Val;
  VectorType                                m_D;
  MatrixType                                m_H;
  ImageType::Pointer                        m_Data;
  ImageType::Pointer                        m_DataWS;
  VectorType                                m_Data1D;
  VectorImageType::Pointer                  m_DataWSX;
  VectorImageType::Pointer                  m_DataWSXX;
  ValueFunctionType::Pointer                m_FuncVal;
  OptimizerValueFunctionType::Pointer       m_OptimizerNDVal;
  OptimizerDerivativeFunctionType::Pointer  m_OptimizerNDDeriv;
  OptimizerND::Pointer                      m_OptimizerND;
  Spline1D::Pointer                         m_Spline1D;

private:

  // Copy constructor not implemented.
  SplineND( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class SplineND

} // End namespace tube

#endif // End !defined( __tubeSplineND_h )
