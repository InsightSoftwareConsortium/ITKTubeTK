/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkSplineND.h,v $
  Language:  C++
  Date:      $Date: 2005/07/29 22:11:14 $
  Version:   $Revision: 1.6 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef itkSplineND_h
#define itkSplineND_h


#include "itkSpline1D.h"
#include "itkOptimizer1D.h"
#include "itkOptimizerND.h"

/************************************************

  NEED TO BE TEMPLATED OVER THE DIMENSION !!!!!
  problems in ValueVDD2

************************************************/

namespace itk 
{

/*! Multidimensional spline class
 *  Evaluates user supplied function to determine values at integers and will then
 *  interpolate (based on derivation used) non-integer values.   Includes methods for
 *  determining value, first derivative, hessian, and local extrema.   Relies on a
 *  derivation of Spline1D to specify how interpolation is performed.
 *  \author Stephen R. Aylward
 *  \date 11/21/99
 */

class SplineND
{
  
public :
   
  /** Typedef for the vector type used */
  typedef vnl_vector< double > VectorType;
  typedef vnl_vector_ref< double > VectorRefType; 
  
  /** Typedef for the vector type used */
  typedef vnl_vector< int > IntVectorType;

  /** Typedef for the matrix type used */
  typedef vnl_matrix< double > MatrixType;

  /** Typedef for the multidimensional structure (ie an image) */
  typedef Image< double, 4 >   ImageType;
 
  //! Default constructor: insufficient for using class
  SplineND();

 /** Constructor produces usable instance of this class
  *  \param newNDims dimensionality of the space being interpolated
  *  \param newFuncVal function evaluated at integer points to determine control point values
  *  \param newSpline1D derivation of Spline1D (e.g., SplApprox1D) used to marginally (i.e., per dimension)
  *      interpolate values
  *  \param newOptND an instance (can be NULL constructed) of OptimizerND that is used to find local extrema
  *  \warning xMin and xMax must be set!
  */
  SplineND(unsigned int newNDims,
    UserFunc< IntVectorType, double > * newFuncVal,
    Spline1D * newSpline1D,
    Optimizer1D * newOptND);

 /*! Destructor
  *  This class is not usually derived
  */
  virtual ~SplineND();
      
 /*! Specify a new spline function
  *  \param newNDims dimensionality of the space being interpolated
  *  \param newFuncVal function evaluated at integer points to determine control point values
  *  \param newSpline1D derivation of Spline1D (e.g., SplApprox1D) used to marginally (i.e., per dimension)
  *      interpolate values
  *  \param newOptND an instance (can be NULL constructed) of OptimizerND that is used to find local extrema
  *  \warning xMin and xMax must be set!
  */
  void use(unsigned int newNDims,
    UserFunc<IntVectorType, double > * newFuncVal,
    Spline1D * newSpline1D,
    Optimizer1D * newOptND);
     
  //! Returns the characteristics of spline evaluations near data bounds (xMin and xMax)
 /*! If true, values beyond edges (xMin and xMax) are set to zero.
  *  If false, values beyond edges are faded to 0 as a function of distance from
  *      edge, squeared.
  */
  bool clipEdge();

  int nDims()
    {
    return cNDims;
    }

  //! User specification of characteristics
  void clipEdge(bool newClip);

  //! Returns the control points' (integer value locations) lower bound
  IntVectorType & xMin();

  //! User specification of lower bound
  void xMin(IntVectorType newXMin);

  //! Sets control points' (integer value locations) upper bound
  IntVectorType & xMax();

  //! User Specification of upper bound
  void xMax(IntVectorType newXMax);

  //! Tracks the validity of internally maintained intermediate calculations and data
  /*! Returns true if a new spline instance has been created (e.g., use has been called)
   */
  bool newData();

  //! User sets to true to force recalcuation of internal data
  /*! For example, use to flag that UserFunc has changed externally
   */
  void newData(bool newNewData);

  //! Returns spline interpolated value at x
  /*! Calculates the values at control (integer) points by calling the UserFunc and
   *  returns the interpolated value between those points.   Type of interpolation is
   *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
   *  calculations and control point evaluations are stored to speed subsequent calls
   */
  const double & value(const VectorType & x);

 //! Returns spline interpolated first derivative at x projected onto dx
 /*! Calculates the values at control (integer) points by calling the UserFunc and
  *  returns the interpolated first derivative between those points.   Type of interpolation is
  *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
  *  calculations and control point evaluations are stored to speed subsequent calls
  */
  double valueD(const VectorType & x, IntVectorType & dx);
  //! Returns spline interpolated first derivative at x
  /*! Calculates the values at control (integer) points by calling the UserFunc and
   *  returns the interpolated first derivative between those points.   Type of interpolation is
   *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
   *  calculations and control point evaluations are stored to speed subsequent calls
   */
  VectorType & valueD(const VectorType & x);

 //! Returns spline interpolated Hessian at x
 /*! Calculates the values at control (integer) points by calling the UserFunc and
  *  returns the interpolated Hessian between those points.   Type of interpolation is
  *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
  *  calculations and control point evaluations are stored to speed subsequent calls
  */
  MatrixType & hessian(VectorType & x);

  //! Returns spline interpolated derivative jet (value, 1st deriv, Hessian) at x
  /*! Calculates the values at control (integer) points by calling the UserFunc and
   *  returns the interpolated derivative jet between those points.   Type of interpolation is
   *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
   *  calculations and control point evaluations are stored to speed subsequent calls
   */
  double valueJet(VectorRefType x, VectorRefType d, MatrixType & h);

  //! Returns spline interpolated 1st derivatives and 2nd derivatives at x
 /*! Calculates the values at control (integer) points by calling the UserFunc and
  *  returns the interpolated 1st derivatives and 2nd derivatives between those points.  
  *  Type of interpolation is
  *  dependent on which spline derivation is used (e.g., SplApprox1D).   Intermediate
  *  calculations and control point evaluations are stored to speed subsequent calls
  */
  double valueVDD2(VectorType & x, VectorType & d, VectorType & d2);
  //! Calculates the local extreme using the supplied instance of a derivation of OptimizerND
 /*! Function returns true on successful local extreme finding, false otherwise.
  *  \param extX User supplied initial point, On return equals location of extreme local to initial point
  *  \param extVal On return equals the value at the local extreme
  */
  bool extreme(VectorRefType extX, double * extVal);

  //! Calculates the local extreme in the direction dir using the supplied instance of a derivation of OptimizerND
  /*! Function returns true on successful local extreme finding, false otherwise.
   *  \param extX User supplied initial point, On return equals location of extreme local to initial point
   *  \param extVal On return equals the value at the local extreme
   *  \param dir Direction to search for local extreme
   */
  bool extreme(VectorRefType extX, double * extVal, VectorType &dir);

 /** Calculates the local extreme in the basis space dirs using the supplied instance of a derivation of OptimizerND
  *  Function returns true on successful local extreme finding, false otherwise.
  *  \param extX User supplied initial point, On return equals location of extreme local to initial point
  *  \param extVal On return equals the value at the local extreme
  *  \param n number of vectors in dirs to use to define the basis space
  *  \param dirs TNT::Vectors that define the basis space to search for local extreme
  */
  bool extreme(VectorRefType extX, double * extVal, unsigned int n, MatrixType &dirs);
  
 /** Calculates the local extreme using an approximation to the conjugate gradient descent method
  * Function returns true on successful local extreme finding, false otherwise.
  *  \param extX User supplied initial point, On return equals location of extreme local to initial point
  *  \param extVal On return equals the value at the local extreme
  */
  bool extremeConjGrad(VectorType & extX, double * extVal);

protected : 
  
  /** Used to enable/disable cout of intermediate calculations */
  bool           m_debug;

  bool           cDefined;
  unsigned int   cNDims;
  bool           cClip;
  IntVectorType  cXMin;
  IntVectorType  cXMax;
  bool           cNewData;
  IntVectorType  cXi;
  VectorType     cD;
  MatrixType     cH;

  ImageType::Pointer  cData;
  ImageType::Pointer  cDataWS;
  VectorType          cData1D;

  UserFunc< IntVectorType, double >   * cFuncVal;

  UserFunc< VectorType, double >      * cOptNDVal;
  UserFunc< VectorType, VectorType >  * cOptNDDeriv;

  OptimizerND       * cOptND;
  Spline1D          * cSpline1D;     

  void  cGetData(const VectorType &x);  
};



} //end namespace itk



#endif

