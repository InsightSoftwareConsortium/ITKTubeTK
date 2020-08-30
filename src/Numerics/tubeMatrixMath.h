/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#ifndef __tubeMatrixMath_h
#define __tubeMatrixMath_h

#include "tubeMacro.h"

#include <itkVector.h>
#include <itkCovariantVector.h>

#include <vnl/vnl_math.h>
#include <vnl/vnl_vector_ref.h>

#define EIGEN_MAX_ITERATIONS 100

namespace tube
{

/**
 * \brief This class provides useful mathematical operations
 */

/** simple function that return the orthogonal vector to one */
template< class T >
vnl_vector<T>
ComputeOrthogonalVector( vnl_vector<T> x );

/** simple function that return the cross vector from two vectors */
template< class T >
vnl_vector<T>
ComputeCrossVector( vnl_vector<T> v1, vnl_vector<T> v2 );

template< class ScalarT=double >
void
ComputeNormalsFromTangents(
  const itk::Vector<ScalarT, 2> & prevT,
  const itk::CovariantVector<ScalarT, 2> & prevN1,
  const itk::Vector<ScalarT, 2> & t,
  itk::CovariantVector<ScalarT, 2> & n1 );

template< class ScalarT=double >
void
ComputeNormalsFromTangents(
  const itk::Vector<ScalarT, 3> & prevT,
  const itk::CovariantVector<ScalarT, 3> & prevN1,
  const itk::CovariantVector<ScalarT, 3> & prevN2,
  const itk::Vector<ScalarT, 3> & t,
  itk::CovariantVector<ScalarT, 3> & n1,
  itk::CovariantVector<ScalarT, 3> & n2 );

/** return the new position following the vector direction */
template< class T >
vnl_vector<T>
ComputeLineStep( vnl_vector<T> x, double a, vnl_vector<T> dir );

/** Compute the Euclidean distance */
template< class T >
double
ComputeEuclideanDistanceVector( vnl_vector<T> x, const vnl_vector<T> y );

/** Compute the Euclidean distance for two points */
template< class TPoint >
double
ComputeEuclideanDistance( TPoint x, TPoint y );

/** Compute Ridgeness measures */
template< class T >
void
ComputeRidgeness( const vnl_matrix<T> & H,
  const vnl_vector<T> & D,
  const vnl_vector<T> & prevTangent,
  double & ridgeness,
  double & roundness,
  double & curvature,
  double & linearity,
  vnl_matrix<T> & HEVect, vnl_vector<T> & HEVal );

/** Compute eigenvalues and vectors  */
template< class T >
void
FixMatrixSymmetry( vnl_matrix<T> & mat );

/** Compute eigenvalues and vectors  */
template< class T >
void
ComputeEigenOfMatrixInvertedTimesMatrix(
  vnl_matrix<T> const & matToInvert, vnl_matrix<T> const & mat,
  vnl_matrix<T> &eVects, vnl_vector<T> &eVals,
  bool orderByAbs, bool minToMax = true );

/** Compute eigenvalues and vectors  */
template< class T >
void
ComputeEigen( vnl_matrix<T> const & mat, vnl_matrix<T> &eVects,
  vnl_vector<T> &eVals, bool orderByAbs = false, bool minToMax = true );

} // End namespace tube


#ifndef TUBE_MANUAL_INSTANTIATION
#include "tubeMatrixMath.hxx"
#endif

#endif // End !defined( __tubeMatrixMath_h )
