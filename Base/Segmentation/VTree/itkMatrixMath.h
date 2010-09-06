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
#ifndef __itkMatrixMath_h
#define __itkMatrixMath_h

#include <vnl/vnl_vector_ref.h>
#define EIGEN_MAX_ITERATIONS 100

namespace itk 
{

/**
 * \brief This class provides usefull mathematical operations
 */

/** simple function that return the orthogonal vector to one*/
template <class T> 
vnl_vector<T>
GetOrthogonalVector(vnl_vector<T> x);

/** simple function that return the cross vector from two vectors*/
template <class T> 
vnl_vector<T>
GetCrossVector(vnl_vector<T> x);

/** return the new position folowing the vector dir */
template <class T> 
vnl_vector<T> 
ComputeLineStep(vnl_vector<T> x, double a, vnl_vector<T> dir);

/** Compute the euclidean distance */
template <class T> 
T
ComputeEuclideanDistanceVector(vnl_vector<T> x, const vnl_vector<T> y);

/** Compute eigen values and vectors  */
template <class T>
void
Eigen(vnl_matrix<T> &mat, vnl_matrix<T> &eVects,
  vnl_vector<T> &eVals, bool orderByAbs);

/** Preform trilinear diagonalisation */
template <class T>
void
TriDiag(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD);

/** Preform trilinear diagonalisation in 2D */
template <class T> 
void
TriDiag2D(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD);

/** Preform trilinear diagonalisation in 3D */
template <class T>
void
TriDiag3D(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD);

/**                                          */
template <class T>
void
Tqli(vnl_vector<T> &diag, vnl_vector<T> &subD, vnl_matrix<T> &mat);

/** Compute the euclidean distance for two points */
template <class PointType>
double
ComputeEuclideanDistance(PointType x, PointType y);

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMatrixMath.txx"
#endif

#endif /* __itkMatrixMath_h */
