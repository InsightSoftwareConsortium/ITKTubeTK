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

#ifndef __tubeMatrixMath_hxx
#define __tubeMatrixMath_hxx

#include "tubeMatrixMath.h"

#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/algo/vnl_matrix_inverse.h>

namespace tube
{

/**
 * Compute the orthogonal vector of one vector
 * Do not check for wrong dimension             */
template< class T >
vnl_vector<T>
ComputeOrthogonalVector(vnl_vector<T>  v)
{
  // create the orthogonal vector
  vnl_vector<T> n(v.size());

  /** if the dimension is 2
   *  purely arbitrary sign */
  if(v.size() == 2)
    {
    n(0)=v(1);
    n(1)=-v(0);
    }
  /** if the dimension is 3 */
  if(v.size() == 3)
    {
    vnl_vector< T > tt(3);

    tt[0] = -v[1];
    tt[1] = v[0];
    tt[2] = 0.5;
    tt.normalize();

    n = ComputeCrossVector( v, tt );
    }

  return n;
}


/**
 * Compute the cross vector in 3D resulting from two vectors
 * Do not check for wrong dimension
 * could use vnl_cross3D instead        */
template< class T >
vnl_vector<T>
ComputeCrossVector(vnl_vector<T> v1, vnl_vector<T> v2)
{
  vnl_vector<T> dest(v1.size());
  dest(0) = ((v1)(1) * (v2)(2)) - ((v1)(2) * (v2)(1));
  dest(1) = ((v1)(2) * (v2)(0)) - ((v1)(0) * (v2)(2));
  dest(2) = ((v1)(0) * (v2)(1)) - ((v1)(1) * (v2)(0));

  return dest;
}


/**
 *return the new position following the vector direction */
template< class T >
vnl_vector<T>
ComputeLineStep(vnl_vector<T> x, double a, vnl_vector<T> dir)
{
  int n = x.size();
  vnl_vector<T> destX(n);
  for(int i=0; i<n; i++)
    {
    destX(i) = (T)(x(i) + a * dir(i));
    }

  return destX;
}


/**
 * Compute the Euclidean distance  */
template< class T >
double
ComputeEuclideanDistanceVector(vnl_vector<T> x, const vnl_vector<T> y)
{
  double s = 0;
  for(unsigned int i=0; i<x.size(); i++)
    {
    s += (x(i)-y(i))*(x(i)-y(i));
    }
  return vcl_sqrt(s);
}

/**
 * Compute the Euclidean distance  */
template< class TPoint >
double
ComputeEuclideanDistance( TPoint x, TPoint y )
{
  double s = 0;
  for( unsigned int i = 0; i < TPoint::PointDimension; i++)
    {
    s += (x[i]-y[i])*(x[i]-y[i]);
    }
  return vcl_sqrt(s);
}

template< class T >
void
ComputeRidgeness( const vnl_matrix<T> & H,
  const vnl_vector<T> & D,
  double & ridgeness,
  double & roundness,
  double & curvature,
  double & levelness,
  vnl_matrix<T> & HEVect, vnl_vector<T> & HEVal )
{
  unsigned int ImageDimension = D.size();

  ::tube::ComputeEigen( H, HEVect, HEVal, true, false );
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    for( unsigned int j=0; j<ImageDimension; j++ )
      {
      std::cout << H(i, j) << " ";
      }
    std::cout << std::endl;
    }
  std::cout << std::endl;

  vnl_vector<T> Dv = D;
  if( Dv.magnitude() > 0 )
    {
    Dv.normalize();
    }
  else
    {
    Dv = HEVect.get_column( ImageDimension-1 );
    }

  double sump = 0;
  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    std::cout << HEVal[i] << " ";
    }
  std::cout << std::endl;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    double tf = 0;
    for( unsigned int j=0; j<ImageDimension; ++j )
      {
      tf += HEVect.get_column( i )[j] * Dv[j];
      }
    sump += tf * tf;

    tf = HEVal[i];
    sumv += tf * tf;
    if( tf >= 0 )
      {
      ridge = -1;
      }
    }
  double avgp = sump / ( ImageDimension - 1 );
  double avgv = sumv / ( ImageDimension - 1 );

  double curv = sumv
    / ( sumv + HEVal[ ImageDimension - 1 ] * HEVal[ ImageDimension - 1] );
  ridgeness = ridge * ( 1 - avgp ) * curv;

  // roundness is 1 - || 1 - v2^2 / v1^2  ||
  // However, since v1^2 can be extremely small, we use the average
  roundness = 0;
  if( avgv != 0 )
    {
    if( ImageDimension > 2 )
      {
      roundness =
        1 - vnl_math_abs( 1 - ( ( HEVal[ ImageDimension-2 ] *
          HEVal[ ImageDimension-2] ) / avgv ) );
      }
    else
      {
      double denom = sumv + ( HEVal[1] * HEVal[1] );
      if( denom != 0 )
        {
        roundness = ridge * (1 - ( HEVal[1] * HEVal[1] ) / denom );
        }
      }
    }

  // curvature is (v1^2 + v2^2) / 2.0
  curvature = 0;
  if( avgv != 0 )
    {
    if( ImageDimension > 2 )
      {
      // Multiplied by 50 assuming the image is from 0 to 1;
      curvature = ridge * vcl_sqrt( avgv ) * 50;
      }
    else
      {
      curvature = ridge * vcl_sqrt( avgv ) / 4;
      }
    }

  // levelness is (v1^2 + v2^2) / (v1^2 + v2^2 + v3^2) = 1 for a flat ridge
  levelness = 0;
  double denom =
     sumv + ( HEVal[ ImageDimension-1 ] * HEVal[ ImageDimension-1] );
  if( denom != 0 )
    {
    levelness = ridge * sumv / denom;
    }
}

/**
 * Compute eigenvalues and vectors from ( W.inv() * B ) */
template< class T >
void
ComputeEigenOfMatrixInvertedTimesMatrix(
  vnl_matrix<T> const & matToInvert, vnl_matrix<T> const & mat,
  vnl_matrix<T> &eVects, vnl_vector<T> &eVals,
  bool orderByAbs, bool minToMax )
{
  vnl_matrix<T> l, u, a;
  l = vnl_cholesky( matToInvert, vnl_cholesky::quiet ).lower_triangle();
  u = l.transpose();
  a = vnl_matrix_inverse<double>( l ) * mat
    * vnl_matrix_inverse<double>( u );
  ::tube::ComputeEigen( a, eVects, eVals, orderByAbs, minToMax );
  eVects = vnl_matrix_inverse<double>( u ) * eVects;
}

/**
 * Ensure the matrix is symmetric  */
template< class T >
void
FixMatrixSymmetry( vnl_matrix<T> & mat )
{
  for( unsigned int r=0; r<mat.rows(); ++r )
    {
    for( unsigned int c=r+1; c<mat.columns(); ++c )
      {
      mat( r, c ) = ( mat( r, c ) + mat( c, r ) ) / 2;
      mat( c, r ) = mat( r, c );
      }
    }
}

/**
 * Compute eigenvalues and vectors  */
template< class T >
void
ComputeEigen( vnl_matrix<T> const & mat,
  vnl_matrix<T> &eVects, vnl_vector<T> &eVals,
  bool orderByAbs, bool minToMax )
{

  unsigned int n = mat.rows();

  vnl_vector<T> subD(n);

  eVects = mat;
  eVals.set_size( n );
  bool symmetric = true;
  for( unsigned int i=0; i<n; ++i )
    {
    for( unsigned int j=i+1; j<n; ++j )
      {
      if( mat( i, j ) != mat( j, i ) )
        {
        i = n;
        std::cout
          << "Calling ComputeEigen with non-symmetric matrix is unstable."
          << std::endl;
        symmetric = false;
        break;
        }
      }
    }
  if( symmetric )
    {
    vnl_symmetric_eigensystem< T > eigen( mat );
    eVects = eigen.V;
    for( unsigned int d=0; d<eVects.columns(); d++ )
      {
      eVals[d] = eigen.get_eigenvalue( d );
      }
    }
  else
    {
    std::cout << "Non-symmetric matrix passed to eign-solver."
      << std::endl;
    vnl_matrix< double > matD( mat.rows(), mat.columns() );
    for( unsigned int c=0; c<mat.columns(); c++ )
      {
      for( unsigned int r=0; r<mat.rows(); r++ )
        {
        matD( r, c ) = mat( r, c );
        }
      }
    vnl_real_eigensystem eigen( matD );
    for( unsigned int c=0; c<mat.columns(); c++ )
      {
      eVals( c ) = eigen.D( c ).real();
      for( unsigned int r=0; r<mat.rows(); r++ )
        {
        eVects( r, c ) = eigen.Vreal(r, c);
        }
      }
    }

  if(orderByAbs)
    {
    for(unsigned int i=0; i<n-1; i++)
      {
      for(unsigned int j=i+1; j<n; j++)
        {
        if( ( vnl_math_abs(eVals(j))>vnl_math_abs(eVals(i)) && !minToMax )
          || ( vnl_math_abs(eVals(j))<vnl_math_abs(eVals(i)) && minToMax ) )
          {
          T tf = eVals(j);
          eVals(j) = eVals(i);
          eVals(i) = tf;
          vnl_vector<T> tv;
          tv = eVects.get_column( j );
          eVects.set_column( j, eVects.get_column( i ) );
          eVects.set_column( i, tv );
          }
        }
      }
    }
  else
    {
    for(unsigned int i=0; i<n-1; i++)
      {
      for(unsigned int j=i+1; j<n; j++)
        {
        if( ( eVals(j)>eVals(i) && !minToMax )
          || ( eVals(j)<eVals(i) && minToMax ) )
          {
          T tf = eVals(j);
          eVals(j) = eVals(i);
          eVals(i) = tf;
          vnl_vector<T> tv;
          tv = eVects.get_column( j );
          eVects.set_column( j, eVects.get_column( i ) );
          eVects.set_column( i, tv );
          }
        }
      }
    }
}

} // End namespace tube

#endif // End !defined(__tubeMatrixMath_hxx)
