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
#include <vnl/algo/vnl_svd.h>
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

  ::tube::ComputeEigen( H, HEVect, HEVal, false, true );

  vnl_vector<T> Dv = D;
  if( Dv.magnitude() > 0 )
    {
    Dv.normalize();
    }
  else
    {
    Dv = HEVect.get_column( ImageDimension-1 );
    }

  double P = dot_product( HEVect.get_column( ImageDimension-1 ), Dv );
  P = P * P;

  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    sumv += HEVal[i]*HEVal[i];
    if( HEVal[i] >= 0 )
      {
      ridge = -1;
      }
    }
  double avgv = sumv / ( ImageDimension - 1 );

  ridgeness = ridge * P;

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
      roundness = 1 - ( ( HEVal[0] * HEVal[0] ) - ( HEVal[1] * HEVal[1] ) )
        / ( ( HEVal[0] * HEVal[0] ) + ( HEVal[1] * HEVal[1] ) );
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
      double denom = sumv + ( HEVal[1] * HEVal[1] );
      if( denom != 0 )
        {
        curvature = ( HEVal[1] * HEVal[1] ) / denom;
        }
      }
    }

  // levelness is (v1^2 + v2^2) / (v1^2 + v2^2 + v3^2) = 1 for a flat ridge
  levelness = 0;
  double denom =
     sumv + ( HEVal[ ImageDimension-1 ] * HEVal[ ImageDimension-1] );
  if( denom != 0 )
    {
    levelness = sumv / denom;
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
  ::tube::ComputeEigen( a, eVects, eVals, orderByAbs, minToMax, true );
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
  bool orderByAbs, bool minToMax, bool useSVD )
{

  unsigned int n = mat.rows();

  vnl_vector<T> subD(n);

  eVects = mat;
  eVals.set_size( n );
  if( !useSVD )
    {
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
      switch(n)
        {
        case 1:
          eVects.set_size(1,1);
          eVects.fill( 1 );
          eVals.set_size(1);
          eVals.fill( mat[0][0] );
          break;
        case 2:
          ComputeTriDiag2D(eVects, eVals, subD);
          ComputeTqli(eVals, subD, eVects);
          break;
        case 3:
          ComputeTriDiag3D(eVects, eVals, subD);
          ComputeTqli(eVals, subD, eVects);
          break;
        default:
          vnl_symmetric_eigensystem< T > eigen( mat );
          eVects = eigen.V;
          for( unsigned int d=0; d<eVects.columns(); d++ )
            {
            eVals[d] = eigen.get_eigenvalue( d );
            }
          break;
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
    }
  else
    {
    vnl_svd< T > svd( mat );
    eVects = svd.U();
    eVals.set_size( eVects.columns() );
    for( unsigned int d=0; d<eVects.columns(); d++ )
      {
      eVals[d] = svd.W( d );
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
          eVects.get_column( j ) = eVects.get_column( i );
          eVects.get_column( i ) = tv;
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
          eVects.get_column( j ) = eVects.get_column( i );
          eVects.get_column( i ) = tv;
          }
        }
      }
    }
}


/**
 * Perform trilinear diagonalization in 2D */
template< class T >
void
ComputeTriDiag2D(vnl_matrix<T> &mat,
  vnl_vector<T> &diag, vnl_vector<T> &subD)
{
  diag(0) = mat(0,0);
  diag(1) = mat(1,1);
  subD(0) = mat(0,1);
  subD(1) = 0;

  mat(0,0) = 1;
  mat(0,1) = 0;

  mat(1,0) = 0;
  mat(1,1) = 1;
}

/**
 * Perform trilinear diagonalization in 3D */
template< class T >
void
ComputeTriDiag3D(vnl_matrix<T> &mat,
  vnl_vector<T> &diag, vnl_vector<T> &subD)
{
  double  a = mat(0,0), b = mat(0,1), c = mat(0,2),
          d = mat(1,1), e = mat(1,2), f = mat(2,2);

  diag(0) = static_cast< T >( a );
  subD(2) = 0;
  if(c != 0)
    {
    const double s = vcl_sqrt(b*b+c*c);
    b /= s;
    c /= s;
    const double q = 2*b*e+c*(f-d);
    diag(1) = static_cast< T >( d+c*q );
    diag(2) = static_cast< T >( f-c*q );
    subD(0) = static_cast< T >( s );
    subD(1) = static_cast< T >( e-b*q );

    mat(0,0) = 1;
    mat(0,1) = 0;
    mat(0,2) = 0;

    mat(1,0) = 0;
    mat(1,1) = static_cast< T >( b );
    mat(1,2) = static_cast< T >( c );

    mat(2,0) = 0;
    mat(2,1) = static_cast< T >( c );
    mat(2,2) = static_cast< T >( -b );
    }
  else
    {
    diag(1) = static_cast< T >( d );
    diag(2) = static_cast< T >( f );
    subD(0) = static_cast< T >( b );
    subD(1) = static_cast< T >( e );

    mat(0,0) = 1;
    mat(0,1) = 0;
    mat(0,2) = 0;

    mat(1,0) = 0;
    mat(1,1) = 1;
    mat(1,2) = 0;

    mat(2,0) = 0;
    mat(2,1) = 0;
    mat(2,2) = 1;
    }
}

template< class T >
void
ComputeTqli (vnl_vector<T> &diag, vnl_vector<T> &subD, vnl_matrix<T> &mat)
{
  int iter;
  int i;
  int k;
  int l;
  int m;

  double dd;
  double g;
  double r;
  double f;
  double s;
  double c;
  double p;
  double b;

  int n = mat.rows();

  for(l=0; l<n; l++)
    {
    for(iter = 0; iter < EIGEN_MAX_ITERATIONS; iter++)
      {
      for(m=l; m<n; m++)
        {
        if(m!=(n-1))
          {
          dd = vnl_math_abs(diag(m))+vnl_math_abs(diag(m+1));
          }
        else
          {
          dd = vnl_math_abs(diag(m));
          }

        if(vnl_math_abs(subD(m))+dd == dd)
          {
          break;
          }
        }
      if(m == l)
        {
        break;
        }
      g = (diag(l+1)-diag(l))/(2*subD(l));
      r = vcl_sqrt(g*g+1);
      if(g<0)
        {
        g = diag(m)-diag(l)+subD(l)/(g-r);
        }
      else
        {
        g = diag(m)-diag(l)+subD(l)/(g+r);
        }
      s = 1;
      c = 1;
      p = 0;
      for(i=m-1; i>=l; i--)
        {
        f = s*subD(i);
        b = c*subD(i);
        if(vnl_math_abs(f)>=vnl_math_abs(g))
          {
          c = g/f;
          r = vcl_sqrt(c*c+1);
          subD(i+1) = static_cast< T >( f*r );
          c *= (s = 1/r);
          }
        else
          {
          s = f/g;
          r = vcl_sqrt(s*s+1);
          subD(i+1) = static_cast< T >( g*r );
          s *= (c = 1/r);
          }
        g = diag(i+1)-p;
        r = (diag(i)-g)*s+2*b*c;
        p = s*r;
        diag(i+1) = static_cast< T >( g+p );
        g = c*r-b;

        for(k=0; k<n; k++)
          {
          f = mat(k,i+1);
          mat(k,i+1) = static_cast< T >( s*mat(k,i)+c*f );
          mat(k,i) = static_cast< T >( c*mat(k,i)-s*f );
          }
        }
      diag(l) -= static_cast< T >( p );
      subD(l) = static_cast< T >( g );
      subD(m) = 0;
      }
    if(iter == EIGEN_MAX_ITERATIONS)
      {
      throw("NR_tqli - exceeded maximum iterations\n");
      }
    }
}

} // End namespace tube

#endif // End !defined(__tubeMatrixMath_hxx)
