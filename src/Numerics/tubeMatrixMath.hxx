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
ComputeOrthogonalVector( vnl_vector<T>  v )
{
  // create the orthogonal vector
  vnl_vector<T> n( v.size() );

  /** if the dimension is 2
   *  purely arbitrary sign */
  if( v.size() == 2 )
    {
    n( 0 )=v( 1 );
    n( 1 )=-v( 0 );
    }
  /** if the dimension is 3 */
  if( v.size() == 3 )
    {
    vnl_vector< T > tt( 3 );

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
ComputeCrossVector( vnl_vector<T> v1, vnl_vector<T> v2 )
{
  vnl_vector<T> dest( v1.size() );
  dest( 0 ) = ( ( v1 )( 1 ) * ( v2 )( 2 ) ) - ( ( v1 )( 2 ) * ( v2 )( 1 ) );
  dest( 1 ) = ( ( v1 )( 2 ) * ( v2 )( 0 ) ) - ( ( v1 )( 0 ) * ( v2 )( 2 ) );
  dest( 2 ) = ( ( v1 )( 0 ) * ( v2 )( 1 ) ) - ( ( v1 )( 1 ) * ( v2 )( 0 ) );

  return dest;
}


template < class TubePointT = itk::TubeSpatialObjectPoint<3> >
void
ComputeNormalsFromTangents( std::vector< TubePointT > )
{
}

template <class ScalarT>
void
ComputeNormalsFromTangents(
  const itk::Vector<ScalarT, 2> & prevT,
  const itk::CovariantVector<ScalarT, 2> & prevN1,
  const itk::Vector<ScalarT, 2> & t,
  itk::CovariantVector<ScalarT, 2> & n1
  )
{
  double l = 0;
  for (unsigned int i = 0; i < 2; i++)
  {
    l = l + t[i] * t[i];
  }
  l = std::sqrt(l);
  if (Math::AlmostEquals(l, 0.0) || std::isnan(l))
  {
    n1 = prevN1;
    return;
  }
    
  if (PointDimensionT == 2)
  {
    n1[0] = t[1];
    n1[1] = -t[0];
    l = 0;
    for (unsigned int i = 0; i < 2; i++)
    {
      l += n1[i] * prevN1[i];
    }
    if (l < 0)
    {
      n1 *= -1;
    }
  }
}

template <class ScalarT>
void
ComputeNormalsFromTangents(
  const itk::Vector<ScalarT, 3> & prevT,
  const itk::CovariantVector<ScalarT, 3> & prevN1,
  const itk::CovariantVector<ScalarT, 3> & prevN2,
  const itk::Vector<ScalarT, 3> & t,
  itk::CovariantVector<ScalarT, 3> & n1,
  itk::CovariantVector<ScalarT, 3> & n2
  )
{
  double l = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    l = l + t[i] * t[i];
  }
  l = std::sqrt(l);
  if (Math::AlmostEquals(l, 0.0) || std::isnan(l))
  {
    n1 = prevN1;
    n2 = prevN2;
    return;
  }
    
  // The normal to the tanget in 3D is the cross product of adjacent
  //   tangent directions.
  n1[0] = t[1] * prevT[2] - t[2] * prevT[1];
  n1[1] = t[2] * prevT[0] - t[0] * prevT[2];
  n1[2] = t[0] * prevT[1] - t[1] * prevT[0];

  if (n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2] == 0.0)
  {
    // if the normal is null, pick an orthogonal direction
    double d = std::sqrt(t[0] * t[0] + t[1] * t[1]);
    if (d != 0)
    {
      n1[0] = t[1] / d;
      n1[1] = -t[0] / d;
      n1[2] = 0;
      l = 0;
      for (unsigned int i = 0; i < 3; i++)
      {
        l += n1[i] * prevN1[i];
      }
      if (l < 0)
      {
        n1 *= -1;
      }
    }
    else
    {
      d = std::sqrt(t[1] * t[1] + t[2] * t[2]);
      if (d != 0)
      {
        n1[0] = 0;
        n1[1] = t[2] / d;
        n1[2] = -t[1] / d;
        l = 0;
        for (unsigned int i = 0; i < 3; i++)
        {
          l += n1[i] * prevN1[i];
        }
        if (l < 0)
        {
          n1 *= -1;
        }
      }
      else
      {
        n1 = prevN1;
      }
    }
  }

  // The second normal is the cross product of the tangent and the
  //   first normal
  n2[0] = t[1] * n1[2] - t[2] * n1[1];
  n2[1] = t[2] * n1[0] - t[0] * n1[2];
  n2[2] = t[0] * n1[1] - t[1] * n1[0];
  l = 0;
  for (unsigned int i = 0; i < 3; i++)
  {
    l += n2[i] * prevN2[i];
  }
  if (l < 0)
  {
    n2 *= -1;
  }
}


/**
 *return the new position following the vector direction */
template< class T >
vnl_vector<T>
ComputeLineStep( vnl_vector<T> x, double a, vnl_vector<T> dir )
{
  int n = x.size();
  vnl_vector<T> destX( n );
  for( int i=0; i<n; i++ )
    {
    destX( i ) = ( T )( x( i ) + a * dir( i ) );
    }

  return destX;
}


/**
 * Compute the Euclidean distance  */
template< class T >
double
ComputeEuclideanDistanceVector( vnl_vector<T> x, const vnl_vector<T> y )
{
  double s = 0;
  for( unsigned int i=0; i<x.size(); i++ )
    {
    s += ( x( i )-y( i ) )*( x( i )-y( i ) );
    }
  return std::sqrt( s );
}

/**
 * Compute the Euclidean distance  */
template< class TPoint >
double
ComputeEuclideanDistance( TPoint x, TPoint y )
{
  double s = 0;
  for( unsigned int i = 0; i < TPoint::PointDimension; i++ )
    {
    s += ( x[i]-y[i] )*( x[i]-y[i] );
    }
  return std::sqrt( s );
}

template< class T >
void
ComputeRidgeness( const vnl_matrix<T> & H,
  const vnl_vector<T> & D,
  const vnl_vector<T> & prevTangent,
  double & ridgeness,
  double & roundness,
  double & curvature,
  double & levelness,
  vnl_matrix<T> & HEVect, vnl_vector<T> & HEVal )
{
  unsigned int ImageDimension = D.size();

  vnl_matrix<T> HSym( H );
  ::tube::FixMatrixSymmetry( HSym );
  ::tube::ComputeEigen( HSym, HEVect, HEVal, true, false );

  vnl_vector<T> Dv = D;
  if( Dv.magnitude() > 0 )
    {
    Dv.normalize();
    }
  else
    {
    Dv = HEVect.get_column( ImageDimension-1 );
    }

  if( !prevTangent.empty() )
    {
    // Check dot product between the previous tangent and every eigenvector
    // to determine which one corresponds to the tangent direction.
    // Then reorder the matrix if needed to always have the tangent in the
    // last column.
    unsigned int closestV = 0;
    double closestVDProd = 0;
    for( unsigned int i=0; i<ImageDimension; i++ )
      {
      double dProd = std::fabs( dot_product( prevTangent,
        HEVect.get_column( i ) ) );
      if( dProd > closestVDProd )
        {
        closestV = i;
        closestVDProd = dProd;
        }
      }
    if( closestV != ImageDimension-1 )
      {
      //std::cout << "***********Mixing things up: Chosen evect = "
        //<< closestV << std::endl;
      //for( unsigned int i=0; i<ImageDimension; i++ )
        //{
        //std::cout << "   dotProd(" << i << ") = "
          //<< std::fabs( dot_product( prevTangent, HEVect.get_column( i ) ) )
          //<< std::endl;
        //}
      double tf = HEVal[closestV];
      HEVal[closestV] = HEVal[ImageDimension-1];
      HEVal[ImageDimension-1] = tf;
      vnl_vector< double > tv;
      tv = HEVect.get_column( closestV );
      HEVect.set_column( closestV, HEVect.get_column( ImageDimension-1 ) );
      HEVect.set_column( ImageDimension-1, tv );
      }
    if( dot_product( prevTangent, HEVect.get_column( ImageDimension-1 ) )
      < 0 )
      {
      for( unsigned int i=0; i<ImageDimension-1; ++i )
        {
        HEVect.get_column( ImageDimension-1 )[i] = -1 * HEVect.get_column(
          ImageDimension-1 )[i];
        }
      }
    }

  double sump = 0;
  double sumv = 0;
  int ridge = 1;
  for( unsigned int i=0; i<ImageDimension-1; i++ )
    {
    double dProd = dot_product( Dv, HEVect.get_column( i ) );
    sump += dProd * dProd;

    double tf = HEVal[i];
    sumv += tf * tf;
    if( tf >= 0 )
      {
      ridge = -1;
      }
    }
  double avgp = sump / ( ImageDimension - 1 );
  double avgv = sumv / ( ImageDimension - 1 );

  ridgeness = 0;
  curvature = 0;
  roundness = 0;
  levelness = 0;
  if( ridge > 0 )
    {
    ridgeness = ( 1 - avgp );

    // roundness is v2^2 / v1^2
    if( ImageDimension > 2 )
      {
      roundness = std::sqrt(
        ( HEVal[ ImageDimension-2 ] * HEVal[ ImageDimension-2] ) / avgv );
      }
    else
      {
      roundness = 1;
      }
  
    // curvature is sqrt( ( v1^2 + v2^2 ) / 2.0 )
    if( avgv != 0 )
      {
      curvature = std::sqrt( avgv );
      }
  
    // levelness is ( v1^2 + v2^2 ) / ( v1^2 + v2^2 + v3^2 )
    //     = 1 for a flat ridge
    double denom =
       sumv + ( HEVal[ ImageDimension-1 ] * HEVal[ ImageDimension-1] );
    if( denom != 0 )
      {
      levelness = sumv / denom;
      }
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
  vnl_matrix<T> l, u, a, li, ui;
  l = vnl_cholesky( matToInvert, vnl_cholesky::quiet ).lower_triangle();
  u = l.transpose();
  li = vnl_matrix_inverse<double>( l ).as_matrix();
  ui = vnl_matrix_inverse<double>( u ).as_matrix();
  a = li * mat * ui;
  ::tube::FixMatrixSymmetry( a );
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
 * Perform trilinear diagonalization in 2D */
template< class T >
void
ComputeTriDiag2D( vnl_matrix<T> &mat,
  vnl_vector<T> &diag, vnl_vector<T> &subD )
{
  diag( 0 ) = mat( 0, 0 );
  diag( 1 ) = mat( 1, 1 );
  subD( 0 ) = mat( 0, 1 );
  subD( 1 ) = 0;

  mat( 0, 0 ) = 1;
  mat( 0, 1 ) = 0;

  mat( 1, 0 ) = 0;
  mat( 1, 1 ) = 1;
}

/**
 * Perform trilinear diagonalization in 3D */
template< class T >
void
ComputeTriDiag3D( vnl_matrix<T> &mat,
  vnl_vector<T> &diag, vnl_vector<T> &subD )
{
  double a = mat( 0, 0 );
  double b = mat( 0, 1 );
  double c = mat( 0, 2 );
  double d = mat( 1, 1 );
  double e = mat( 1, 2 );
  double f = mat( 2, 2 );

  diag( 0 ) = static_cast< T >( a );
  subD( 2 ) = 0;
  if( c != 0 )
    {
    const double s = std::sqrt( b*b+c*c );
    b /= s;
    c /= s;
    const double q = 2*b*e+c*( f-d );
    diag( 1 ) = static_cast< T >( d+c*q );
    diag( 2 ) = static_cast< T >( f-c*q );
    subD( 0 ) = static_cast< T >( s );
    subD( 1 ) = static_cast< T >( e-b*q );

    mat( 0, 0 ) = 1;
    mat( 0, 1 ) = 0;
    mat( 0, 2 ) = 0;

    mat( 1, 0 ) = 0;
    mat( 1, 1 ) = static_cast< T >( b );
    mat( 1, 2 ) = static_cast< T >( c );

    mat( 2, 0 ) = 0;
    mat( 2, 1 ) = static_cast< T >( c );
    mat( 2, 2 ) = static_cast< T >( -b );
    }
  else
    {
    diag( 1 ) = static_cast< T >( d );
    diag( 2 ) = static_cast< T >( f );
    subD( 0 ) = static_cast< T >( b );
    subD( 1 ) = static_cast< T >( e );

    mat( 0, 0 ) = 1;
    mat( 0, 1 ) = 0;
    mat( 0, 2 ) = 0;

    mat( 1, 0 ) = 0;
    mat( 1, 1 ) = 1;
    mat( 1, 2 ) = 0;

    mat( 2, 0 ) = 0;
    mat( 2, 1 ) = 0;
    mat( 2, 2 ) = 1;
    }
}

template< class T >
void
ComputeTqli ( vnl_vector<T> &diag, vnl_vector<T> &subD, vnl_matrix<T> &mat )
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

  for( l=0; l<n; l++ )
    {
    for( iter = 0; iter < EIGEN_MAX_ITERATIONS; iter++ )
      {
      for( m=l; m<n; m++ )
        {
        if( m!=( n-1 ) )
          {
          dd = std::fabs( diag( m ) )+std::fabs( diag( m+1 ) );
          }
        else
          {
          dd = std::fabs( diag( m ) );
          }

        if( std::fabs( subD( m ) )+dd == dd )
          {
          break;
          }
        }
      if( m == l )
        {
        break;
        }
      g = ( diag( l+1 )-diag( l ) )/( 2*subD( l ) );
      r = std::sqrt( g*g+1 );
      if( g<0 )
        {
        g = diag( m )-diag( l )+subD( l )/( g-r );
        }
      else
        {
        g = diag( m )-diag( l )+subD( l )/( g+r );
        }
      s = 1;
      c = 1;
      p = 0;
      for( i=m-1; i>=l; i-- )
        {
        f = s*subD( i );
        b = c*subD( i );
        if( std::fabs( f )>=std::fabs( g ) )
          {
          c = g/f;
          r = std::sqrt( c*c+1 );
          subD( i+1 ) = static_cast< T >( f*r );
          c *= ( s = 1/r );
          }
        else
          {
          s = f/g;
          r = std::sqrt( s*s+1 );
          subD( i+1 ) = static_cast< T >( g*r );
          s *= ( c = 1/r );
          }
        g = diag( i+1 )-p;
        r = ( diag( i )-g )*s+2*b*c;
        p = s*r;
        diag( i+1 ) = static_cast< T >( g+p );
        g = c*r-b;

        for( k=0; k<n; k++ )
          {
          f = mat( k, i+1 );
          mat( k, i+1 ) = static_cast< T >( s*mat( k, i )+c*f );
          mat( k, i ) = static_cast< T >( c*mat( k, i )-s*f );
          }
        }
      diag( l ) -= static_cast< T >( p );
      subD( l ) = static_cast< T >( g );
      subD( m ) = 0;
      }
    if( iter == EIGEN_MAX_ITERATIONS )
      {
      throw( "NR_tqli - exceeded maximum iterations\n" );
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

  vnl_vector<T> subD( n );

  eVects = mat;
  eVals.set_size( n );
  bool symmetric = true;
  for( unsigned int i=0; i<n-1; ++i )
    {
    for( unsigned int j=i+1; j<n; ++j )
      {
      if( mat( i, j ) != mat( j, i ) )
        {
        std::cout << "Non-symmetric matrix passed to eign-solver."
          << std::endl;
        std::cout << "   " << mat( i, j ) << " != " << mat( j, i )
          << std::endl;
        symmetric = false;
        i = n - 1;
        break;
        }
      }
    }
  if( symmetric )
    {
    switch( n )
      {
      case 1:
        eVects.set_size( 1, 1 );
        eVects.fill( 1 );
        eVals.set_size( 1 );
        eVals.fill( mat[0][0] );
        break;
      case 2:
        ComputeTriDiag2D( eVects, eVals, subD );
        ComputeTqli( eVals, subD, eVects );
        break;
      case 3:
        ComputeTriDiag3D( eVects, eVals, subD );
        ComputeTqli( eVals, subD, eVects );
        break;
      default:
        vnl_symmetric_eigensystem< T > eigen( mat );
        eVects = eigen.V;
        eVals.set_size( n );
        for( unsigned int d=0; d<n; d++ )
          {
          eVals[d] = eigen.get_eigenvalue( d );
          }
        break;
      }
    }
  else
    {
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
        eVects( r, c ) = eigen.Vreal( r, c );
        }
      }
    }

  if( orderByAbs )
    {
    for( unsigned int i=0; i<n-1; i++ )
      {
      for( unsigned int j=i+1; j<n; j++ )
        {
        if( ( std::fabs( eVals( j ) )>std::fabs( eVals( i ) )
            && !minToMax )
          || ( std::fabs( eVals( j ) )<std::fabs( eVals( i ) )
            && minToMax ) )
          {
          T tf = eVals( j );
          eVals( j ) = eVals( i );
          eVals( i ) = tf;
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
    for( unsigned int i=0; i<n-1; i++ )
      {
      for( unsigned int j=i+1; j<n; j++ )
        {
        if( ( eVals( j )>eVals( i ) && !minToMax )
          || ( eVals( j )<eVals( i ) && minToMax ) )
          {
          T tf = eVals( j );
          eVals( j ) = eVals( i );
          eVals( i ) = tf;
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

#endif // End !defined( __tubeMatrixMath_hxx )
