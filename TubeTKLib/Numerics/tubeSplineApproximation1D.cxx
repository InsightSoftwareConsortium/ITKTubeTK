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

#include "tubeSplineApproximation1D.h"

namespace tube
{

SplineApproximation1D
::SplineApproximation1D( void )
  : Spline1D()
{
  m_SplineApproximation1DMatrixConst = ( double )( 1.0/6.0 );
  m_SplineApproximation1DMatrix( 0, 0 ) = 1;
  m_SplineApproximation1DMatrix( 0, 1 ) = 0;
  m_SplineApproximation1DMatrix( 0, 2 ) = 0;
  m_SplineApproximation1DMatrix( 0, 3 ) = 0;
  m_SplineApproximation1DMatrix( 1, 0 ) = -3;
  m_SplineApproximation1DMatrix( 1, 1 ) = 3;
  m_SplineApproximation1DMatrix( 1, 2 ) = 3;
  m_SplineApproximation1DMatrix( 1, 3 ) = 1;
  m_SplineApproximation1DMatrix( 2, 0 ) = 3;
  m_SplineApproximation1DMatrix( 2, 1 ) = -6;
  m_SplineApproximation1DMatrix( 2, 2 ) = 0;
  m_SplineApproximation1DMatrix( 2, 3 ) = 4;
  m_SplineApproximation1DMatrix( 3, 0 ) = -1;
  m_SplineApproximation1DMatrix( 3, 1 ) = 3;
  m_SplineApproximation1DMatrix( 3, 2 ) = -3;
  m_SplineApproximation1DMatrix( 3, 3 ) = 1;
}


SplineApproximation1D
::SplineApproximation1D( ValueFunctionType::Pointer funcVal,
  Optimizer1D::Pointer optimizer1D )
  : Spline1D( funcVal, optimizer1D )
{
  m_SplineApproximation1DMatrixConst = ( double )( 1.0/6.0 );
  m_SplineApproximation1DMatrix( 0, 0 ) = 1;
  m_SplineApproximation1DMatrix( 0, 1 ) = 0;
  m_SplineApproximation1DMatrix( 0, 2 ) = 0;
  m_SplineApproximation1DMatrix( 0, 3 ) = 0;
  m_SplineApproximation1DMatrix( 1, 0 ) = -3;
  m_SplineApproximation1DMatrix( 1, 1 ) = 3;
  m_SplineApproximation1DMatrix( 1, 2 ) = 3;
  m_SplineApproximation1DMatrix( 1, 3 ) = 1;
  m_SplineApproximation1DMatrix( 2, 0 ) = 3;
  m_SplineApproximation1DMatrix( 2, 1 ) = -6;
  m_SplineApproximation1DMatrix( 2, 2 ) = 0;
  m_SplineApproximation1DMatrix( 2, 3 ) = 4;
  m_SplineApproximation1DMatrix( 3, 0 ) = -1;
  m_SplineApproximation1DMatrix( 3, 1 ) = 3;
  m_SplineApproximation1DMatrix( 3, 2 ) = -3;
  m_SplineApproximation1DMatrix( 3, 3 ) = 1;
}


SplineApproximation1D
::~SplineApproximation1D( void )
{
}


double
SplineApproximation1D
::DataValue( const VectorType & y, double x )
{
  double u[4];
  u[3] = 1.0;
  u[2] = x-( int )x;
  u[1] = u[2]*u[2];
  u[0] = u[1]*u[2];

  double s = 0;
  for( unsigned int i=0; i<4; i++ )
    {
    double b = 0;
    for( unsigned int p=0; p<4; p++ )
      {
      b += m_SplineApproximation1DMatrix( i, p ) * u[p];
      }
    s += y( 3-i ) * b * m_SplineApproximation1DMatrixConst;
    }

  return s;
}


double
SplineApproximation1D
::DataValueD( const VectorType & y, double x )
{
  double u[3];
  u[2] = 1.0;
  u[1] = x-( int )x;
  u[0] = u[1]*u[1];

  double s = 0;
  for( unsigned int i=0; i<4; i++ )
    {
    double b = 0;
    for( unsigned int p=0; p<3; p++ )
      {
      b += ( 3-p )*m_SplineApproximation1DMatrix( i, p ) * u[p];
      }
    s += y( 3-i ) * b * m_SplineApproximation1DMatrixConst;
    }

  return s;
}


double
SplineApproximation1D
::DataValueD2( const VectorType & y, double x )
{
  double u[2];
  u[1] = 1.0;
  u[0] = x-( int )x;

  double s = 0;
  for( unsigned int i=0; i<4; i++ )
    {
    double b = 0;
    for( unsigned int p=0; p<2; p++ )
      {
      b += ( 2-p ) * m_SplineApproximation1DMatrix( i, p ) * u[p];
      }
    s += y( 3-i ) * b * m_SplineApproximation1DMatrixConst;
    }
  return s;
}


double
SplineApproximation1D
::DataValueJet( const VectorType & y, double x, double * d, double * d2 )
{

  double u[4];
  u[3] = 1.0;
  u[2] = x-( int )x;
  u[1] = u[2]*u[2];
  u[0] = u[1]*u[2];

  double s = 0;
  *d = 0;
  *d2 = 0;
  for( unsigned int i=0; i<4; i++ )
    {
    double b = 0;
    double bD = 0;
    double bD2 = 0;
    for( unsigned int p=0; p<4; p++ )
      {
      b += m_SplineApproximation1DMatrix( i, p ) * u[p];
      }
    for( unsigned int p=0; p<3; p++ )
      {
      bD += ( 3-p ) * m_SplineApproximation1DMatrix( i, p ) * u[p+1];
      }
    for( unsigned int p=0; p<2; p++ )
      {
      bD2 += ( 2-p ) * m_SplineApproximation1DMatrix( i, p ) * u[p+2];
      }
    s += y( 3-i ) * b * m_SplineApproximation1DMatrixConst;
    *d += y( 3-i ) * bD * m_SplineApproximation1DMatrixConst;
    *d2 += y( 3-i ) * bD2 * m_SplineApproximation1DMatrixConst;
    }

  return s;
}


void
SplineApproximation1D
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "SplineApproximation1DMatrixConst: "
     << m_SplineApproximation1DMatrixConst << std::endl;
  os << indent << "SplineApproximation1DMatrix:      "
     << m_SplineApproximation1DMatrix << std::endl;
}

} // End namespace tube
