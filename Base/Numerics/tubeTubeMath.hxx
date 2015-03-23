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

#ifndef __tubeTubeMath_hxx
#define __tubeTubeMath_hxx

#include "tubeMatrixMath.h"
#include "tubeTubeMath.h"

namespace tube
{

/** Compute the tangent of the centerline of the tube */
template< class TTubePoint >
void
ComputeNormalsFromTangent( TTubePoint & tubePoint,
  const typename TTubePoint::VectorType & prevT )
{
  typedef typename TTubePoint::PointType             PointType;
  typedef typename TTubePoint::VectorType            VectorType;
  typedef typename TTubePoint::CovariantVectorType   CovariantVectorType;

  unsigned int dimension = tubePoint.GetPosition().GetPointDimension();

  VectorType t = tubePoint.GetTangent();

  if( t.GetNorm() == 0 )
    {
    if( prevT.GetNorm() != 0 )
      {
      t = prevT;
      t.Normalize();
      }
    else
      {
      t.Fill( 0 );
      t[0] = 1;
      }
    tubePoint.SetTangent(t);

    CovariantVectorType n;

    n.Fill( 0.0 );
    n[0] = 1;
    tubePoint.SetNormal1(n);

    if( dimension == 3 )
      {
      n.Fill( 0.0 );
      n[1] = 1;
      tubePoint.SetNormal2(n);
      }

    return;
    }

  // Compute the normal
  CovariantVectorType n1;
  CovariantVectorType n2;

  t = tubePoint.GetTangent();

  if( dimension == 2 )
    {
    n1.SetVnlVector( ::tube::ComputeOrthogonalVector( t.GetVnlVector() )
      .normalize() );
    if( n1 * prevT < 0 )
      {
      n1 *= -1;
      }
    tubePoint.SetNormal1( n1 );
    }
  else if( dimension == 3 )
    {
    CovariantVectorType tt;

    tt[0] = prevT[0];
    tt[1] = prevT[1];
    tt[2] = prevT[2];
    tt.Normalize();

    vnl_vector< double > vv = ::tube::ComputeCrossVector(
      t.GetVnlVector(), tt.GetVnlVector() );

    n2[0] = vv[0];
    n2[1] = vv[1];
    n2[2] = vv[2];
    n2.Normalize();

    vv = ::tube::ComputeCrossVector( t.GetVnlVector(),
      n2.GetVnlVector() );
    n1[0] = vv[0];
    n1[1] = vv[1];
    n1[2] = vv[2];
    n1.Normalize();

    VectorType tDiff = t - prevT;
    if( n1 * tDiff < 0 )
      {
      n1 *= -1;
      n2 *= -1;
      }
    tubePoint.SetNormal1( n1 );
    tubePoint.SetNormal2( n2 );
    }

  return;
}

/** Compute the tangent of the centerline of the tube */
template< class TTube >
bool
ComputeTubeTangentsAndNormals( typename TTube::Pointer & tube )
{
  return ComputeVectorTangentsAndNormals( tube->GetPoints() );
}

/** Compute the tangent of the centerline of the tube */
template< class TTubePoint >
bool
ComputeVectorTangentsAndNormals( std::vector< TTubePoint > & tubeV )
{
  typedef typename TTubePoint::PointType             PointType;
  typedef typename TTubePoint::VectorType            VectorType;
  typedef typename TTubePoint::CovariantVectorType   CovariantVectorType;

  unsigned int dimension = tubeV[0].GetPosition().GetPointDimension();

  int length = tubeV.size();
  if( length == 0 )
    {
    return false;
    }

  PointType x1, x3;
  VectorType t;
  t.Fill(0.0);

  if( length == 1 )
    {
    t[0] = 1;
    ::tube::ComputeNormalsFromTangent( tubeV[0], t );
    return true;
    }

  int it1 = 0;
  int it2 = 1;
  int it3 = 2;

  while(it3 < length)
    {
    x1 = tubeV[it1].GetPosition();
    x3 = tubeV[it3].GetPosition();
    double l=0;
    for(unsigned int i=0; i<dimension; i++)
      {
      t[i] = (x3[i] - x1[i]);
      l = l + t[i]*t[i];
      }

    l = vcl_sqrt(l);
    if(l < 0.0001)
      {
      std::cerr << "TubeSpatialObject::ComputeTangentAndNormals() : ";
      std::cerr << "length between two consecutive points is 0";
      std::cerr << " (use RemoveDuplicatePoints())" << std::endl;
      std::cerr << "   p1 = " << x1 << std::endl;
      std::cerr << "   p3 = " << x3 << std::endl;
      return false;
      }
    for(unsigned int i=0; i<dimension; i++)
      {
      t[i] /= l;
      }

    tubeV[it2].SetTangent(t);
    it1++;
    it2++;
    it3++;
    }

  it1 = 0;
  it2 = 1;
  t = tubeV[it2].GetTangent();
  tubeV[it1].SetTangent(t);
  it1 = length-1;
  it2 = length-2;
  t = tubeV[it2].GetTangent();
  tubeV[it1].SetTangent(t);


  // Compute the normal
  CovariantVectorType n1;
  CovariantVectorType n2;

  it1 = 0;
  VectorType prevT = tubeV[it1].GetTangent();
  while(it1 < length-1)
    {
    t = tubeV[it1].GetTangent();

    VectorType tt;
    it2 = it1+1;
    while( it2 < length )
      {
      tt = tubeV[it2].GetTangent();
      double diff = 0;
      for( unsigned int i=0; i<dimension; ++i )
        {
        diff += vnl_math_abs( tt[i] - t[i] );
        }
      if( diff > 0.0001 )
        {
        break;
        }
      ++it2;
      }
    if( it2 >= length )
      {
      tt[0] = t[1]+0.5;
      tt[1] = t[0]-0.5;
      tt[2] = t[2];
      tt.Normalize();
      std::cout << "Warning: Ilconditioned frenet frame, linear tube."
        << std::endl;
      }

    ::tube::ComputeNormalsFromTangent( tubeV[it1], tt );

    ++it1;
    }

  it1 = length-1;
  it2 = length-2;
  ::tube::ComputeNormalsFromTangent( tubeV[it1], tubeV[it2].GetTangent() );

  return true;
}

/**
 * Smooth a tube */
template< class TTube >
typename TTube::Pointer
SmoothTube( const typename TTube::Pointer & tube, double h,
  SmoothTubeFunctionEnum smoothFunction )
{
  typename TTube::PointType avg;

  typename TTube::PointListType & pointList = tube->GetPoints();
  typename TTube::PointListType::iterator pointItr;
  typename TTube::PointListType::iterator tmpPointItr;

  typename TTube::PointListType::iterator beginItr = pointList.begin();
  typename TTube::PointListType::iterator endItr = pointList.end();

  typename TTube::Pointer newTube = TTube::New();
  newTube->CopyInformation( tube );

  typename TTube::PointListType newPointList;

  if( h == 0 )
    {
    return newTube;
    }

  std::vector< double > w;
  int wSize = 0;

  if(smoothFunction == SMOOTH_TUBE_USING_INDEX_AVERAGE ||
     smoothFunction == SMOOTH_TUBE_USING_INDEX_GAUSSIAN )
    {
    // Calculate the weighing window w
    if( smoothFunction == SMOOTH_TUBE_USING_INDEX_AVERAGE)
      {
      int maxIndex = static_cast< int >( h );
      wSize = 2 * maxIndex + 1;
      w.resize(wSize, 1.0);
      }
    else
      {
      // Standard Deviation
      double sigma = h;
      // Consider the points until 3*sigma
      int maxIndex = static_cast< int >( 3*sigma );
      wSize = 2 * maxIndex + 1;
      w.resize(wSize, 0.0);
      for(int i = 0; i <= maxIndex; i++)
        {
        // The multiplication term 1/sigma*sqrt(2*pi) isn't necessary
        // since we normalize at the end by the sum of w
        w[maxIndex+i] = exp(-i*i/(2.0*sigma*sigma));
        w[maxIndex-i] = w[maxIndex+i];
        }
      }

    // Apply the weighing window
    int count = 0;
    unsigned int pointDimension = TTube::ObjectDimension;
    for( pointItr = beginItr; pointItr != endItr; ++pointItr )
      {
      typename TTube::TubePointType newPoint = *pointItr;
      double wTotal = 0;
      avg.Fill( 0 );
      tmpPointItr = pointItr;
      int wCenter = (wSize-1)/2;

      // Place the tmpPointItr at the beginning of the window
      tmpPointItr -= std::min( count, wCenter );

      // Place the window iterator so that the window center
      // is aligned to the current point.
      int pos = std::max(wCenter - count, 0);

      // Compute the average over the window, weighing with w
      while( pos < wSize && tmpPointItr != endItr )
        {
        for( unsigned int j=0; j<pointDimension; ++j)
          {
          avg[j] += w[pos] * tmpPointItr->GetPosition()[j];
          }
        wTotal += w[pos];
        ++pos;
        ++tmpPointItr;
        }

      // Divide by sum of weights -> finish average
      if( wTotal > 0 )
        {
        for( unsigned int i=0; i<pointDimension; ++i )
          {
          avg[i] /= wTotal;
          }
        // Update the new point coordinates
        newPoint.SetPosition( avg );
        }

      newPointList.push_back( newPoint );
      ++count;
      }
    }

  if( smoothFunction == SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN )
    {
    // // Set w for Gaussians
    // double sigma = h;
    // double dist = (pointItr->GetPosition()-tmpPointItr->GetPosition()).GetNorm();
    // w[pos] = 1/(sigma*2.50663)*exp(-dist*dist/(2.0*sigma*sigma));

    // TODO : Finish implementation
    std::cerr<<" Function not yet implemented, results you might get are false.\n";

    }

  newTube->SetPoints( newPointList );
  tube::ComputeTubeTangentsAndNormals< TTube >( newTube );

  return newTube;
}

} // End namespace tube

#endif // End !defined(__tubeTubeMath_hxx)
