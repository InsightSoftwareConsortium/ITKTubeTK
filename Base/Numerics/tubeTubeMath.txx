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

#ifndef __tubeTubeMath_txx
#define __tubeTubeMath_txx

#include <iostream>

#include "tubeMatrixMath.h"
#include "tubeTubeMath.h"

namespace tube
{

/** Compute the tangent of the centerline of the tube */
template< class TubeT >
bool  ComputeTubeTangentsAndNormals( TubeT * tube )
{
  typedef typename TubeT::PointType             PointType;
  typedef typename TubeT::TubePointType         TubePointType;
  typedef typename TubeT::VectorType            VectorType;
  typedef typename TubeT::CovariantVectorType   CovariantVectorType;

  unsigned int dimension = tube->GetObjectDimension();

  int length = tube->GetNumberOfPoints();
  if( length == 0 )
    {
    return false;
    }

  PointType x1, x3;
  VectorType t;
  t.Fill( 0.0 );

  if( length == 1 )
    {
    t[dimension-1] = 1;
    ((TubePointType *)(tube->GetPoint(0)))->SetTangent(t);
    t[dimension-1] = 0;
    CovariantVectorType n;
    n.Fill( 0.0 );
    n[0] = 1;
    ((TubePointType *)(tube->GetPoint(0)))->SetNormal1(n);
    if( dimension == 3 )
      {
      n.Fill( 0.0 );
      n[1] = 1;
      ((TubePointType *)(tube->GetPoint(0)))->SetNormal2(n);
      }
    return true;
    }

  unsigned int it1 = 0;
  unsigned int it2 = 1;
  unsigned int it3 = 2;

  while(it3 < (unsigned int)length)
    {
    x1 = tube->GetPoint(it1)->GetPosition();
    x3 = tube->GetPoint(it3)->GetPosition();
    double l=0;
    for(unsigned int i=0; i<dimension; i++)
      {
      t[i] = (x3[i] - x1[i]);
      l = l + t[i]*t[i];
      }

    l = vcl_sqrt(l);
    if(l <= 0.0001)
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

    ((TubePointType*)(tube->GetPoint(it2)))->SetTangent(t);
    it1++;
    it2++;
    it3++;
    }

  it1 = 0;
  it2 = 1;
  t = ((TubePointType*)(tube->GetPoint(it2)))->GetTangent();
  ((TubePointType*)(tube->GetPoint(it1)))->SetTangent(t);
  it1 = length-1;
  it2 = length-2;
  t = ((TubePointType*)(tube->GetPoint(it2)))->GetTangent();
  ((TubePointType*)(tube->GetPoint(it1)))->SetTangent(t);

  // Compute the normal
  CovariantVectorType n1;
  CovariantVectorType n2;

  it1 = 0;
  while( it1 < (unsigned int)length-1 )
    {
    t = ((TubePointType*)(tube->GetPoint(it1)))->GetTangent();

    if( dimension == 2 )
      {
      t = ((TubePointType*)(tube->GetPoint(it1)))->GetTangent();
      n1[0] = -t[1];
      n1[1] = t[0];
      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
      }
    else if( dimension == 3 )
      {
      CovariantVectorType tt;

      it2 = it1+1;
      while( it2 < (unsigned int)length )
        {
        tt[0] = ((TubePointType*)(tube->GetPoint(it2)))->GetTangent()[0];
        tt[1] = ((TubePointType*)(tube->GetPoint(it2)))->GetTangent()[1];
        tt[2] = ((TubePointType*)(tube->GetPoint(it2)))->GetTangent()[2];
        tt.Normalize();
        if( tt[0] != t[0] || tt[1] != t[1] || tt[2] != t[2] )
          {
          break;
          }
        ++it2;
        }
      if( it2 == (unsigned int)length )
        {
        tt[0] = t[0]+0.5;
        tt[1] = t[1]+0.5;
        tt[2] = t[2];
        tt.Normalize();
        }

      vnl_vector< double > vv = ::tube::ComputeCrossVector( t.GetVnlVector(),
        tt.GetVnlVector() );
      n2[0] = vv[0];
      n2[1] = vv[1];
      n2[2] = vv[2];
      n2.Normalize();

      vv = ::tube::ComputeCrossVector( t.GetVnlVector(), n2.GetVnlVector() );
      n1[0] = vv[0];
      n1[1] = vv[1];
      n1[2] = vv[2];
      n1.Normalize();

      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal2(n2);
      }

    it1++;
    }

  it1 = length-1;
  if( dimension == 2 )
    {
    t = ((TubePointType*)(tube->GetPoint(it1)))->GetTangent();
    n1[0] = -t[1];
    n1[1] = t[0];
    ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
    }
  else if( dimension == 3 )
    {
    it2 = length-2;
    n1 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal1();
    n2 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal2();
    ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
    ((TubePointType*)(tube->GetPoint(it1)))->SetNormal2(n2);
    }

  return true;
}

/** Compute the tangent of the centerline of the tube */
template< class TubePointT >
bool  ComputeVectorTangentsAndNormals( std::vector< TubePointT > & tubeV )
{
  typedef TubePointT                                 TubePointType;
  typedef typename TubePointT::PointType             PointType;
  typedef typename TubePointT::VectorType            VectorType;
  typedef typename TubePointT::CovariantVectorType   CovariantVectorType;

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
    t[dimension-1] = 1;
    tubeV[0].SetTangent(t);
    t[dimension-1] = 0;
    CovariantVectorType n;
    n.Fill( 0.0 );
    n[0] = 1;
    tubeV[0].SetNormal1(n);
    if( dimension == 3 )
      {
      n.Fill( 0.0 );
      n[1] = 1;
      tubeV[0].SetNormal2(n);
      }
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
  while(it1 < length-1)
    {
    t = tubeV[it1].GetTangent();

    if( dimension == 2 )
      {
      t = tubeV[it1].GetTangent();
      n1[0] = -t[1];
      n1[1] = t[0];
      tubeV[it1].SetNormal1(n1);
      }
    else if( dimension == 3 )
      {
      CovariantVectorType tt;

      it2 = it1+1;
      while( it2 < length )
        {
        tt[0] = tubeV[it2].GetTangent()[0];
        tt[1] = tubeV[it2].GetTangent()[1];
        tt[2] = tubeV[it2].GetTangent()[2];
        tt.Normalize();
        if( vnl_math_abs(tt[0] - t[0]) > 0.0001
          || vnl_math_abs(tt[1] - t[1]) > 0.0001
          || vnl_math_abs(tt[2] - t[2]) > 0.0001 )
          {
          break;
          }
        ++it2;
        }
      if( it2 >= length )
        {
        it2 = it1-1;
        while( it2 >= 0 )
          {
          tt[0] = tubeV[it2].GetTangent()[0];
          tt[1] = tubeV[it2].GetTangent()[1];
          tt[2] = tubeV[it2].GetTangent()[2];
          tt.Normalize();
          if( vnl_math_abs(tt[0] - t[0]) > 0.0001
            || vnl_math_abs(tt[1] - t[1]) > 0.0001
            || vnl_math_abs(tt[2] - t[2]) > 0.0001 )
            {
            break;
            }
          --it2;
          }
        if( it2 < 0 )
          {
          tt[0] = t[1]+0.5;
          tt[1] = t[0]-0.5;
          tt[2] = t[2];
          tt.Normalize();
          std::cout << "Warning: Ilconditioned frenet frame, linear tube."
            << std::endl;
          }
        }

      vnl_vector< double > vv = ::tube::ComputeCrossVector( t.GetVnlVector(),
        tt.GetVnlVector() );
      n2[0] = vv[0];
      n2[1] = vv[1];
      n2[2] = vv[2];
      n2.Normalize();
      if( vnl_math_abs(n2.GetNorm() - 1) > 0.0001 )
        {
        std::cout << "Warning: Normal direction no unit."
          << std::endl;
        std::cout << "  n2 = " << n2 << std::endl;
        }

      vv = ::tube::ComputeCrossVector( t.GetVnlVector(), n2.GetVnlVector() );
      n1[0] = vv[0];
      n1[1] = vv[1];
      n1[2] = vv[2];
      n1.Normalize();
      if( vnl_math_abs(n1.GetNorm() - 1) > 0.0001 )
        {
        std::cout << "Warning: Normal direction no unit."
          << std::endl;
        std::cout << "  n1 = " << n1 << std::endl;
        }

      tubeV[it1].SetNormal1(n1);
      tubeV[it1].SetNormal2(n2);
      }

    ++it1;
    }

  it1 = length-1;
  if( dimension == 2 )
    {
    t = tubeV[it1].GetTangent();
    n1[0] = -t[1];
    n1[1] = t[0];
    tubeV[it1].SetNormal1(n1);
    }
  else if( dimension == 3 )
    {
    it2 = length-2;
    n1 = tubeV[it2].GetNormal1();
    n2 = tubeV[it2].GetNormal2();
    tubeV[it1].SetNormal1(n1);
    tubeV[it1].SetNormal2(n2);
    }

  return true;
}

} // End namespace tube

#endif // End !defined(__tubeTubeMath_txx)
