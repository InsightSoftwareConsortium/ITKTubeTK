/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itkTubeMath_txx
#define __itkTubeMath_txx

#include <cstdlib>
#include <iostream>
#include "itkTubeMath.h"


namespace itk
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
  double l;
  t.Fill(0.0);

  if( length == 1 )
    {
    ((TubePointType *)(tube->GetPoint(0)))->SetTangent(t);
    return true;
    }

  unsigned int it1 = 0;
  unsigned int it2 = 1;
  unsigned int it3 = 2;

  while(it3 < (unsigned int)length)
    {
    x1 = tube->GetPoint(it1)->GetPosition();
    x3 = tube->GetPoint(it3)->GetPosition();
    l=0;
    for(unsigned int i=0; i<dimension; i++)
      {
      t[i] = (x3[i] - x1[i]);
      l = l + t[i]*t[i];
      }

    l = vcl_sqrt(l);
    if(l == 0)
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
  while(it1 < (unsigned int)length)
    {
    t = ((TubePointType*)(tube->GetPoint(it1)))->GetTangent();

    if (dimension == 2)
      {
      t = ((TubePointType*)(tube->GetPoint(it1)))->GetTangent();
      n1[0] = -t[1];
      n1[1] = t[0];
      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
      }
    else if (dimension == 3)
      {
      CovariantVectorType tt;

      tt[0] = -t[1];
      tt[1] = t[0];
      tt[2] = 0.5;
      tt.Normalize();

      vnl_vector< double > vv = GetCrossVector( t.GetVnlVector(),
        tt.GetVnlVector() );
      n1[0] = vv[0];
      n1[1] = vv[1];
      n1[2] = vv[2];
      n1.Normalize();

      vv = GetCrossVector( t.GetVnlVector(), n1.GetVnlVector() );
      n2[0] = vv[0];
      n2[1] = vv[1];
      n2[2] = vv[2];
      n2.Normalize();

      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);
      ((TubePointType*)(tube->GetPoint(it1)))->SetNormal2(n2);
      }

    it1++;
    }

  it1 = 0;
  it2 = 1;
  n1 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal1();
  ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);

  if (dimension == 3)
    {
    n2 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal2();
    ((TubePointType*)(tube->GetPoint(it1)))->SetNormal2(n2);
    }

  it1 = length-1;
  it2 = length-2;
  n1 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal1();
  ((TubePointType*)(tube->GetPoint(it1)))->SetNormal1(n1);

  if (dimension == 3)
    {
    n2 = ((TubePointType*)(tube->GetPoint(it2)))->GetNormal2();
    ((TubePointType*)(tube->GetPoint(it1)))->SetNormal2(n2);
    }

  return true;
}


}; // end namespace itk

#endif
