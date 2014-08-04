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

#ifndef __itktubeTortuositySpatialObjectFilter_hxx
#define __itktubeTortuositySpatialObjectFilter_hxx

#include "itktubeTortuositySpatialObjectFilter.h"

namespace itk
{

namespace tube
{

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject >
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::TortuositySpatialObjectFilter( void )
{
  this->m_MeasureFlag = ALL;
  this->m_DistanceMetric = -1.0;
  this->m_InflectionCountMetric = -1.0;
  this->m_SumOfAnglesMetric = -1.0;
  this->m_InflectionPoints = itk::Array<double>();
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject >
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::~TortuositySpatialObjectFilter( void )
{
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > bool
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::IsUniqueMeasure(int flag)
{
  return flag == DISTANCE_METRIC ||
    flag == INFLECTION_COUNT_METRIC ||
    flag == SUM_OF_ANGLES_METRIC;
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > double
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GetMetric(int flag) const
{
  switch (flag)
    {
    case DISTANCE_METRIC:
      {
      return this->m_DistanceMetric;
      break;
      }
    case INFLECTION_COUNT_METRIC:
      {
      return this->m_InflectionCountMetric;
      break;
      }
    case SUM_OF_ANGLES_METRIC:
      {
      return this->m_SumOfAnglesMetric;
      break;
      }
    default:
      {
      return -1.0;
      break;
      }
    }
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > double
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GetDistanceMetric(int flag ) const
{
  return this->GetMetric(flag & DISTANCE_METRIC);
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > double
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GetInflectionCountMetric(int flag ) const
{
  return this->GetMetric(flag & INFLECTION_COUNT_METRIC);
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > double
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GetSumOfAnglesMetric(int flag ) const
{
  return this->GetMetric(flag & SUM_OF_ANGLES_METRIC);
}

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject > double
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GetInflectionPointValue( int i, int flag ) const
{
  if (! (flag & INFLECTION_POINTS) )
    {
    return -1.0;
    }
  return this->m_InflectionPoints.GetElement(i);
}

//------------------------------------------------------------------------------
namespace
{

double SafeAcos(double x)
{
  if (x < -1.0)
    {
    x = -1.0;
    }
  else if (x > 1.0)
    {
    x = 1.0;
    }
  return acos(x);
}

} // end namespace

//----------------------------------------------------------------------------
template< class TPointBasedSpatialObject >
void
TortuositySpatialObjectFilter< TPointBasedSpatialObject >
::GenerateData( void )
{
  PointBasedSpatialObject* output = this->GetOutput();
  const PointBasedSpatialObject* input = this->GetInput();

  typedef typename PointBasedSpatialObject::PointType PointType;
  typedef typename PointBasedSpatialObject::VectorType VectorType;

  if ( input->GetNumberOfPoints() < 2 )
    {
    itkExceptionMacro( << "Cannot run Tortuosity on input. "
                       << "Input has less than 2 points.");
    return;
    }

  bool noProblem = true;
  bool dm = this->m_MeasureFlag & DISTANCE_METRIC;
  bool icm = this->m_MeasureFlag & INFLECTION_COUNT_METRIC;
  bool ip = this->m_MeasureFlag & INFLECTION_POINTS;
  bool soam = this->m_MeasureFlag & SUM_OF_ANGLES_METRIC;

  //
  // DM variables
  VectorType start, end;
  double pathLength = 0.0;

  // ICM variables
  VectorType previousN(0.0);
  VectorType T(0.0), N(0.0), B(0.0); // for the Frenet frame
  bool vectorBIsValid = false;
  int inflectionCount = 1;
  
  // SOAM variables
  double totalCurvature = 0.0;

  size_t numberOfPoints = input->GetPoints().size();
  if (ip)
    {
    this->m_InflectionPoints.SetSize(numberOfPoints);
    }

  for(size_t index = 0; index < numberOfPoints; ++index)
    {
    const PointType* point =
      dynamic_cast<const PointType*>( input->GetPoint( index ) );
    assert( point );

    VectorType currentPoint =
      input->GetPoint( index )->GetPosition().GetVectorFromOrigin();

    //
    // General variables
    bool nextPointAvailable = ( index < numberOfPoints - 1 );
    VectorType nextPoint(0.0);
    if ( nextPointAvailable )
      {
      nextPoint =
        input->GetPoint( index + 1 )->GetPosition().GetVectorFromOrigin();
      }
    bool previousPointAvailable = ( index > 0 );
    VectorType previousPoint(0.0);
    if ( previousPointAvailable )
      {
      previousPoint =
        input->GetPoint( index - 1 )->GetPosition().GetVectorFromOrigin();
      }
    // t1 and t2, used both in icm and soam
    VectorType t1(0.0), t2(0.0);
    if ( previousPointAvailable && nextPointAvailable )
      {
      t1 = currentPoint - previousPoint;
      t2 = nextPoint - currentPoint;
      }

    bool nPlus2PointAvailable = ( index < numberOfPoints - 2 );
    VectorType nPlus2Point(0.0);
    if ( nPlus2PointAvailable )
      {
      nPlus2Point =
        input->GetPoint( index + 2 )->GetPosition().GetVectorFromOrigin();
      }

    //
    // DM Computations
    if (index == 0)
      {
      start = currentPoint;
      currentPoint = start;
      }
    if  (index == numberOfPoints - 1 )
      {
      end = currentPoint;
      }

    if ( ( dm || icm || ip || soam ) && nextPointAvailable )
      {
      pathLength += ( nextPoint - currentPoint ).GetNorm();
      }

    //
    // ICM Computations
    double inflectionValue = 0.0;
    if ( ( icm || ip ) && previousPointAvailable && nextPointAvailable )
      {
      // Compute velocity and acceleration
      VectorType v = nextPoint - previousPoint;
      VectorType a = t2 - t1;

      // Compute the Frenet frame
      // 1 - T = v / |v|
      T = v;
      T.Normalize();

      // 2 - N = v × a × v / |v × a × v|
      bool canCheckForinflection = a.GetNorm() > 1e-6;
      if ( canCheckForinflection )
        {
        N = CrossProduct( v, a );
        N = CrossProduct( v, N );
        N.Normalize();
        vectorBIsValid = true;
        }
      else if ( vectorBIsValid ) // 2nd chance
        {
        // Acceleration can be null when the curve approximates a straight
        // line (sin around pi for example). Unfortunately that could happen
        // when the curve is crossing the straight line and the inflection
        // would be missed...
        // This assumes that no pure torsion along the N vector happened.
        // Note that this is only valid is B was already computed at least
        // once.
        N = CrossProduct( B, T );
        N.Normalize();
        canCheckForinflection = true;
        }

      // 3 - B = T x N (in case of null acceleration. See above)
      B = CrossProduct( T, N );
      B.Normalize();

      if ( canCheckForinflection )
        {
        // Check for inflection
        VectorType deltaN = N - previousN;

        inflectionValue = deltaN * deltaN;
        if ( inflectionValue > 1.0 + 1e-6 )
          {
          inflectionCount += 1;
          }
        }

      N = previousN;
      }

    // Set the inflection value for this point
    if (ip)
      {
      this->m_InflectionPoints.SetElement(index, inflectionValue);
      }

    //
    // SOAM Computations
    if ( soam &&
       previousPointAvailable && nextPointAvailable && nPlus2PointAvailable )
      {
      VectorType t3 = nPlus2Point - nextPoint;

      // Compute in-plane angle
      VectorType normT1 = t1;
      normT1.Normalize();
      VectorType normT2 = t2;
      normT2.Normalize();

      double inPlaneAngle = SafeAcos( normT1 * normT2 );

      // Compute torsionnal angle
      VectorType t1t2Cross = CrossProduct( t1, t2 );
      VectorType t2t3Cross = CrossProduct( t2, t3 );

      double torsionAngle = 0.0;
      double t1t2Norm = t1t2Cross.GetNorm();
      double t2t3Norm = t2t3Cross.GetNorm();
      if ( t1t2Norm > 1e-6 || t2t3Norm > 1e-6 )
        {
        t1t2Cross *= 1.0 / t1t2Norm;
        t2t3Cross *= 1.0 / t2t3Norm;
        // We need to make sure we don't artificially scale those vectors
        // back to life. Otherwise we end up with a torsion angles where
        // there should not be one
        torsionAngle = SafeAcos( t1t2Cross * t2t3Cross );
        }
      // Finally get the curvature
      totalCurvature +=
        sqrt( inPlaneAngle*inPlaneAngle + torsionAngle*torsionAngle );
      }
    }

  //
  // DM final calculation
  double dmResult = -1.0;
  if ( dm || icm )
    {
    VectorType startToEnd = start- end;

    double straighLineLength = startToEnd.GetNorm();
    if ( straighLineLength > 0.0 )
      {
      this->m_DistanceMetric = pathLength / straighLineLength;
      }

    if ( dmResult < 1.0 )
      {
      std::cerr << "Error while computing the distance metric."
        << "DM (="<<dmResult<<") > 1.0" << std::endl;
      noProblem = false;
      }
    }

  // ICM final calculation
  this->m_InflectionCountMetric = inflectionCount * this->m_DistanceMetric;

  // SOAM final calculation
  double soamResult = 0.0;
  if ( soam && pathLength > 0.0 )
    {
    this->m_SumOfAnglesMetric = totalCurvature / pathLength;
    }
  else if ( soam && pathLength <= 0.0 )
    {
    std::cerr<<"Cannot compute SOAM, total tube path (="
      <<pathLength<<") <= 0.0"<<std::endl;
    noProblem = false;
    }

  if ( ! noProblem )
    {
    itkExceptionMacro( << "Problem while computing tortuosity.");
    return;
    }

  output->CopyInformation(input);
}

} // End namespace tube

} // End namespace itk

#endif // End !defined(__itktubeTortuositySpatialObjectFilter_hxx)
