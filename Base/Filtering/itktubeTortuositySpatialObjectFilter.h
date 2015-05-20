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

#ifndef __itktubeTortuositySpatialObjectFilter_h
#define __itktubeTortuositySpatialObjectFilter_h

#include "itktubeSpatialObjectToSpatialObjectFilter.h"

namespace itk
{

namespace tube
{

/** \class TortuositySpatialObjectFilter
 * \brief Compute tortuosity on a PointBasedSpatialObject.
 *
 * Compute tortuosity metrics on the given PointBasedSpatialObject. This does
 * not check for any child of the input. Different king of tortuosity can be
 * ran on the input, depending on the MeasureFlag.
 *
 * This filter does not modify the output so the filter can be part of a
 * pipeline and automatically recomputes itself when necessary.
 *
 */
template< class TPointBasedSpatialObject >
class TortuositySpatialObjectFilter
  : public SpatialObjectToSpatialObjectFilter< TPointBasedSpatialObject,
                                               TPointBasedSpatialObject >
{
public:
  /** Standard class typedefs. */
  typedef TortuositySpatialObjectFilter  Self;
  typedef SpatialObjectToSpatialObjectFilter< TPointBasedSpatialObject, TPointBasedSpatialObject >
    Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  typedef TPointBasedSpatialObject PointBasedSpatialObject;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( TortuositySpatialObjectFilter,
    SpatialObjectToSpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** The kind of metrics that can be ran on the input spatial object.
    * In details:
    *
    *  - Distance metric: Ratio of the spatial object path length
    * (sum of the line's length between each point) and the object direct
    * length (line length between the first and last point)
    *
    *  - Inflection count metric: Distance metric multiplied by the number of
    * times the object path crossed the first to last point line + 1.
    *  - Inflection points: Values checked to know when an inflection
    * occurs. Mostly useful for debugging.
    *
    *  - Sum of angles metric: Sum of the angles found along the object's path
    *
    *  - All: Compute all the metrics.
    *
    * NOTE: Contrary to the others, the Inflection points has one value for
    * each point of the vessels.
    */
  enum MeasureType
    {
    DISTANCE_METRIC = 0x01,
    INFLECTION_COUNT_METRIC = 0x02,
    INFLECTION_POINTS = 0x04,
    SUM_OF_ANGLES_METRIC = 0x08,
    ALL = 0xFF,
    };
  
  /** Set/Get the measure flag. This flag governs what metric is computed. */
  itkSetMacro( MeasureFlag, int )
  itkGetConstMacro( MeasureFlag, int )

  /** Return whether the given flag is only on measure or not. */
  bool IsUniqueMeasure( int flag );

  /** Return the given metric value. Before the filter is ran, those values
    * are invalid (-1.0). If a flag is given the value will be returned
    * only if the flag has the corresponding metric in it.
    */
  double GetMetric( int flag ) const;
  double GetDistanceMetric( int flag = DISTANCE_METRIC ) const;
  double GetInflectionCountMetric( int flag = INFLECTION_COUNT_METRIC ) const;
  double GetSumOfAnglesMetric( int flag = SUM_OF_ANGLES_METRIC ) const;

  /** Same than GetMetric for the Inflection Point metric. Since there are
    * values for each point of the object, an index must be given. No range
    * check is performed on that index.
    */
  double GetInflectionPointValue( int i, int flag = INFLECTION_POINTS ) const;


  /** Set/Get the sensibility of the filter. This is used in internal
    * computation when comparing whether a value should be considered null or
    * not. Default value (1e-6) should be acceptable for most cases.
    */
  itkSetMacro( Epsilon, double )
  itkGetConstMacro( Epsilon, double )

protected:
  TortuositySpatialObjectFilter( void );
  virtual ~TortuositySpatialObjectFilter( void );

  virtual void GenerateData( void );

private:
  TortuositySpatialObjectFilter( const Self & ); // purposely not implemented
  void operator=( const Self & ); // purposely not implemented

  int m_MeasureFlag;
  Array<double> m_InflectionPoints;
  double m_DistanceMetric;
  double m_InflectionCountMetric;
  double m_SumOfAnglesMetric;
  double m_Epsilon;

}; // End class TortuositySpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTortuositySpatialObjectFilter.hxx"
#endif

#endif // End !defined(__itktubeTortuositySpatialObjectFilter_h)
