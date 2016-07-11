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

#include "tubeTubeMath.h"

namespace itk
{

namespace tube
{

/** \class TortuositySpatialObjectFilter
* \brief Compute tortuosity on a PointBasedSpatialObject.
*
* Compute tortuosity metrics on the given PointBasedSpatialObject. This
* does not check for any child of the input. Different king of tortuosity
* can be run on the input, depending on the MeasureFlag.
*
* This filter does not modify the output so the filter can be part of a
* pipeline and automatically recomputes itself when necessary.
*/
template< class TPointBasedSpatialObject >
class TortuositySpatialObjectFilter : public
  SpatialObjectToSpatialObjectFilter< TPointBasedSpatialObject,
  TPointBasedSpatialObject >
{
public:
  /** Standard class typedefs. */
  typedef TortuositySpatialObjectFilter     Self;
  typedef SpatialObjectToSpatialObjectFilter< TPointBasedSpatialObject,
    TPointBasedSpatialObject >              Superclass;
  typedef SmartPointer< Self >              Pointer;
  typedef SmartPointer< const Self >        ConstPointer;

  typedef TPointBasedSpatialObject          PointBasedSpatialObject;

  typedef typename TPointBasedSpatialObject::VectorType SOVectorType;

  typedef typename PointBasedSpatialObject::Pointer
    PointBasedSpatialObjectPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro( TortuositySpatialObjectFilter,
    SpatialObjectToSpatialObjectFilter );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** The kind of metrics that can be ran on the input spatial object.
  * In details:
  *
  *
  *  **** VESSEL-WISE METRICS ****
  *
  *  - Average radius metric: average of the radius at each point of the
  *  vessel
  *
  *  - Distance metric: Ratio of the spatial object path length
  * (sum of the line's length between each point) and the object direct
  * length (length of the line between the first and last point)
  *
  *  - Inflection count metric: Distance metric multiplied by the number of
  * times the object path crossed the first to last point line + 1.
  *
  *  - Inflection count 1 metric: This inflection count method is computed
  * by thresholding the curvature scalar function along the vessel, and
  * counting the number of times the curve goes over the threshold. For
  * this 1st method, the threshold is set to lambda/sigma. sigma is the
  * smoothing scale (standard * deviation for gaussian smoothing), and
  * lambda is an arbitrary value, set to 1.5, because it works. See
  * Curvature scalar metric below.
  *
  *  - Inflection count 2 metric: This inflection countis computed by the
  * same method as Inflection count 1 metric, except that the threshold
  * value is the Percentile 95 metric. See Percentile 95 metric below.
  *
  *  - Chord length metric: Length of the line between the first and last
  * point
  *
  *  - Path length metric: Sum of distances between each point
  *
  *  - Percentile 95 metric: The 95th percentile of the histogram of the
  * curvature scalars along the vessel. It is the value for which 95% of
  * the points have curvature scalars below it. See the Curvature scalar
  * metric below.
  *
  *  - Sum of angles metric: Sum of the angles found along the object's path
  *
  *  - Sum of torsion metric: Sum of the torsion angles found along the
  * object's path
  *
  *  - Total curvature metric: Sum of the curvature scalar along the vessel.
  * See the Curvature scalar metric below.
  *
  *  - Total squared curvature metric: Sum of the squared curvature scalar
  * along the vessel. See the Curvature scalar metric below.
  *
  *
  *
  *   **** POINT-WISE METRICS ****
  *
  *  - Curvature scalar metric: Defined at each point as the norm of the
  * curvature vector. See the Curvature vector metric below. This metric
  * is used to calculate a lot of other metrics.
  *
  *  - Curvature vector metric: Defined at each point as the cross product
  * between the first and the second derivative of the spatial line. It is
  * a 3-components vector (x, y, z). Mostly used for debugging
  *
  *  - Inflection points metric: Values checked to know when an inflection
  * occurs. Mostly useful for debugging.
  *
  *
  *
  *   **** OTHER-WISE METRICS ****
  *
  *  - Curvature histogram metrics: This metric contains multiple metrics.
  * The histogram of the curvature scalar along the vessel is computed with
  * the number of bins, minimum and maximum passed in parameters. The
  * population in each bin is considered as a metric. This is intended to
  * be fed to a machine learning algorithm.
  *
  *
  *
  *   **** METRIC BITMASKS ****
  *
  *  Helpful bitmasks used when iterating through the enum to
  *  know of what type a metric is.
  *
  *  - All metrics: All metrics are in this bitmask
  *
  *  - Point-wise metrics: All point-wise metrics are in this bitmask
  *
  *  - Other-wise metrics: All other-wise metrics are in this bitmask
  *
  *
  * NOTE: For the metrics (not the bitmasks), the numbers are given by
  * date of addition, but the names are sorted by metric type, then
  * alphabetically.
  *
  */
  enum MeasureType
    {
    // Vessel-wise metrics
    AVERAGE_RADIUS_METRIC =               0x4000,
    CHORD_LENGTH_METRIC =                 0x0020,
    DISTANCE_METRIC =                     0x0001,
    INFLECTION_COUNT_METRIC =             0x0002,
    INFLECTION_COUNT_1_METRIC =           0x0400,
    INFLECTION_COUNT_2_METRIC =           0x0800,
    PATH_LENGTH_METRIC =                  0x0010,
    PERCENTILE_95_METRIC =                0x1000,
    SUM_OF_ANGLES_METRIC =                0x0008,
    SUM_OF_TORSION_METRIC =               0x8000,
    TOTAL_CURVATURE_METRIC =              0x0100,
    TOTAL_SQUARED_CURVATURE_METRIC =      0x0200,

    // Point-wise metrics
    CURVATURE_SCALAR_METRIC =             0x0040,
    CURVATURE_VECTOR_METRIC =             0x0080,
    INFLECTION_POINTS_METRIC =            0x0004,

    // Other-wise metrics
    CURVATURE_HISTOGRAM_METRICS =         0x2000,

    // Metric bitmasks
    BITMASK_ALL_METRICS =                 0xFFFF,
    BITMASK_VESSEL_WISE_METRICS =         0xDF3B,
    BITMASK_POINT_WISE_METRICS =          0x00C4,
    BITMASK_OTHER_WISE_METRICS =          0x2000,
    BITMASK_CURVATURE_METRICS =           0x3FC0
    };
  
  /** Set/Get the measure flag. This flag governs what metric computed. */
  itkSetMacro( MeasureFlag, int );
  itkGetConstMacro( MeasureFlag, int );

  /** Getters for vessel-wise metrics */
  double GetAverageRadiusMetric() const;
  double GetDistanceMetric() const;
  double GetInflectionCountMetric() const;
  double GetInflectionCount1Metric() const;
  double GetInflectionCount2Metric() const;
  double GetChordLengthMetric() const;
  double GetPathLengthMetric() const;
  double GetPercentile95Metric() const;
  double GetSumOfAnglesMetric() const;
  double GetSumOfTorsionMetric() const;
  double GetTotalCurvatureMetric() const;
  double GetTotalSquaredCurvatureMetric() const;

  /** Getters for point-wise metrics */
  double GetCurvatureScalarMetric( unsigned int i ) const;

  /** Same than GetMetric for the Inflection Point metric. Since there are
  * values for each point of the object, an index must be given.
  */
  double GetInflectionPointValue( unsigned int i ) const;

  /** Getters for other-wise metrics */
  int GetCurvatureHistogramMetric( unsigned int bin ) const;

  /** Set/Get the sensibility of the filter. Spacing values are usually
  * around 0.1. On some vessels, it can be lower, and it causes the
  * CURVATURE_SCALAR_METRIC to go abnormally high. We don't normalize by
  * spacing when it is lower than EpsilonForSpacing. Default value (1e-2)
  * should be acceptable for most cases.
  */
  itkSetClampMacro( EpsilonForSpacing, double, 0,
    std::numeric_limits<double>::max())
  itkGetConstMacro( EpsilonForSpacing, double )

  /** Set/Get the sensibility of the filter. This is used in internal
  * computation when comparing whether a value should be considered null or
  * not. Default value (1e-6) should be acceptable for most cases.
  */
  itkSetClampMacro( EpsilonForZero, double, 0,
    std::numeric_limits<double>::max())
  itkGetConstMacro( EpsilonForZero, double )

  /** Set/Get the maximum of the range of values the histogram must be
  * computed on. Used for CURVATURE_HISTOGRAM_METRICS. Default value is 1.
  */
  itkSetClampMacro( HistogramMax, double, 0,
    std::numeric_limits<double>::max() )
  itkGetConstMacro( HistogramMax, double )

  /** Set/Get the minimum of the range of values the histogram must be
  * computed on. Used for CURVATURE_HISTOGRAM_METRICS. Default value is 0.
  */
  itkSetClampMacro( HistogramMin, double, 0,
    std::numeric_limits<double>::max() )
  itkGetConstMacro( HistogramMin, double )

  /** Set/Get the arbitrary parameter lambda used to compute Inflection
  * count 1 metric. By default, it is set to 1.5, because it works well.
  * These methods are mainly for experimentation. See
  * INFLECTION_COUNT_1_METRIC
  */
  itkSetMacro( Lambda, double )
  itkGetConstMacro( Lambda, double )

  /** Set/Get the number of bins that the histogram must have.
  * Used for CURVATURE_HISTOGRAM_METRICS. Default value is 20.
  */
  itkSetClampMacro( NumberOfBins, size_t, 1,
    std::numeric_limits<size_t>::max() )
  itkGetConstMacro( NumberOfBins, size_t )

  /** Set/Get the smoothing method to be applied to the vessel before
  * computing the metrics. Default method is
  * SMOOTH_TUBE_USING_INDEX_GAUSSIAN.
  */
  itkSetMacro( SmoothingMethod, ::tube::SmoothTubeFunctionEnum  )
  itkGetConstMacro( SmoothingMethod, ::tube::SmoothTubeFunctionEnum  )

  /** Set/Get the smoothing scale to be applied to the vessel before
  * computing the metrics. This value has different meanings when using
  * different smoothing methods.
  * See ::tube::SmoothTube() for more details. Default value is 5.0.
  */
  itkSetClampMacro( SmoothingScale, double, 0,
    std::numeric_limits<double>::max())
  itkGetConstMacro( SmoothingScale, double )

  /** Set/Get the subsampling scale to be applied to the vessel before
  * computing the metrics.
  */
  itkSetClampMacro( SubsamplingScale, int, 0,
    std::numeric_limits<int>::max())
  itkGetConstMacro( SubsamplingScale, int )

protected:
  TortuositySpatialObjectFilter( void );
  virtual ~TortuositySpatialObjectFilter( void );

  virtual void GenerateData( void );

private:
  // purposely not implemented
  TortuositySpatialObjectFilter( const Self & );

  /** Input parameters */
  double                         m_EpsilonForSpacing;
  double                         m_EpsilonForZero;
  double                         m_HistogramMax;
  double                         m_HistogramMin;
  double                         m_Lambda;
  int                            m_MeasureFlag;
  size_t                         m_NumberOfBins;
  ::tube::SmoothTubeFunctionEnum m_SmoothingMethod;
  double                         m_SmoothingScale;
  int                            m_SubsamplingScale;

  /** Vessel parameters */
  size_t                         m_NumberOfPoints;
  int                            m_TubeID;

  /** Vessel-wise Metrics */
  double                         m_AverageRadiusMetric;
  double                         m_ChordLengthMetric;
  double                         m_DistanceMetric;
  double                         m_InflectionCountMetric;
  double                         m_InflectionCount1Metric;
  double                         m_InflectionCount2Metric;
  double                         m_PathLengthMetric;
  double                         m_Percentile95Metric;
  double                         m_SumOfAnglesMetric;
  double                         m_SumOfTorsionMetric;
  double                         m_TotalCurvatureMetric;
  double                         m_TotalSquaredCurvatureMetric;

  /** Point-wise Metrics */
  std::vector<double>            m_CurvatureScalar;
  std::vector<SOVectorType>      m_CurvatureVector;
  Array<double>                  m_InflectionPoints;

  /** Other-wise metrics */
  std::vector<int>               m_CurvatureHistogramMetrics;

}; // End class TortuositySpatialObjectFilter

} // End namespace tube

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itktubeTortuositySpatialObjectFilter.hxx"
#endif

#endif // End !defined(__itktubeTortuositySpatialObjectFilter_h)
