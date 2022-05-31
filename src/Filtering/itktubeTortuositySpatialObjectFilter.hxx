/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __itktubeTortuositySpatialObjectFilter_hxx
#define __itktubeTortuositySpatialObjectFilter_hxx


#include "itkSampleToHistogramFilter.h"
#include "itkListSample.h"
#include "itkHistogram.h"
#include "itkVector.h"

#include "tubeTubeMathFilters.h"

namespace itk
{

namespace tube
{

//-------------------------------------------------------------------------
template< class TTubeSpatialObject >
TortuositySpatialObjectFilter< TTubeSpatialObject >
::TortuositySpatialObjectFilter( void )
{
  // Setting default parameters
  this->m_EpsilonForSpacing = 1e-2;
  this->m_EpsilonForZero = 1e-6;
  this->m_HistogramMin = 0;
  this->m_HistogramMax = 1;
  this->m_Lambda = 1.5;
  this->m_MeasureFlag = BITMASK_ALL_METRICS;
  this->m_NumberOfBins = 20;
  this->m_SmoothingScale = 5.0;
  this->m_SubsamplingScale = 1;

  // Setting vessel-wise metrics to -1.0
  this->m_AverageRadiusMetric = -1.0;
  this->m_ChordLengthMetric = -1.0;
  this->m_DistanceMetric = -1.0;
  this->m_InflectionCountMetric = -1.0;
  this->m_InflectionCount1Metric = -1.0;
  this->m_InflectionCount1Metric = -1.0;
  this->m_PathLengthMetric = -1.0;
  this->m_Percentile95Metric = -1.0;
  this->m_SumOfAnglesMetric = -1.0;
  this->m_SumOfTorsionMetric = -1.0;
  this->m_TotalCurvatureMetric = -1.0;
  this->m_TotalSquaredCurvatureMetric = -1.0;

  // Initializing point-wise metric arrays
  this->m_InflectionPoints = itk::Array<double>();
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject >
TortuositySpatialObjectFilter< TTubeSpatialObject >
::~TortuositySpatialObjectFilter( void )
{
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetAverageRadiusMetric() const
{
  return this->m_AverageRadiusMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetChordLengthMetric() const
{
  return this->m_ChordLengthMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetDistanceMetric() const
{
  return this->m_DistanceMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetInflectionCountMetric() const
{
  return this->m_InflectionCountMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetInflectionCount1Metric() const
{
  return this->m_InflectionCount1Metric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetInflectionCount2Metric() const
{
  return this->m_InflectionCount2Metric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetPathLengthMetric() const
{
  return this->m_PathLengthMetric;
}


//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetPercentile95Metric() const
{
  return this->m_Percentile95Metric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetSumOfAnglesMetric() const
{
  return this->m_SumOfAnglesMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetSumOfTorsionMetric() const
{
  return this->m_SumOfTorsionMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetTotalCurvatureMetric() const
{
  return this->m_TotalCurvatureMetric;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetTotalSquaredCurvatureMetric() const
{
  return this->m_TotalSquaredCurvatureMetric;

}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetCurvatureScalarMetric( unsigned int i ) const
{
  if( this->m_CurvatureScalar.size() == 0 )
    {
    std::cerr << "CurvatureScalarMetric not computed" << std::endl;
    }
  else if( i >= this->m_CurvatureScalar.size() )
    {
    std::cerr << "GetCurvatureScalarMetric( int ): Index " << i
      << " out of bounds" << std::endl;
    }
  else
    {
    return this->m_CurvatureScalar.at( i );
    }
  return -1;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > double
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetInflectionPointValue( unsigned int i ) const
{
  if( this->m_InflectionPoints.Size() == 0 )
    {
    std::cerr << "InflectionPointMetric not computed" << std::endl;
    }
  else if( i >= this->m_InflectionPoints.Size() )
    {
    std::cerr << "GetInflectionPointValue( int ): Index " << i
      << " out of bounds" <<std::endl;
    }
  else
    {
    return this->m_InflectionPoints.GetElement( i );
    }
  return -1;
}

//--------------------------------------------------------------------------
template< class TTubeSpatialObject > int
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GetCurvatureHistogramMetric( unsigned int bin ) const
{
  if( this->m_CurvatureHistogramMetrics.size() == 0 )
    {
    std::cerr << "CurvatureHistogramMetric not computed" << std::endl;
    return 0;
    }
  else if( bin >= this->m_CurvatureHistogramMetrics.size() )
    {
    std::cerr << "GetHistogramMetric( int ): Index " << bin
      << " out of bounds" << std::endl;
    return 0;
    }
  else
    {
    return this->m_CurvatureHistogramMetrics.at( bin );
    }
}

//--------------------------------------------------------------------------
namespace
{

double SafeAcos( double x )
{
  if( x < -1.0 )
    {
    x = -1.0;
    }
  else if( x > 1.0 )
    {
    x = 1.0;
    }
  return acos( x );
}

template<typename VectorType >
typename VectorType::RealValueType SafeNormalize( VectorType& v )
{
  typedef typename VectorType::RealValueType  RealType;
  typedef typename VectorType::ValueType      ValueType;

  typename VectorType::RealValueType norm = v.GetNorm();
  if( norm != 0.0 )
    {
    for( unsigned int i = 0; i < VectorType::Dimension; ++i )
      {
      v[i] = static_cast< ValueType >( static_cast< RealType >( v[i] )
        / norm );
      }
    }
  return norm;
}

} // end namespace

//--------------------------------------------------------------------------
template< class TTubeSpatialObject >
void
TortuositySpatialObjectFilter< TTubeSpatialObject >
::GenerateData( void )
{
  // Get I/O
  TubeSpatialObjectPointer output = this->GetOutput();

  TubeSpatialObjectPointer originalInput = TubeSpatialObject::New();
  originalInput = const_cast< TubeSpatialObject * >( this->GetInput() );

  // Safety check
  if( originalInput->GetNumberOfPoints() < 2 )
    {
    itkExceptionMacro( << "Cannot run Tortuosity on input. "
                       << "Input has less than 2 points." );
    return;
    }

  // Determine metrics to compute
  bool ipm = this->m_MeasureFlag & INFLECTION_POINTS_METRIC;
  bool arm = this->m_MeasureFlag & AVERAGE_RADIUS_METRIC;
  bool clm = this->m_MeasureFlag & CHORD_LENGTH_METRIC;
  bool dm = this->m_MeasureFlag & DISTANCE_METRIC;
  bool icm = this->m_MeasureFlag & INFLECTION_COUNT_METRIC;
  bool ic1m = this->m_MeasureFlag & INFLECTION_COUNT_1_METRIC;
  bool ic2m = this->m_MeasureFlag & INFLECTION_COUNT_2_METRIC;
  bool plm = this->m_MeasureFlag & PATH_LENGTH_METRIC;
  bool p95m = this->m_MeasureFlag & PERCENTILE_95_METRIC;
  bool soam = this->m_MeasureFlag & SUM_OF_ANGLES_METRIC;
  bool sot = this->m_MeasureFlag & SUM_OF_TORSION_METRIC;
  bool tcm = this->m_MeasureFlag & TOTAL_CURVATURE_METRIC;
  bool tscm = this->m_MeasureFlag & TOTAL_SQUARED_CURVATURE_METRIC;
  bool cvm = this->m_MeasureFlag & CURVATURE_VECTOR_METRIC;
  bool chm = this->m_MeasureFlag & CURVATURE_HISTOGRAM_METRICS;

  // DM variables
  SOVectorType start;
  SOVectorType end;
  double pathLength = 0.0;

  // ICM variables
  SOVectorType previousN( 0.0 );
  SOVectorType T( 0.0 ), N( 0.0 ), B( 0.0 ); // for the Frenet frame
  bool vectorBIsValid = false;
  int inflectionCount = 1;

  // SOAM variables
  double sumOfAngles = 0.0;
  double sumOfTorsion = 0.0;

  // Other metrics variables
  double totalCurvature = 0.0;
  double totalSquaredCurvature = 0.0;
  double sumOfRadius = 0.0;

  // Smooth the vessel
  ::tube::TubeMathFilters<TubeSpatialObject::ObjectDimension> filter;
  filter.SetInputTube( originalInput );
  filter.SmoothTube( this->m_SmoothingScale );

  // Subsample the vessel
  filter.SubsampleTube( this->m_SubsamplingScale );

  // Make the measurements on the pre-processed tube.
  TubeSpatialObjectPointer processedInput = filter.GetOutputTube();

  if( processedInput->GetNumberOfPoints() < 2 )
    {
    itkExceptionMacro( << "Cannot run Tortuosity on input. "
                       << "Processed has less than 2 points." );
    return;
    }

  this->m_NumberOfPoints = processedInput->GetPoints().size();
  this->m_TubeId = processedInput->GetId();
  if( ipm )
    {
    this->m_InflectionPoints.SetSize( this->m_NumberOfPoints );
    }
  if( m_MeasureFlag & BITMASK_CURVATURE_METRICS )
    {
    this->m_CurvatureScalar.resize( this->m_NumberOfPoints );
    }
  if( cvm )
    {
    this->m_CurvatureVector.resize( this->m_NumberOfPoints );
    }

  for( size_t index = 0; index < this->m_NumberOfPoints; ++index )
    {
    SOVectorType currentPoint = processedInput->GetPoint( index )->
      GetPositionInObjectSpace().GetVectorFromOrigin();

    // General variables
    bool nextPointAvailable = ( index < this->m_NumberOfPoints - 1 );
    SOVectorType nextPoint( 0.0 );
    if( nextPointAvailable )
      {
      nextPoint = processedInput->GetPoint( index + 1 )->
        GetPositionInObjectSpace().GetVectorFromOrigin();
      }
    bool previousPointAvailable = ( index > 0 );
    SOVectorType previousPoint( 0.0 );
    if( previousPointAvailable )
      {
      previousPoint = processedInput->GetPoint( index - 1 )->
        GetPositionInObjectSpace().GetVectorFromOrigin();
      }
    // t1 and t2, used both in icm and soam
    SOVectorType t1( 0.0 ), t2( 0.0 );
    if( previousPointAvailable && nextPointAvailable )
      {
      t1 = currentPoint - previousPoint;
      t2 = nextPoint - currentPoint;
      }

    bool nPlus2PointAvailable = ( index < this->m_NumberOfPoints - 2 );
    SOVectorType nPlus2Point( 0.0 );
    if( nPlus2PointAvailable )
      {
      nPlus2Point = processedInput->GetPoint( index + 2 )->
        GetPositionInObjectSpace().GetVectorFromOrigin();
      }

    //
    // DM Computations
    if( index == 0 )
      {
      start = currentPoint;
      currentPoint = start;
      }
    if( index == this->m_NumberOfPoints - 1 )
      {
      end = currentPoint;
      }

    if( ( dm || icm || ipm || soam || plm ) && nextPointAvailable )
      {
      pathLength += ( nextPoint - currentPoint ).GetNorm();
      }

    //
    // ICM Computations
    double inflectionValue = 0.0;
    if( ( icm || ipm ) && previousPointAvailable && nextPointAvailable )
      {
      // Compute velocity and acceleration
      SOVectorType v = nextPoint - previousPoint;
      SOVectorType a = t2 - t1;

      // Compute the Frenet frame
      // 1 - T = v / |v|
      T = v;
      SafeNormalize( T );

      // 2 - N = v x a x v / | v x a x v |
      bool canCheckForinflection = a.GetNorm() > this->m_EpsilonForZero;
      if( canCheckForinflection )
        {
        N = CrossProduct( v, a );
        N = CrossProduct( v, N );
        SafeNormalize( N );
        vectorBIsValid = true;
        }
      else if( vectorBIsValid ) // 2nd chance
        {
        // Acceleration can be null when the curve approximates a straight
        // line ( sin around pi for example ). Unfortunately that could
        // happen when the curve is crossing the straight line and the
        // inflection would be missed...
        // This assumes that no pure torsion along the N vector happened.
        // Note that this is only valid is B was already computed at least
        // once.
        N = CrossProduct( B, T );
        SafeNormalize( N );
        canCheckForinflection = true;
        }

      // 3 - B = T x N ( in case of null acceleration. See above )
      B = CrossProduct( T, N );
      SafeNormalize( B );

      if( canCheckForinflection )
        {
        // Check for inflection
        SOVectorType deltaN = N - previousN;

        inflectionValue = deltaN * deltaN;
        if( inflectionValue > 1.0 + this->m_EpsilonForZero )
          {
          inflectionCount += 1;
          }
        }

      previousN = N;
      }

    // Set the inflection value for this point
    if( ipm )
      {
      this->m_InflectionPoints.SetElement( index, inflectionValue );
      }

    if( ( soam || sot ) && previousPointAvailable && nextPointAvailable &&
      nPlus2PointAvailable )
      {
      // Compute in-plane angle
      SOVectorType normT1 = t1;
      SOVectorType normT2 = t2;

      double inPlaneAngle = 0.0;
      if( SafeNormalize( normT1 ) > this->m_EpsilonForZero
        && SafeNormalize( normT2 ) > this->m_EpsilonForZero )
        {
        inPlaneAngle = SafeAcos( normT1 * normT2 );
        }

      // Compute torsionnal angle
      SOVectorType t3 = nPlus2Point - nextPoint;
      SOVectorType t1t2Cross = CrossProduct( t1, t2 );
      SOVectorType t2t3Cross = CrossProduct( t2, t3 );

      double torsionAngle = 0.0;
      if( SafeNormalize( t1t2Cross ) > this->m_EpsilonForZero
        && SafeNormalize( t2t3Cross ) > this->m_EpsilonForZero )
        {
        double t1t2t2t3Dot = t1t2Cross * t2t3Cross;
        // It is confusing to include points with torsional angles of
        // 180 degrees when analyzing a planar curve, so we set the angle
        // to 0.
        if( t1t2t2t3Dot < -1 + this->m_EpsilonForZero )
          {
          t1t2t2t3Dot = 1;
          }
        torsionAngle = SafeAcos( t1t2t2t3Dot );
        }

      // Finally add the angle to the sum
      sumOfAngles += sqrt( inPlaneAngle*inPlaneAngle );
      sumOfTorsion += sqrt( torsionAngle*torsionAngle );
      }

    if( arm )
      {
      // Average radius computation
      sumOfRadius += processedInput->GetPoints()[index].GetRadiusInObjectSpace();
      }

    // Metrics that require curvature computation
    if( m_MeasureFlag & BITMASK_CURVATURE_METRICS )
      {
      SOVectorType dg;
      SOVectorType d2g;
      SOVectorType curvatureVector;
      double dtl;
      double dtr;
      double curvatureScalar;

      // Calculate first derivative
      dg = t2;
      SafeNormalize( dg );
      // Calculate second derivative
      d2g = t2-t1;
      dtl = t1.GetNorm();
      dtr = t2.GetNorm();
      if( dtl > this->m_EpsilonForSpacing &&
        dtr > this->m_EpsilonForSpacing )
        {
        d2g /= ( dtl * dtr );
        }

      // Calculate curvature vector
      // No need to divide by cube of dg norm, since dg is normalized
      curvatureVector = CrossProduct( dg, d2g );
      if( cvm )
        {
        this->m_CurvatureVector.at( index ) = curvatureVector;
        }

      // Calculate curvature scalar
      curvatureScalar = curvatureVector.GetNorm();
      // If any curvature metric is asked, we want the curvature scalars
      // to be computed, so the "csm" flag is not necessary
      this->m_CurvatureScalar.at( index ) = curvatureScalar;

      if( tcm )
        {
        totalCurvature += curvatureScalar;
        }
      if( tscm )
        {
        totalSquaredCurvature += curvatureScalar*curvatureScalar;
        }
      }
    }

  // Metrics final calculation

  if( plm )
    {
    this->m_PathLengthMetric = pathLength;
    }

  if( dm || icm || clm )
    {
    double straightLineLength = ( start - end ).GetNorm();

    if( straightLineLength > 0.0 )
      {
      if( clm )
        {
        this->m_ChordLengthMetric = straightLineLength;
        }
      if( pathLength / straightLineLength < 1.0 )
        {
        itkExceptionMacro( << "Error while computing the distance metric."
          << "DM ( =" << pathLength / straightLineLength << " ) < 1.0" );
        }
      else
        {
        this->m_DistanceMetric = pathLength / straightLineLength;
        }
      }
    }

  if( icm )
    {
    this->m_InflectionCountMetric = inflectionCount *
      this->m_DistanceMetric;
    }

  if( soam || sot )
    {
    if( pathLength > 0.0 )
      {
      this->m_SumOfAnglesMetric = sumOfAngles / pathLength;
      this->m_SumOfTorsionMetric = sumOfTorsion / pathLength;
      }
    else
      {
      itkExceptionMacro( <<"Cannot compute angle metrics, path length ( ="
        << pathLength << " ) <= 0.0" );
      }
    }

  if( arm )
    {
    this->m_AverageRadiusMetric = sumOfRadius/this->m_NumberOfPoints;
    }
  if( tcm )
    {
    this->m_TotalCurvatureMetric = totalCurvature;
    }
  if( tscm )
    {
    this->m_TotalSquaredCurvatureMetric = totalSquaredCurvature;
    }

  if( p95m || chm || ic1m || ic2m )
    {
    // Histogram computation
    typedef itk::Vector<double, 1> MeasurementVectorType;
    typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;

    typedef itk::Statistics::Histogram< double,
      itk::Statistics::DenseFrequencyContainer2 > HistogramType;

    SampleType::Pointer sample = SampleType::New();
    for( size_t i = 0; i < this->m_CurvatureScalar.size(); ++i )
      {
      MeasurementVectorType mv;
      mv[0] = this->m_CurvatureScalar[i];
      sample->PushBack( mv );
      }

    typedef itk::Statistics::SampleToHistogramFilter
      <SampleType, HistogramType> SampleToHistogramFilterType;
    SampleToHistogramFilterType::Pointer sampleToHistogramFilter =
      SampleToHistogramFilterType::New();
    sampleToHistogramFilter->SetInput( sample );

    SampleToHistogramFilterType::HistogramSizeType histogramSize( 1 );
    histogramSize.Fill( 10 );
    sampleToHistogramFilter->SetHistogramSize( histogramSize );

    sampleToHistogramFilter->Update();

    const HistogramType* histogram = sampleToHistogramFilter->GetOutput();

    if( p95m || ic2m )
      {
      this->m_Percentile95Metric = histogram->Quantile( 0, 0.95 );
      }

    if( ic1m || ic2m )
      {
      // Compute Inflection Count with 2 different methods

      // 1st Method: Every blob of curvature curve > lambda/sigma
      // is considered as an inflection.
      int inflectionCount1 = 0;
      double threshold1 = this->m_Lambda/this->m_SmoothingScale;

      // 2nd Method: Every blob of curvature curve > 95 percentile
      // is considered as an inflection.
      int inflectionCount2 = 0;
      double threshold2 = this->m_Percentile95Metric;

      for( size_t i = 1; i <this->m_CurvatureScalar.size(); ++i )
        {
        if( ic1m && ( this->m_CurvatureScalar[i-1]-threshold1 )*
            ( this->m_CurvatureScalar[i]-threshold1 ) < 0 )
          {
          ++inflectionCount1;
          }
        if( ic2m && ( this->m_CurvatureScalar[i-1]-threshold2 )*
            ( this->m_CurvatureScalar[i]-threshold2 ) < 0 )
          {
          ++inflectionCount2;
          }
        }
      if( ic1m )
        {
        this->m_InflectionCount1Metric = inflectionCount1 / 2;
        }
      if( ic2m )
        {
        this->m_InflectionCount2Metric = inflectionCount2 / 2;
        }
      }

    if( chm )
      {
      // Compute 2nd histogram: for features, not for percentiles
      // This one has a range specified
      SampleToHistogramFilterType::Pointer sampleToHistogramFilter2 =
        SampleToHistogramFilterType::New();
      sampleToHistogramFilter2->SetInput( sample );

      // Set Number of bins
      SampleToHistogramFilterType::HistogramSizeType histogramSize2( 1 );
      histogramSize2.Fill( this->m_NumberOfBins );
      sampleToHistogramFilter2->SetHistogramSize( histogramSize2 );

      // Set Minimum of the histogram
      HistogramType::MeasurementVectorType minimum;
      minimum.SetSize( 1 );
      minimum[0] = this->m_HistogramMin;
      sampleToHistogramFilter2->SetHistogramBinMinimum( minimum );

      // Set Maximum of the histogram
      HistogramType::MeasurementVectorType maximum;
      maximum.SetSize( 1 );
      maximum[0] = this->m_HistogramMax;
      sampleToHistogramFilter2->SetHistogramBinMaximum( maximum );

      // Don't set automatic range
      sampleToHistogramFilter2->SetAutoMinimumMaximum( false );

      // Compute the histogram
      sampleToHistogramFilter2->Update();
      const HistogramType* histogramForFeatures =
        sampleToHistogramFilter2->GetOutput();

      if( histogramForFeatures->Size() != this->m_NumberOfBins )
        {
        itkExceptionMacro( << "The histogram is not of expected size:"
          << std::endl << "Expected = " << m_NumberOfBins << std::endl
          << "Real size = " << histogramForFeatures->Size() << std::endl );
        }
      else
        {
        // Add bin values to the metric array
        for( unsigned int i=0; i < this->m_NumberOfBins; ++i )
          {
          this->m_CurvatureHistogramMetrics.push_back(
            histogramForFeatures->GetFrequency( i ) );
          }
        }
      }
    }

  output->CopyInformation( processedInput );
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTortuositySpatialObjectFilter_hxx )
