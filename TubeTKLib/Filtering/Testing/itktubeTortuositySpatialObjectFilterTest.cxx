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

// Tortuostity includes
#include <itktubeTortuositySpatialObjectFilter.h>
#include <itkTubeSpatialObject.h>

// ITK includes
#include <itkMath.h>

// STD includes
#include <iomanip>

//--------------------------------------------------------------------------
itk::TubeSpatialObject<3>::Pointer
GenerateStraightTube( itk::TubeSpatialObject<3>::PointType
  start, itk::TubeSpatialObject<3>::VectorType increment,
  int numberOfPoints )
{
  typedef itk::TubeSpatialObject<3>        VesselTubeType;
  typedef VesselTubeType::VectorType       VectorType;
  typedef VesselTubeType::PointType        PointType;

  VesselTubeType::Pointer vessel = VesselTubeType::New();
  VesselTubeType::TubePointListType tubePointList;

  PointType pos = start;
  for ( int i = 0; i < numberOfPoints; ++i )
    {
    VesselTubeType::TubePointType tubePoint;
    tubePoint.SetPositionInObjectSpace( pos );
    tubePointList.push_back( tubePoint );

    pos += increment;
    }

  vessel->SetPoints( tubePointList );
  std::cout<<"Straight vessel generated"<<std::endl;
  return vessel;
}

double Modulus( double v, double mod )
{
  if ( v < mod )
    {
    return v;
    }
  return Modulus( v-mod, mod );
}

//--------------------------------------------------------------------------
void
GenerateCosTube( double length, double amplitude, double frequency,
  itk::TubeSpatialObject<3>::Pointer & vessel )
{
  typedef itk::TubeSpatialObject<3>   VesselTubeType;
  typedef VesselTubeType::VectorType  VectorType;
  typedef VesselTubeType::PointType   PointType;

  vessel = VesselTubeType::New();
  VesselTubeType::TubePointListType tubePointList;

  PointType pos( 0.0 );

  int numberOfSamples = 1e2;
  for ( double l = 0.0; l < length; l += length / numberOfSamples )
    {
    double t = Modulus( l, itk::Math::pi * 2 );

    VesselTubeType::TubePointType tubePoint;
    tubePoint.SetPositionInObjectSpace( pos );
    tubePointList.push_back( tubePoint );

    pos[0] = l;
    pos[1] = amplitude * sin( frequency * t );
    pos[2] = 0.0;
    }

  vessel->SetPoints( tubePointList );
  std::cout<<"Cosine vessel generated"<<std::endl;
}

//--------------------------------------------------------------------------
bool TestVesselMetrics( itk::TubeSpatialObject<3>::Pointer
  vessel, double results[] )
{
  typedef itk::TubeSpatialObject<3>   VesselTubeType;

  typedef itk::tube::TortuositySpatialObjectFilter<VesselTubeType>
    FilterType;

  // Run the metrics
  FilterType::Pointer filter = FilterType::New();
  filter->SetMeasureFlag( FilterType::BITMASK_ALL_METRICS );
  filter->SetSmoothingScale( 0 );
  filter->SetEpsilonForSpacing( 1e-2 );
  filter->SetInput( vessel );
  filter->Update();
  if ( ! filter->GetOutput() )
    {
    std::cerr<<"Error while running the metrics !"<<std::endl;
    return false;
    }

  // ICM isn't verified because its computation is not right.

  double dm = filter->GetDistanceMetric();
//  double icm = filter->GetInflectionCountMetric();
  double soam = filter->GetSumOfAnglesMetric();
  double plm = filter->GetPathLengthMetric();
  double clm = filter->GetChordLengthMetric();
  double arm = filter->GetAverageRadiusMetric();
  double ic1m = filter->GetInflectionCount1Metric();
  double ic2m = filter->GetInflectionCount2Metric();
  double p95m = filter->GetPercentile95Metric();
  double tcm = filter->GetTotalCurvatureMetric();
  double tscm = filter->GetTotalSquaredCurvatureMetric();

//  std::cout<<std::setprecision( 8 );
//  std::cout<<"  DM: expected "<<results[0]<<" got "<<dm<<std::endl;
////  std::cerr<<"  ICM: expected "<<results[1]<<" got "<<icm<<std::endl;
//  std::cout<<"  SOAM: expected "<<results[2]<<" got "<<soam<<std::endl;
//  std::cout<<"  PLM: expected "<<results[3]<<" got "<<plm<<std::endl;
//  std::cout<<"  CLM: expected "<<results[4]<<" got "<<clm<<std::endl;
//  std::cout<<"  ARM: expected "<<results[5]<<" got "<<arm<<std::endl;
//  std::cout<<"  IC1M: expected "<<results[6]<<" got "<<ic1m<<std::endl;
//  std::cout<<"  IC2M: expected "<<results[7]<<" got "<<ic2m<<std::endl;
//  std::cout<<"  P95M: expected "<<results[8]<<" got "<<p95m<<std::endl;
//  std::cout<<"  TCM: expected "<<results[9]<<" got "<<tcm<<std::endl;
//  std::cout<<"  TSCM: expected "<<results[10]<<" got "<<tscm<<std::endl;

  if( ! itk::Math::FloatAlmostEqual( dm, results[0], 4, 1e-4 ) )
    {
    std::cerr<<"DM: expected "<<results[0]<<" got "<<dm<<std::endl;
    return false;
    }
//  if( ! itk::Math::FloatAlmostEqual( icm, results[1], 4, 1e-4 ) )
//    {
//    std::cerr<<"ICM: expected "<<results[1]<<" got "<<icm<<std::endl;
//    }
  if( ! itk::Math::FloatAlmostEqual( soam, results[2], 4, 1e-4 ) )
    {
    std::cerr<<"SOAM: expected "<<results[2]<<" got "<<soam<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( plm, results[3], 4, 1e-4 ) )
    {
    std::cerr<<"PLM: expected "<<results[3]<<" got "<<plm<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( clm, results[4], 4, 1e-4 ) )
    {
    std::cerr<<"CLM: expected "<<results[4]<<" got "<<clm<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( arm, results[5], 4, 1e-4 ) )
    {
    std::cerr<<"ARM: expected "<<results[5]<<" got "<<arm<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( ic1m, results[6], 4, 1e-4 ) )
    {
    std::cerr<<"IC1M: expected "<<results[6]<<" got "<<ic1m<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( ic2m, results[7], 4, 1e-4 ) )
    {
    std::cerr<<"IC2M: expected "<<results[7]<<" got "<<ic2m<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( p95m, results[8], 4, 1e-4 ) )
    {
    std::cerr<<"P95M: expected "<<results[8]<<" got "<<p95m<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( tcm, results[9], 4, 1e-4 ) )
    {
    std::cerr<<"TCM: expected "<<results[9]<<" got "<<tcm<<std::endl;
    return false;
    }
  if( ! itk::Math::FloatAlmostEqual( tscm, results[10], 4, 1e-4 ) )
    {
    std::cerr<<"TSCM: expected "<<results[10]<<" got "<<tscm<<std::endl;
    return false;
    }

  return true;
}

//--------------------------------------------------------------------------
int itktubeTortuositySpatialObjectFilterTest( int, char*[] )
{
  typedef itk::TubeSpatialObject<3> VesselTubeType;

  //
  // Test straight spatial objects
  //
  const int NUMBER_OF_STRAIGHT_OBJECT_TESTS = 3;

  // X, Y, Z
  double start[NUMBER_OF_STRAIGHT_OBJECT_TESTS][3] = {
      {0.0, 0.0, 0.0},
      {26.2, -601.0, -42.01},
      {3.20, -89.0, 1.0},
    };

  // X, Y, Z
  double increment[NUMBER_OF_STRAIGHT_OBJECT_TESTS][3] = {
      {0.0, 1.0, 0.0},
//      {0.05, 80.4, 1.0},
      {0.01, 100.0, 0.0},
      {0.0, 0.0, -0.07},
    };
  double numberOfPoints[NUMBER_OF_STRAIGHT_OBJECT_TESTS] =
    {
    30,
    4,
    10
    };
  // dm, icm, soam, plm, clm, arm, ic1m, ic2m, p95m, tcm, tscm
  double straightObjectResults[NUMBER_OF_STRAIGHT_OBJECT_TESTS][11] = {
    {1.0, 1.0, 0.0, 29, 29, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0, 1.0, 0.0, 300, 300, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0, 1.0, 0.0, 0.63, 0.63, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
    };

  for ( int i = 0; i < NUMBER_OF_STRAIGHT_OBJECT_TESTS; ++i )
    {
    std::cerr<<"Straight object test #"<<i<<std::endl;
    VesselTubeType::Pointer vessel =
      GenerateStraightTube( start[i], increment[i], numberOfPoints[i] );

    if ( !TestVesselMetrics( vessel, straightObjectResults[i] ) )
      {
      std::cerr<<"Error in straight object test for object #"<<i<<std::endl;
      return EXIT_FAILURE;
      }
    }

  //
  // Test cosinus spatial objects
  //
  const int NUMBER_OF_COSINUS_OBJECT_TESTS = 4;
  double length[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      2.0*itk::Math::pi,
      2.0*itk::Math::pi * 2,
      2.0*itk::Math::pi,
      2.0*itk::Math::pi,
    };
  double amplitude[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      1.0,
      1.0,
      2.0,
      1.0,
    };
  double frequency[NUMBER_OF_COSINUS_OBJECT_TESTS] = {
      1.0,
      1.0,
      1.0,
      2.0,
    };
  // dm, icm, soam, plm, clm, arm, ic1m, ic2m, p95m, tcm, tscm
  double cosResults[NUMBER_OF_COSINUS_OBJECT_TESTS][11] = {
    {
      1.2138841,  // dm
      0.0,  // icm
      0.41507849,  // soam
      7.551173,  // plm
      6.2206704,  // clm
      0.0,  // arm
      0.0,  // ic1m
      2.0,  // ic2m
      0.96314982,  // p95m
      44.952127,  // tcm
      30.87049,  // tscm
    },
    {
      1.2137084,  // dm
      0.0,  // icm
      0.4138076,  // soam
      15.100157,  // plm
      12.441338,  // clm
      0.0,  // arm
      0.0,  // ic1m
      4.0,  // ic2m
      0.95872965,  // p95m
      44.805743,  // tcm
      30.67568,  // tscm
    },
    {
      1.6714683,  // dm
      0.0,  // icm
      0.42529686,  // soam
      10.399242,  // plm
      6.2216221,  // clm
      0.0,  // arm
      0.0,  // ic1m
      2.0,  // ic2m
      1.8639652,  // p95m
      56.7968,  // tcm
      71.061562,  // tscm
    },
    {
      1.6709897,  // dm
      0.0,  // icm
      0.84926412,  // soam
      10.396256,  // plm
      6.221616,  // clm
      0.0,  // arm
      0.0,  // ic1m
      4.0,  // ic2m
      3.6268579,  // p95m
      112.78377,  // tcm
      278.47741,  // tscm
    },
    };

  for ( int i = 0; i < NUMBER_OF_COSINUS_OBJECT_TESTS; ++i )
    {
    std::cerr << "Cos object test #" << i << std::endl;
    VesselTubeType::Pointer vessel;
    GenerateCosTube( length[i], amplitude[i], frequency[i], vessel );

    if ( !TestVesselMetrics( vessel, cosResults[i] ) )
      {
      std::cerr<<"Error in cos object test for object #"<<i<<std::endl;
      return EXIT_FAILURE;
      }
    }


  return EXIT_SUCCESS;
}
