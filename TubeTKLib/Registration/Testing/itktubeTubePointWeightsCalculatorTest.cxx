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

#include "itktubeTubeExponentialResolutionWeightFunction.h"
#include "itktubeTubePointWeightsCalculator.h"

#include <itkOptimizerParameters.h>
#include <itkTubeSpatialObject.h>
#include <itkGroupSpatialObject.h>

#include <fstream>

int itktubeTubePointWeightsCalculatorTest( int argc, char * argv[] )
{
  if( argc < 2 )
    {
    std::cerr << "Usage: "
      << argv[0]
      << " output.csv"
      << std::endl;
    return EXIT_FAILURE;
    }
  const char * outputCSV = argv[1];

  enum { Dimension = 2 };
  typedef itk::GroupSpatialObject< Dimension >      GroupType;
  typedef itk::TubeSpatialObject< Dimension >       TubeType;
  typedef TubeType::TubePointType                   TubePointType;
  typedef TubeType::PointListType                   TubePointContainerType;

  // Create tube points with increasing radius.
  const unsigned int numberOfPoints = 100;
  TubePointContainerType tubePointContainer( numberOfPoints );
  const double startingRadius = 0.00;
  const double radiusIncrement = 0.01;
  double radius = startingRadius;
  for( unsigned int ii = 0; ii < numberOfPoints; ++ii )
    {
    tubePointContainer[ii].SetRadius( radius );
    radius += radiusIncrement;
    }
  TubeType::Pointer tubes = TubeType::New();
  tubes->SetPoints( tubePointContainer );

  GroupType::Pointer group = GroupType::New();
  GroupType::ChildrenListType children;
  children.push_back( tubes.GetPointer() );
  group->SetChildren( children );

  // Write the tube point weights to file.
  std::ofstream outputFile( outputCSV );
  if( !outputFile.is_open() )
    {
    std::cerr << "Cannot open output file." << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::tube::Function::TubeExponentialResolutionWeightFunction<
    TubePointType, float > WeightFunctionType;
  WeightFunctionType::Pointer weightFunction = WeightFunctionType::New();

  typedef itk::OptimizerParameters< float > PointWeightsType;

  typedef itk::tube::TubePointWeightsCalculator< Dimension,
    TubeType, WeightFunctionType,
    PointWeightsType > PointWeightsCalculatorType;
  PointWeightsCalculatorType::Pointer resolutionWeightsCalculator
    = PointWeightsCalculatorType::New();
  resolutionWeightsCalculator->SetTubeTreeSpatialObject( group );
  resolutionWeightsCalculator->SetPointWeightFunction( weightFunction );
  resolutionWeightsCalculator->Compute();
  const PointWeightsType & resolutionWeights
    = resolutionWeightsCalculator->GetPointWeights();

  for( unsigned int ii = 0; ii < numberOfPoints; ++ii )
    {
    const float pointWeight = resolutionWeights[ii];
    outputFile << pointWeight << '\n';
    }

  outputFile.flush();
  outputFile.close();

  // test PrintSelf
  std::cout << resolutionWeightsCalculator << std::endl;

  return EXIT_SUCCESS;
}
