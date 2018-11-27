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

#include "itktubeTubeParametricExponentialResolutionWeightFunction.h"

#include <itkTubeSpatialObjectPoint.h>

#include <fstream>

template< class TWeightFunction, class TTubePoint >
void
writeElement( std::ofstream & outputFile,
  TWeightFunction & weightFunction,
  const TTubePoint & tubePoint,
  float alpha,
  float delta,
  char delim=',' )
{
  weightFunction.SetAlpha( alpha );
  weightFunction.SetDelta( delta );
  float pointWeight = weightFunction( tubePoint );
  outputFile << pointWeight << delim;
}

int itktubeTubeParametricExponentialResolutionWeightFunctionTest( int argc, char * argv[] )
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
  typedef itk::TubeSpatialObjectPoint< Dimension >  TubePointType;
  typedef std::vector< TubePointType >              TubePointContainerType;

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

  // Write the tube point weights to file.
  std::ofstream outputFile( outputCSV );
  if( !outputFile.is_open() )
    {
    std::cerr << "Cannot open output file." << std::endl;
    return EXIT_FAILURE;
    }

  typedef itk::tube::Function::TubeParametricExponentialResolutionWeightFunction< TubePointType, float >
    WeightFunctionType;
  WeightFunctionType weightFunction;

  weightFunction.SetAlpha( 2.0f );
  float alpha = weightFunction.GetAlpha();
  if( alpha != 2.0f )
    {
    std::cerr << "Failure in GetAlpha." << std::endl;
    return EXIT_FAILURE;
    }

  weightFunction.SetDelta( 2.0f );
  float delta = weightFunction.GetDelta();
  if( delta != 2.0f )
    {
    std::cerr << "Failure in GetDelta." << std::endl;
    return EXIT_FAILURE;
    }

  unsigned int numberOfAlphas = 5;
  float alphas[] = { 5.0, 2.0, 1.0, 0.5, 0.0 };
  unsigned int numberOfDeltas = 5;
  float deltas[] = { 0.7, 0.5, 0.0, -0.5, -0.7 };

  for( unsigned int ii = 0; ii < numberOfPoints; ++ii )
    {
    for( unsigned int jj = 0; jj < numberOfAlphas - 1; ++jj )
      {
      for( unsigned int kk = 0; kk < numberOfDeltas; ++kk )
        {
        writeElement( outputFile, weightFunction,
          tubePointContainer[ii], alphas[jj], deltas[kk] );
        }
      }
    for( unsigned int kk = 0; kk < numberOfDeltas - 1; ++kk )
      {
      writeElement( outputFile, weightFunction,
        tubePointContainer[ii], alphas[numberOfAlphas-1], deltas[kk] );
      }
    writeElement( outputFile, weightFunction,
      tubePointContainer[ii], alphas[numberOfAlphas-1], deltas[numberOfDeltas-1], '\n' );
    }

  outputFile.flush();
  outputFile.close();

  return EXIT_SUCCESS;
}
