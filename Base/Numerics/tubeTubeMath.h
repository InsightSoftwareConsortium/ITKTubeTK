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
#ifndef __tubeTubeMath_h
#define __tubeTubeMath_h

#include "tubeMacro.h"

#include <vector>

namespace tube
{

template< class TTubePoint >
void
ComputeNormalsFromTangent( TTubePoint & tubePoint,
  const typename TTubePoint::VectorType & prevT );

template< class TTube >
bool
ComputeTubeTangentsAndNormals( typename TTube::Pointer & tube );

template< class TTubePoint >
bool
ComputeVectorTangentsAndNormals( std::vector< TTubePoint > & tube );

enum SmoothTubeFunctionEnum { SMOOTH_TUBE_USING_INDEX_AVERAGE,
  SMOOTH_TUBE_USING_INDEX_GAUSSIAN, SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN };

/** Smooth a tube
 * The parameter h has different meanings when using different smoothing
 * functions:
 *
 * smoothFunction = SMOOTH_TUBE_USING_INDEX_AVERAGE:
 *    h is half of the window size
 * smoothFunction = SMOOTH_TUBE_USING_INDEX_GAUSSIAN:
 *    h is the gaussian's standard deviation
 * smoothFunction = SMOOTH_TUBE_USING_DISTANCE_GAUSSIAN:
 *    h is the gaussian's standard deviation
 */
template< class TTube >
typename TTube::Pointer
SmoothTube( const typename TTube::Pointer & tube, double h = 2,
  SmoothTubeFunctionEnum smoothFunction = SMOOTH_TUBE_USING_INDEX_AVERAGE );

template< class TTube >
int
RemoveDuplicateTubePoints( typename TTube::Pointer & tube );

template< class TTube >
typename TTube::Pointer
SubsampleTube( const typename TTube::Pointer & tube, int N = 2 );

template< class TTube >
double
ComputeTubeLength( const typename TTube::Pointer & tube );

} // End namespace tube

#include "tubeTubeMath.hxx"

#endif // End !defined( __tubeTubeMath_h )
