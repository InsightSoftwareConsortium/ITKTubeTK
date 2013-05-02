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

#ifndef __itkTubeExponentialResolutionWeightFunction_h
#define __itkTubeExponentialResolutionWeightFunction_h

#include "vnl/vnl_math.h"

namespace itk
{
namespace Function
{

/** \class TubeExponentialResolutionWeightFunction
 *
 * \brief Weight tube points exponentially by their radius.
 *
 * \f$ w_i = \frac{2}{1 + e^{-2 r_i}} \f$
 *
 * As in Eqn. 2. Alyward, S. Weeks, S. and Bullitt, E.  Analysis of the
 * Parameter Space of a Metric for Registering 3D Vascular Images.  MICCAI, 2001.
 *
 * \sa TubeExponentialResolutionWeightFunction
 * \sa TubeExponentialWithBoundsResolutionWeightFunction
 */
template< class TTubePoint, class TOperatorValue=double >
class TubeExponentialResolutionWeightFunction
{
public:
  typedef TOperatorValue OperatorValueType;
  typedef TTubePoint     TubePointType;

  TubeExponentialResolutionWeightFunction( void )
    {}
  ~TubeExponentialResolutionWeightFunction( void )
    {}

  inline OperatorValueType operator()( const TubePointType & tubePoint )
    {
    const OperatorValueType radius = tubePoint.GetRadius();
    return static_cast< OperatorValueType >( 2.0 /
      (1.0 + vcl_exp( -2.0 * radius ) ));
    }

}; // End class TubeExponentialResolutionWeightFunction

} // End namespace Function

} // End namespace itk

#endif // End !defined(__itkTubeExponentialResolutionWeightFunction_h)
