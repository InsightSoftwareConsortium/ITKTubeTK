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
#ifndef itkTubeExponentialResolutionWeightFunction_h
#define itkTubeExponentialResolutionWeightFunction_h

#include "vnl/vnl_math.h"

namespace itk
{
namespace Function
{

/**
 * \class TubeExponentialResolutionWeightFunction
 *
 * \brief Weight tube points exponentially by their radius.
 *
 * \f$ w_i = \frac{1}{1 + \delta (e^{-\alpha r_i} - 1)} \f$
 *
 * where \f$\delta\ \in (-1, 1)f$ controls the amount of weighting.  A
 * positive Delta weighs large tubes higher than small tubes, and a
 * negative Delta weighs small tubes higher than large tubes.
 * \f$\alpha \in [0, \infty)\f$ controls transition in weights.
 * \f$r\f$ is the tube radius at a point.
 *
 * */
template< class TTubePoint, class TOperatorValue=double >
class TubeExponentialResolutionWeightFunction
{
public:
  typedef TOperatorValue OperatorValueType;
  typedef TTubePoint     TubePointType;

  TubeExponentialResolutionWeightFunction():
    m_Alpha( 2.0 )
    {}
  ~TubeExponentialResolutionWeightFunction()
    {}

  inline OperatorValueType operator()( const TubePointType & tubePoint )
    {
    const OperatorValueType radius = tubePoint.GetRadius();
    return static_cast< OperatorValueType >( 1.0 /
      (1.0 + this->m_Delta * (vcl_exp( -this->m_Alpha * radius ) - 1.0)));
    }

  void SetDelta( const OperatorValueType delta )
    {
    this->m_Delta = delta;
    }

  OperatorValueType GetDelta( void ) const
    {
    return this->m_Delta;
    }

  void SetAlpha( const OperatorValueType alpha )
    {
    this->m_Alpha = alpha;
    }

  OperatorValueType GetAlpha( void ) const
    {
    return this->m_Alpha;
    }

private:
  OperatorValueType m_Delta;
  OperatorValueType m_Alpha;
};

} // end namespace Function
} // end namespace itk

#endif
