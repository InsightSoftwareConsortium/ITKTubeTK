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

#ifndef __itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction_h
#define __itktubeTubeParametricExponentialWithBoundsResolutionWeightFunction_h

#include "itktubeTubeParametricExponentialResolutionWeightFunction.h"

#include <itkNumericTraits.h>

namespace itk
{

namespace tube
{

namespace Function
{

/** \class TubeParametricExponentialWithBoundsResolutionWeightFunction
 *
 * \brief Weight tube points exponentially by their radius if within bounds.
 *
 * This is similar to TubeParametricExponentialResolutionWeightFunction except that values
 * outside the UpperBound or LowerBound are given a weight of zero.
 */
template< class TTubePoint, class TOperatorValue = double >
class TubeParametricExponentialWithBoundsResolutionWeightFunction
  : public TubeParametricExponentialResolutionWeightFunction< TTubePoint,
                                                              TOperatorValue >
{
public:
  typedef TubeParametricExponentialWithBoundsResolutionWeightFunction
    Self;
  typedef TubeParametricExponentialResolutionWeightFunction< TTubePoint, TOperatorValue >
    Superclass;

  typedef typename Superclass::OperatorValueType OperatorValueType;
  typedef typename Superclass::TubePointType     TubePointType;

  TubeParametricExponentialWithBoundsResolutionWeightFunction( void )
    : m_LowerBound( NumericTraits< OperatorValueType >::min() ),
      m_UpperBound( NumericTraits< OperatorValueType >::max() )
    {}
  ~TubeParametricExponentialWithBoundsResolutionWeightFunction( void )
    {}

  inline OperatorValueType operator()( const TubePointType & tubePoint )
    {
    const OperatorValueType radius = tubePoint.GetRadius();
    if( radius < this->m_LowerBound || radius > this->m_UpperBound )
      {
      return NumericTraits< OperatorValueType >::Zero;
      }
    return static_cast< OperatorValueType >( 1.0 /
      ( 1.0 + this->GetDelta() * ( std::exp( -this->GetAlpha() * radius ) - 1.0 ) ) );
    }

  void SetLowerBound( const OperatorValueType lowerBound )
    {
    this->m_LowerBound = lowerBound;
    }

  OperatorValueType GetLowerBound( void ) const
    {
    return this->m_LowerBound;
    }

  void SetUpperBound( const OperatorValueType upperBound )
    {
    this->m_UpperBound = upperBound;
    }

  OperatorValueType GetUpperBound( void ) const
    {
    return this->m_UpperBound;
    }

  bool operator!=( const Self & other ) const
    {
    if( Superclass::operator!=( other ) ||
        this->m_LowerBound != other.m_LowerBound ||
        this->m_UpperBound != other.m_UpperBound )
      {
      return true;
      }
    return false;
    }

  bool operator==( const Self & other ) const
    {
    return !( *this != other );
    }


private:
  OperatorValueType m_LowerBound;
  OperatorValueType m_UpperBound;

}; // End class TubeParametricExponentialWithBoundsResolutionWeightFunction

} // End namespace Function

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktub...
