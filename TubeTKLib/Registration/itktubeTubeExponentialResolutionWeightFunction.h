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

#ifndef __itktubeTubeExponentialResolutionWeightFunction_h
#define __itktubeTubeExponentialResolutionWeightFunction_h

#include <vnl/vnl_math.h>
#include <itkFunctionBase.h>

namespace itk
{

namespace tube
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
 * Parameter Space of a Metric for Registering 3D Vascular Images.
 * MICCAI, 2001.
 *
 * \sa TubeExponentialResolutionWeightFunction
 * \sa TubeExponentialWithBoundsResolutionWeightFunction
 */
template< class TTubePoint, class TWeight = double >
class TubeExponentialResolutionWeightFunction:
  public FunctionBase< TTubePoint, TWeight >
{
public:
  /** Standard class typedefs. */
  typedef TubeExponentialResolutionWeightFunction< TTubePoint, TWeight >
                                               Self;
  typedef FunctionBase< TTubePoint, TWeight >  Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( TubeExponentialResolutionWeightFunction, FunctionBase );

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  typedef TWeight    WeightType;
  typedef TTubePoint TubePointType;

  WeightType Evaluate( const TubePointType & tubePoint ) const
    {
    const WeightType radius = tubePoint.GetRadiusInObjectSpace();
    return static_cast< WeightType >( 2.0 /
      ( 1.0 + std::exp( -2.0 * radius ) ) );
    }

protected:
  TubeExponentialResolutionWeightFunction( void )
    {}
  ~TubeExponentialResolutionWeightFunction( void )
    {}

private:
  // purposely not implemented
  TubeExponentialResolutionWeightFunction( const Self & );
  // purposely not implemented
  void operator=( const Self & );
}; // End class TubeExponentialResolutionWeightFunction

} // End namespace Function

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeTubeExponentialResolutionWeightFunction_h )
