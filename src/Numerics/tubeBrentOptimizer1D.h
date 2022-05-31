/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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

#ifndef __tubeBrentOptimizer1D_h
#define __tubeBrentOptimizer1D_h

#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

namespace tube
{

class BrentOptimizer1D : public Optimizer1D
{
public:

  typedef BrentOptimizer1D                    Self;
  typedef Optimizer1D                         Superclass;
  typedef Self *                              Pointer;
  typedef const Self *                        ConstPointer;

  typedef Superclass::ValueFunctionType       ValueFunctionType;
  typedef Superclass::DerivativeFunctionType  DerivativeFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( BrentOptimizer1D );

  /** Constructor. */
  BrentOptimizer1D( void );

  /** Constructor. */
  BrentOptimizer1D( ValueFunctionType::Pointer funcVal,
    DerivativeFunctionType::Pointer funcDeriv );

  /** Destructor. */
  ~BrentOptimizer1D( void );

  tubeGetMacro( Epsilon, double );

  tubeSetMacro( Epsilon, double );

  void Use( ValueFunctionType::Pointer funcVal,
    DerivativeFunctionType::Pointer funcDeriv );

protected:

  /** Print out information about this object. */
  void PrintSelf( std::ostream & os, Indent indent ) const override;

  void m_Move( double & a, double & b, double & c,
    double d, double e, double f );

  bool m_Extreme( double * x, double * xVal ) override;

private:

  // Copy constructor not implemented.
  BrentOptimizer1D( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

  double m_Epsilon;

}; // End class BrentOptimizer1D

} // End namespace tube

#endif // End !defined( __tubeBrentOptimizer1D_h )
