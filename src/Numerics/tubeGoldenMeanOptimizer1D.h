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

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __tubeGoldenMeanOptimizer1D_h
#define __tubeGoldenMeanOptimizer1D_h

#include "tubeOptimizer1D.h"
#include "tubeUserFunction.h"

namespace tube
{

class GoldenMeanOptimizer1D : public Optimizer1D
{
public:

  typedef GoldenMeanOptimizer1D          Self;
  typedef Optimizer1D                    Superclass;
  typedef Self *                         Pointer;
  typedef const Self *                   ConstPointer;

  typedef Superclass::ValueFunctionType  ValueFunctionType;

  /** Return the type of this object. */
  tubeTypeMacro( GoldenMeanOptimizer1D );

  /** Constructor. */
  GoldenMeanOptimizer1D( void );

  /** Constructor. */
  GoldenMeanOptimizer1D( ValueFunctionType::Pointer funcVal );

  /** Destructor. */
  ~GoldenMeanOptimizer1D( void );

  void Use( ValueFunctionType::Pointer funcVal );

protected:

  bool m_Extreme( double * x, double * xVal ) override;

private:

  // Copy constructor not implemented.
  GoldenMeanOptimizer1D( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class GoldenMeanOptimizer1D

} // End namespace tube

#endif // End !defined( __tubeGoldenMeanOptimizer1D_h )
