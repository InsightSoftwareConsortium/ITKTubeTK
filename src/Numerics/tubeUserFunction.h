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

#ifndef __tubeUserFunction_h
#define __tubeUserFunction_h

#include "tubeMacro.h"

namespace tube
{

/** Derive this class to pass functions to optimizer and spline and classes. */
template< class TInput, class TOutput >
class UserFunction
{
public:

  typedef UserFunction  Self;
  typedef Self *        Pointer;
  typedef const Self *  ConstPointer;

  typedef TInput        InputType;
  typedef TOutput       OutputType;

  /** Return the type of this object. */
  tubeTypeMacroNoOverride( UserFunction );

  /** Constructor. */
  UserFunction( void )
    {
    }

  /** Destructor. */
  virtual ~UserFunction( void )
    {
    }

  /** Derive this function. */
  virtual const OutputType & Value( const InputType & input ) = 0;

private:

  // Copy constructor not implemented.
  UserFunction( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class UserFunction

} // End namespace tube

#endif // End !defined( __tubeUserFunction_h )
