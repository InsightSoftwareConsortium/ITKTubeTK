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

#ifndef __tubeUserFunc_h
#define __tubeUserFunc_h

/** UserFunc Derivation Examples
 *  \example TestOptimizerND/testOptimizerND.cpp
 */

namespace tube
{

/** Derive this class to pass functions to Spline and Optimization Classes
 * \class UserFunc
 * \author Stephen R. Aylward
 * \date 11/22/99
 */
template <class InVarT, class OutVarT>
class UserFunc
{
public:
  virtual ~UserFunc( void ) = 0;

  /** Derive this function */
  virtual const OutVarT & value( const InVarT & x ) = 0;

}; // End class UserFunc

template <class InVarT, class OutVarT>
inline UserFunc< InVarT, OutVarT >::~UserFunc( void )
{
}

} // End namespace tube

#endif // End !defined(__tubeUserFunc_h)
