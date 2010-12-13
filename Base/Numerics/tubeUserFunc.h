/*=========================================================================

Library:   TubeTK/VTree

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef USERFUNC_H
#define USERFUNC_H

#include <vnl/vnl_vector.h>

/*! UserFunc Derivation Examples
 *  \example TestOptimizerND/testOptimizerND.cpp
 */

namespace tube {

/*! Derive this class to pass functions to Spline and Optimization Classes
 * \author Stephen R. Aylward
 * \date 11/22/99
 */
template <class InVarT, class OutVarT>
class UserFunc
{
public :

  virtual ~UserFunc() = 0;

  /** Derive this function */
  virtual const OutVarT & value( const InVarT & x ) = 0;

};

template <class InVarT, class OutVarT>
inline UserFunc< InVarT, OutVarT >::~UserFunc()
{
}

}; // namespace tube

#endif
