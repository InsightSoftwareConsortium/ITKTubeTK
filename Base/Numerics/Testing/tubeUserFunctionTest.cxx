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

#include "tubeMacro.h"

#include <vnl/vnl_vector.h>

namespace tube
{

/** Derive this class to pass functions to Spline and Optimization Classes
 * \author Stephen R. Aylward
 * \date 11/22/99
 */
class UserFunction2
{
public:

  virtual ~UserFunction2( void ) = 0;

  /** Derive this function */
  virtual const vnl_vector<double> & Value( const vnl_vector<double> & x )
    = 0;

}; // End class UserFunction2

inline UserFunction2::~UserFunction2( void )
{
}

} // End namespace tube

class MyFunc2 : public tube::UserFunction2
{
public:

  MyFunc2( void )
    {
    cVal.set_size( 1 );
    }
  const vnl_vector<double> & Value( const vnl_vector<double> & x )
    {
    std::cout << "func:x = " << x[0] << ", " << x[1] << std::endl;
    cVal[0] = std::sin( x[0] ) + std::cos( x[1]/2 );
    std::cout << "  val = " << cVal << std::endl;
    return cVal;
    }
private:
  vnl_vector<double> cVal;

}; // End class MyFunc2

int tubeUserFunctionTest( int tubeNotUsed( argc ), char * tubeNotUsed( argv )[] )
{
  MyFunc2 myFunc;

  vnl_vector<double> xTest( 2 );
  xTest[0] = 0.01;
  xTest[1] = 0.01;

  tube::UserFunction2 * op = &myFunc;
  std::cout << "test:func( 0.01,0.01 ) = " << op->Value( xTest )[0]
    << std::endl;

  std::cout << "test:func( 0.01,0.01 ) = " << myFunc.Value( xTest )[0]
    << std::endl;

  return EXIT_SUCCESS;
}
