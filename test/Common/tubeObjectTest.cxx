/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

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

#include "tubeObject.h"

#include <sstream>

int tubeObjectTest( int argc, char * argv[] )
{
  if( argc > 1 )
    {
    tubeStandardErrorMacro( << "Usage: " << argv[0] );

    return EXIT_FAILURE;
    }

  typedef tube::Object ObjectType;

  ObjectType::Pointer object = new ObjectType();
  std::ostringstream oss;

  object->Print( oss );

  tubeStandardOutputMacro( << "Object print method test:" << std::endl
                           << oss.str() );
  tubeStandardOutputMacro( << "Object stream operator test:" << std::endl
                           << *object );

  delete object;

  return EXIT_SUCCESS;
}
