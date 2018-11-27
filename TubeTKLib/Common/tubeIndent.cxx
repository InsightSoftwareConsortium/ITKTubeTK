/*=========================================================================

Library:   TubeTKLib

Copyright Kitware Inc.

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

#include "tubeIndent.h"

namespace
{

const unsigned int TUBE_NUMBER_OF_SPACES_PER_INDENT = 2;
const unsigned int TUBE_MAXIMUM_NUMBER_OF_SPACES = 40;

const char spaces[TUBE_MAXIMUM_NUMBER_OF_SPACES + 1] =
  "                                        ";

} // End namespace

namespace tube
{

// Determine and return the next indentation.
Indent Indent::GetNextIndent( void )
{
  unsigned int indent = m_Indent + TUBE_NUMBER_OF_SPACES_PER_INDENT;

  if( indent > TUBE_MAXIMUM_NUMBER_OF_SPACES )
    {
    indent = TUBE_MAXIMUM_NUMBER_OF_SPACES;
    }

  return indent;
}

// Print out the specified indentation.
std::ostream & operator<<( std::ostream & os, const Indent & indent )
{
  os << spaces + ( TUBE_MAXIMUM_NUMBER_OF_SPACES - indent.m_Indent );

  return os;
}

} // End namespace tube
