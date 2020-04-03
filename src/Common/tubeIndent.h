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

#ifndef __tubeIndent_h
#define __tubeIndent_h

#include "tubeMacro.h"

namespace tube
{

/**
 * Controls indentation during printing for classes not derived from classes in
 * ITK. For ITK-derived classes use itk::Indent.
 *
 * \ingroup  Common
 */
class Indent
{
public:

  typedef Indent        Self;
  typedef Self *        Pointer;
  typedef const Self *  ConstPointer;

  /** Constructor with initial level of indentation. */
  Indent( unsigned int indent = 0 )
    {
    m_Indent = indent;
    }

  /** Determine and return the next indentation. */
  Indent GetNextIndent( void );

  friend std::ostream & operator<<( std::ostream & os, const Indent & indent );

private:

  unsigned int m_Indent;  // Current indentation level.

}; // End class Indent

/** Print out the specified indentation. */
std::ostream & operator<<( std::ostream & os, const Indent & indent );

} // End namespace tube

#endif // End !defined( __tubeIndent_h )
