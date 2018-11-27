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

#ifndef __tubeObject_h
#define __tubeObject_h

#include "tubeIndent.h"
#include "tubeMacro.h"

namespace tube
{

/**
 * Lightweight base class for classes not derived from classes in ITK or VTK.
 * For ITK-derived classes use itk::Object and for VTK-derived classes use
 * vtkObject.
 *
 * \ingroup  Common
 */
class Object
{
public:

  typedef Object        Self;
  typedef Self *        Pointer;
  typedef const Self *  ConstPointer;

  tubeTypeMacro( Object );

  /** Constructor. */
  Object( void );

  /** Destructor. */
  virtual ~Object( void );

  /** Print out information about this object. */
  void Print( std::ostream & os, Indent indent = 0 ) const;

protected:

  /** Header for when printing out information about this object. */
  virtual void PrintHeader( std::ostream & os, Indent indent ) const;

  /** Print out information about the member variables of this object. */
  virtual void PrintSelf( std::ostream & os, Indent indent ) const;

  /** Footer for when printing out information about this object. */
  virtual void PrintFooter( std::ostream & os, Indent indent ) const;

private:

  // Copy constructor not implemented.
  Object( const Self & self );

  // Copy assignment operator not implemented.
  void operator=( const Self & self );

}; // End class Object

/** Print out information about the specified object. */
std::ostream & operator<<( std::ostream & os, const Object & object );

} // End namespace tube

#endif // End !defined( __tubeObject_h )
