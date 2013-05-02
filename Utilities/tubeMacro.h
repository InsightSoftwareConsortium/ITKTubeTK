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
/**
 * tubeMacro.h defines standard system-wide macros, constants, and other
 * parameters. One of its most important functions is to define macros used
 * to interface to instance variables in a standard fashion. For example,
 * these macros manage modified time, debugging information, and provide a
 * standard interface to set and get instance variables.  Macros are
 * available for built-in types; for string classe; vector arrays;
 * object pointers; and debug, warning, and error printout information.
 */

#ifndef __tubeMacro_h
#define __tubeMacro_h

#include <string>
#include <cstdlib>

#ifndef NDEBUG
#include <cassert>
#endif

// Determine type of string stream to use.
#if !defined(CMAKE_NO_ANSI_STRING_STREAM)
#include <sstream>
#elif !defined(CMAKE_NO_ANSI_STREAM_HEADERS)
#include <strstream>
#else
#include <strstream.h>
#endif

/** A convenience macro marks variables as not being used by a method,
 * avoiding compile-time warnings. */
#define tubeNotUsed(x)

namespace tube // This is here for documentation purposes.
{
} // End namespace tube

#endif // End !defined(__tubeMacro_h)
