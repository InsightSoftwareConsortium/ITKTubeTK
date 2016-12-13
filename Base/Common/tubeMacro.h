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

/** Defines standard macros for classes not derived from classes in ITK or VTK.
    For ITK-derived classes use <itkMacro.h> and for VTK-derived classes use
    <vtkSetGet.h>. */

#ifndef __tubeMacro_h
#define __tubeMacro_h

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

// Compiler support, based on ITK and VTK requirements.

// Borland C++ compiler is not supported.
#ifdef __BORLANDC__
#error "__BORLANDC__ is not supported by TubeTK"
#endif

// Cygwin compiler is not supported.
#ifdef __CYGWIN__
#error "__CYGWIN__ is not supported by TubeTK"
#endif

// GNU compiler, version less than 3, is not supported.
#if defined( __GNUC__ ) && ( __GNUC__ < 3 )
#error "__GNUC__ < 3 is not supported by TubeTK"
#endif

// Microsoft Visual C++ compiler, version less than 7.1, is not supported.
#if defined( _MSC_VER ) && ( _MSC_VER < 1310 )
#error "_MSC_VER < 1310 is not supported by TubeTK"
#endif

// Metrowerks compiler is not supported
#ifdef __MWERKS__
#error "__MWERKS__ is not supported by TubeTK"
#endif

// SGI compiler is not supported.
#ifdef __sgi
#error "__sgi is not supported by TubeTK"
#endif

// Sun Pro C++ compiler, version less than 0x590, is not supported.
#if defined( __SUNPRO_CC ) && ( __SUNPRO_CC < 0x590 )
#error "__SUNPRO_CC < 0x590 is not supported by TubeTK"
#endif

// Disable some warnings for the Microsoft Visual C++ compiler, based on ITK and VTK.

#ifdef _MSC_VER
// Conditional expression is constant.
#pragma warning ( disable : 4127 )

// Conversion from type T1 to type T2, possible loss of data.
#pragma warning ( disable : 4244 )

// Class T1 needs to have DLL-interface to be used by clients of class T2.
#pragma warning ( disable : 4251 )

// Inconsistent DLL linkage.
#pragma warning ( disable : 4273 )

// Non-DLL-interface class key used as base for DLL-interface class key.
#pragma warning ( disable : 4275 )

// Exception specification ignored except to indicate a function
// is not __declspec( nothrow ).
#pragma warning ( disable : 4290 )

// Truncation from type T1 to type T2.
#pragma warning ( disable : 4305 )

// Truncation of constant value.
#pragma warning ( disable : 4309 )

// Decorated name length exceeded, name was truncated.
#pragma warning ( disable : 4503 )

// Unreferenced local function has been removed.
#pragma warning ( disable : 4505 )

// Forcing value to type bool, performance warning.
#pragma warning ( disable : 4800 )
#endif

// Convenience macros, based on ITK and VTK convenience macros.

/** Mark that a function argument is unused. */
#define tubeNotUsed( x )

/** Return the name of a class. */
#define tubeTypeMacro( classname )                 \
virtual const char * GetNameOfClass( void ) const  \
  {                                                \
  return #classname;                               \
  }

/** Return a member variable of fundamental type. */
#define tubeGetMacro( variablename, variabletype )  \
virtual variabletype Get##variablename( void )      \
  {                                                 \
  return this->m_##variablename;                    \
  }

/** Return a member variable of fundamental type. */
#define tubeGetConstMacro( variablename, variabletype )  \
virtual variabletype Get##variablename( void ) const     \
  {                                                      \
  return this->m_##variablename;                         \
  }

/** Return a const reference to a member variable of object type. */
#define tubeGetConstReferenceMacro( variablename, variabletype )  \
virtual const variabletype & Get##variablename( void ) const      \
  {                                                               \
  return this->m_##variablename;                                  \
  }

/** Return a member variable of C++ string type as a C string. */
#define tubeGetStringMacro( variablename )            \
virtual const char * Get##variablename( void ) const  \
  {                                                   \
  return this->m_##variablename.c_str();              \
  }

/** Return a pointer to a member variable of object type. */
#define tubeGetModifiableObjectMacro( variablename, variabletype )  \
virtual variabletype * GetModifiable##variablename( void )          \
  {                                                                 \
  return &( this->m_##variablename );                               \
  }

/** Return a const pointer to a member variable of object type. */
#define tubeGetConstObjectMacro( variablename, variabletype )  \
virtual const variabletype * Get##variablename( void ) const   \
  {                                                            \
  return &( this->m_##variablename );                          \
  }

/** Set a member variable of fundamental type. */
#define tubeSetMacro( variablename, variabletype )                   \
virtual void Set##variablename( const variabletype _variablevalue )  \
  {                                                                  \
  this->m_##variablename = _variablevalue;                           \
  }

/** Set a member variable of object type. */
#define tubeSetConstReferenceMacro( variablename, variabletype )       \
virtual void Set##variablename( const variabletype & _variablevalue )  \
  {                                                                    \
  this->m_##variablename = _variablevalue;                             \
  }

/** Set a member variable of C++ string type. */
#define tubeSetStringMacro( variablename )                     \
virtual void Set##variablename( const char * _variablevalue )  \
  {                                                            \
  if( _variablevalue )                                         \
    {                                                          \
    this->m_##variablename = _variablevalue;                   \
    }                                                          \
  else                                                         \
    {                                                          \
    this->m_##variablename = "";                               \
    }                                                          \
  }                                                            \
virtual void Set##variablename( const std::string & _variablevalue )  \
  {                                                                   \
  this->Set##variablename( _variablevalue.c_str() );                  \
  }

/** Set a pointer to a member variable of object type. */
#define tubeSetObjectMacro( variablename, variabletype )         \
virtual void Set##variablename( variabletype * _variablevalue )  \
  {                                                              \
  this->m_##variablename = _variablevalue;                       \
  }

/** Set a const pointer to a member variable of object type. */
#define tubeSetConstObjectMacro( variablename, variabletype )          \
virtual void Set##variablename( const variabletype * _variablevalue )  \
  {                                                                    \
  this->m_##variablename = _variablevalue;                             \
  }

/** Create on and off member functions for a member variable of type bool. */
#define tubeBooleanMacro( variablename )  \
virtual void variablename##On( void )     \
  {                                       \
  this->Set##variablename( true );        \
  }                                       \
virtual void variablename##Off( void )    \
  {                                       \
  this->Set##variablename( false );       \
  }

#ifdef NDEBUG
/** Print debug information to the standard error stream. */
#define tubeDebugMacro( statement )
#else
/** Print debug information to the standard error stream. */
#define tubeDebugMacro( statement )                                 \
  {                                                                 \
  std::cerr << "Debug: In " __FILE__ ", line " << __LINE__ << ": "  \
               statement << std::endl;                              \
  }
#endif

/** Print warning information to the standard error stream. */
#define tubeWarningMacro( statement )                                 \
  {                                                                   \
  std::cerr << "Warning: In " __FILE__ ", line " << __LINE__ << ": "  \
               statement << std::endl;                                \
  }

/** Print error information to the standard error stream. */
#define tubeErrorMacro( statement )                                 \
  {                                                                 \
  std::cerr << "Error: In " __FILE__ ", line " << __LINE__ << ": "  \
               statement << std::endl;                              \
  }

/** Print to the standard error stream. */
#define tubeStandardErrorMacro( statement )  \
  {                                          \
  std::cerr statement << std::endl;          \
  }

/** Print to the standard output stream. */
#define tubeStandardOutputMacro( statement )  \
  {                                           \
  std::cout statement << std::endl;           \
  }

/** Contains all TubeTK classes. */
namespace tube
{
}

#endif // End !defined( __tubeMacro_h )
