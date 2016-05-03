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
 *
 * Defines macros for defining set/get functions of inputs of shadow classes
 * in TubeTK/ITKModukes that do python wrapping of TubeTK/Applications
 */

#ifndef __tubeWrappingMacros_h
#define __tubeWrappingMacros_h

/** Get input of fundamental type */
#define tubeWrapGetMacro( name, type, wrapfiltername )   \
  type Get##name( void ) const                           \
    {                                                    \
    return this->m_##wrapfiltername->Get##name();        \
    }

/** Get pointer to input of object type. */
#define tubeWrapGetObjectMacro( name, type, wrapfiltername )    \
  type * Get##name( void )                                      \
    {                                                           \
    return this->m_##wrapfiltername->Get##name();               \
    }

/** Get a const pointer to input of object type. */
#define tubeWrapGetConstObjectMacro( name, type, wrapfiltername )    \
  const type * Get##name( void ) const                               \
    {                                                                \
    return this->m_##wrapfiltername->Get##name();                    \
    }

/** Set input of fundamental type */
#define tubeWrapSetMacro( name, type, wrapfiltername )   \
  void Set##name( const type value )                     \
    {                                                    \
    this->m_##wrapfiltername->Set##name( value );        \
    }

/** Set input of fundamental type */
#define tubeWrapSetObjectMacro( name, type, wrapfiltername )   \
  void Set##name( type * value )                               \
    {                                                          \
    this->m_##wrapfiltername->Set##name( value );              \
    }

/** Set input of fundamental type */
#define tubeWrapSetConstObjectMacro( name, type, wrapfiltername )   \
  void Set##name( const type * value )                              \
    {                                                               \
    this->m_##wrapfiltername->Set##name( value );                   \
    }

#endif
