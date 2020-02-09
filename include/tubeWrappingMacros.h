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

/**
 *
 * Defines macros for defining set/get functions of inputs of shadow classes
 * in TubeTK/ITKModukes that do python wrapping of TubeTK/Applications
 *
 *  \ingroup TubeTK
 */

#ifndef __tubeWrappingMacros_h
#define __tubeWrappingMacros_h

/** Boolean macro */
#define tubeWrapBooleanMacro( name, wrap_filter_object_name )   \
  void name##On( void ) const                            \
    {                                                    \
    this->m_##wrap_filter_object_name->name##On();       \
    }                                                    \
  void name##Off( void ) const                           \
    {                                                    \
    this->m_##wrap_filter_object_name->name##Off();      \
    }

/** Get input of fundamental type */
#define tubeWrapGetMacro( name, type, wrap_filter_object_name )   \
  type Get##name( void ) const                            \
    {                                                             \
    return this->m_##wrap_filter_object_name->Get##name();        \
    }

/** Get pointer to input of object type. */
#define tubeWrapGetObjectMacro( name, type, wrap_filter_object_name )    \
  type * Get##name( void )                                       \
    {                                                                    \
    return this->m_##wrap_filter_object_name->Get##name();               \
    }

/** Get a const pointer to input of object type. */
#define tubeWrapGetConstObjectMacro( name, type, wrap_filter_object_name )   \
  const type * Get##name( void ) const                               \
    {                                                                        \
    return this->m_##wrap_filter_object_name->Get##name();                   \
    }

/** Get a const reference to input of object type. */
#define tubeWrapGetConstReferenceMacro( name, type, wrap_filter_object_name ) \
  const type & Get##name( void ) const                                        \
    {                                                                         \
    return this->m_##wrap_filter_object_name->Get##name();                    \
    }

/** Set input of fundamental type */
#define tubeWrapSetMacro( name, type, wrap_filter_object_name )   \
  void Set##name( type value )                              \
    {                                                             \
    if( this->m_##wrap_filter_object_name->Get##name() != value ) \
      {                                                           \
      this->m_##wrap_filter_object_name->Set##name( value );      \
      this->Modified();                                           \
      }                                                           \
    }

/** Set input of fundamental type */
#define tubeWrapForceSetMacro( name, type, wrap_filter_object_name )   \
  void Set##name( type value )                              \
    {                                                             \
    this->m_##wrap_filter_object_name->Set##name( value );      \
    this->Modified();                                           \
    }

/** Set input using pointer to object type */
#define tubeWrapSetObjectMacro( name, type, wrap_filter_object_name )   \
  void Set##name( type * value )                                \
    {                                                                   \
    if( this->m_##wrap_filter_object_name->Get##name() != value )       \
      {                                                                 \
      this->m_##wrap_filter_object_name->Set##name( value );            \
      this->Modified();                                                 \
      }                                                                 \
    }

/** Set input using pointer to object type, without testing if changed */
#define tubeWrapForceSetObjectMacro( name, type, wrap_filter_object_name ) \
  void Set##name( type * value )                                           \
    {                                                                      \
    this->m_##wrap_filter_object_name->Set##name( value );                 \
    this->Modified();                                                      \
    }

/** Set input using reference to object type */
#define tubeWrapSetReferenceMacro( name, type, wrap_filter_object_name )  \
  void Set##name( type & value )                                          \
    {                                                                     \
    if( this->m_##wrap_filter_object_name->Get##name() != value )         \
      {                                                                   \
      this->m_##wrap_filter_object_name->Set##name( value );              \
      this->Modified();                                                   \
      }                                                                   \
    }

/** Set input using reference to object type */
#define tubeWrapForceSetReferenceMacro( name, type, wrap_filter_object_name )  \
  void Set##name( const type & value )                                        \
    {                                                                   \
    this->m_##wrap_filter_object_name->Set##name( value );              \
    this->Modified();                                                   \
    }

/** Add input using reference to object type */
#define tubeWrapAddMacro( name, type, wrap_filter_object_name )  \
  void Add##name( type value )                                        \
    {                                                                   \
    this->m_##wrap_filter_object_name->Add##name( value );              \
    this->Modified();                                                   \
    }

/** Set Nth in an object list */
#define tubeWrapSetNthMacro( name, type, wrap_filter_object_name )  \
  void Set##name( unsigned int i, type value )                        \
    {                                                                   \
    this->m_##wrap_filter_object_name->Set##name( i, value );           \
    this->Modified();                                                   \
    }

/** Set input using const pointer to object type */
#define tubeWrapSetConstObjectMacro( name, type, wrap_filter_object_name ) \
  void Set##name( const type * value )                                     \
    {                                                                      \
    if( this->m_##wrap_filter_object_name->Get##name() != value )          \
      {                                                                    \
      this->m_##wrap_filter_object_name->Set##name( value );               \
      this->Modified();                                                    \
      }                                                                    \
    }

#define tubeWrapAddConstObjectMacro( name, type, wrap_filter_object_name ) \
  void Add##name( const type * value )                                     \
    {                                                                      \
    this->m_##wrap_filter_object_name->Add##name( value );               \
    this->Modified();                                                    \
    }

/** Set Nth in an const pointer object list */
#define tubeWrapSetNthObjectMacro( name, type, wrap_filter_object_name )  \
  void Set##name( unsigned int i, type * value )                        \
    {                                                                   \
    this->m_##wrap_filter_object_name->Set##name( i, value );           \
    this->Modified();                                                   \
    }

/** Set Nth in an const pointer object list */
#define tubeWrapSetNthConstObjectMacro( name, type, wrap_filter_object_name )  \
  void Set##name( unsigned int i, const type * value )                        \
    {                                                                   \
    this->m_##wrap_filter_object_name->Set##name( i, value );           \
    this->Modified();                                                   \
    }

/** Get Nth in an object list */
#define tubeWrapGetNthMacro( name, type, wrap_filter_object_name )  \
  type Get##name( unsigned int i )                        \
    {                                                                   \
    return this->m_##wrap_filter_object_name->Get##name( i );           \
    }

/** Get Nth in an object list */
#define tubeWrapGetNthObjectMacro( name, type, wrap_filter_object_name )  \
  type * Get##name( unsigned int i )                        \
    {                                                                   \
    return this->m_##wrap_filter_object_name->Get##name( i );           \
    }

/** Get Nth in an object list */
#define tubeWrapGetNthConstObjectMacro( name, type, wrap_filter_object_name )  \
  const type * Get##name( unsigned int i )                        \
    {                                                                   \
    return this->m_##wrap_filter_object_name->Get##name( i );           \
    }

/** Get Nth in an object list */
#define tubeWrapGetNthConstReferenceMacro( name, type, wrap_filter_object_name )  \
  const type & Get##name( unsigned int i )                        \
    {                                                                   \
    return this->m_##wrap_filter_object_name->Get##name( i );           \
    }

/** Set input using const pointer to object type */
#define tubeWrapForceSetConstObjectMacro( name, type,                    \
  wrap_filter_object_name )                                              \
  void Set##name( const type * value )                                   \
    {                                                                    \
    this->m_##wrap_filter_object_name->Set##name( value );               \
    this->Modified();                                                    \
    }

/** Set input using const reference to object type */
#define tubeWrapSetConstReferenceMacro( name, type, wrap_filter_object_name ) \
  void Set##name( const type & value )                                \
    {                                                                 \
    if( this->m_##wrap_filter_object_name->Get##name() != value )     \
      {                                                               \
      this->m_##wrap_filter_object_name->Set##name( value );          \
      this->Modified();                                               \
      }                                                               \
    }

/** Set input using const reference to object type */
#define tubeWrapForceSetConstReferenceMacro( name, type,              \
  wrap_filter_object_name )                                           \
  void Set##name( const type & value )                                \
    {                                                                 \
    this->m_##wrap_filter_object_name->Set##name( value );            \
    this->Modified();                                                 \
    }

/** Redirect call to a function of the same named in the wrapped filter */
#define tubeWrapCallMacro( name, wrap_filter_object_name )   \
  void name()                                                \
    {                                                        \
    this->m_##wrap_filter_object_name->name();               \
    }

/** Redirect call to a function of the same named in the wrapped filter */
#define tubeWrapCallWithConstReferenceArgMacro( name, type, wrap_filter_object_name )   \
  void name( const type & value )                                        \
    {                                                              \
    this->m_##wrap_filter_object_name->name( value );              \
    }

/** Redirect call to a function of the same named in the wrapped filter */
#define tubeWrapCallOverrideMacro( name, wrap_filter_object_name )   \
  void name() override                                       \
    {                                                        \
    this->m_##wrap_filter_object_name->name();               \
    }

/** Redirect call to Update() wrapped filter's Update() */
#define tubeWrapUpdateMacro( wrap_filter_object_name )                   \
  tubeWrapCallOverrideMacro( Update, wrap_filter_object_name )                   \

#endif
