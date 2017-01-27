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
#ifndef __metaUtilsTemp_h
#define __metaUtilsTemp_h

#include "metaTypes.h"

#ifndef ITKMetaIO_METAUTILS_TEMP_h
#define ITKMetaIO_METAUTILS_TEMP_h

namespace tube
{

template <class T>
bool MET_InitWriteField_Temp( MET_FieldRecordType * _mf,
                                   const char *_name,
                                   MET_ValueEnumType _type,
                                   size_t _length,
                                   T *_v )
{
  strncpy( _mf->name, _name, 254 );
  _mf->name[254] = '\0';
  _mf->type = _type;
  _mf->defined = true;
  _mf->length = static_cast<int>( _length );
  _mf->dependsOn = -1;
  _mf->required = false;
  _mf->terminateRead = false;
  if( _type == MET_FLOAT_MATRIX )
    {
    size_t i;
    for( i=0; i < 4096 && i < _length*_length; i++ )
      {
      _mf->value[i] = ( double )( _v[i] );
      }
    }
  else if( _type != MET_STRING )
    {
    size_t i;
    for( i=0; i < 4096 && i < _length; i++ )
      {
      _mf->value[i] = ( double )( _v[i] );
      }
    }
  else
    {
    strncpy( ( char * )( _mf->value ), ( const char * )_v,
            ( sizeof( _mf->value )-1 ) );
    ( ( char * )( _mf->value ) )[( sizeof( _mf->value )-1 )] = '\0';
    }
  return true;
}


} // End namespace tube

#endif

#endif // __metaUtilsTemp_h
