#include "metaTypes.h"

#ifndef ITKMetaIO_METAUTILS_TEMP_H
#define ITKMetaIO_METAUTILS_TEMP_H

template <class T>
bool MET_InitWriteField_Temp(MET_FieldRecordType * _mf,
                                   const char *_name,
                                   MET_ValueEnumType _type,
                                   size_t _length,
                                   T *_v)
  {
  strncpy(_mf->name, _name,254);
  _mf->name[254] = '\0';
  _mf->type = _type;
  _mf->defined = true;
  _mf->length = static_cast<int>(_length);
  _mf->dependsOn = -1;
  _mf->required = false;
  _mf->terminateRead = false;
  if(_type == MET_FLOAT_MATRIX)
    {
    size_t i;
    for(i=0; i < 4096 && i < _length*_length; i++)
      {
      _mf->value[i] = (double)(_v[i]);
      }
    }
  else if(_type != MET_STRING)
    {
    size_t i;
    for(i=0; i < 4096 && i < _length; i++)
      {
      _mf->value[i] = (double)(_v[i]);
      }
    }
  else
    {
    strncpy((char *)(_mf->value), (const char *)_v,
            (sizeof(_mf->value)-1));
    ((char *)(_mf->value))[(sizeof(_mf->value)-1)] = '\0';
    }
  return true;
  }

#endif
