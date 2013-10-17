/*=========================================================================
  MetaIO
  Copyright 2000-2010 Insight Software Consortium

  Distributed under the OSI-approved BSD License (the "License");
  see accompanying file Copyright.txt for details.

  This software is distributed WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the License for more information.
=========================================================================*/
#include "tubeMetaUtils.h"

#include "metaUtils.h"

//
// Read the type of the object
//
METAIO_STL::string MET_ReadFormTypeName(METAIO_STREAM::istream &_fp)
  {
  METAIO_STL::streampos pos = _fp.tellg();
  METAIO_STL::vector<MET_FieldRecordType *> fields;
  MET_FieldRecordType* mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "FormTypeName", MET_STRING, false);
  mF->required = false;
  mF->terminateRead = true;
  fields.push_back(mF);

  MET_Read(_fp, &fields, '=', true);
  _fp.seekg(pos);

  if(mF->defined)
    {
    METAIO_STL::string value = (char *)(mF->value);
    delete mF;
    return value;
    }

  delete mF;
  return METAIO_STL::string();
  }
