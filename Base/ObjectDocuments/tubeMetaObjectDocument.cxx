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


#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <metaUtils.h>
#include "tubeMetaObjectDocument.h"

namespace tube {


const std::string MetaObjectDocument::LABEL_NOBJECTS   = "NumberOfObjects";
const std::string MetaObjectDocument::LABEL_TYPE       = "Type";
const std::string MetaObjectDocument::LABEL_NAME       = "Name";
const std::string MetaObjectDocument::LABEL_NUM_TRANS  = "NumberOfTransforms";
const std::string MetaObjectDocument::LABEL_TRANSFORM  = "Transform";


/** Object Types *** NOTE:  MUST match the corresponding object values */
const std::string MetaObjectDocument::ID_LABEL_BLOBTYPE  = "Blob";
const std::string MetaObjectDocument::ID_LABEL_IMAGETYPE = "Image";
const std::string MetaObjectDocument::ID_LABEL_SPATIALOBJTYPE = "SpatialObject";


MetaObjectDocument::
MetaObjectDocument()
: m_NObjects( 0 ),
  m_MaxNumTransforms( 20 )
{
  if(META_DEBUG)
    {
    std::cout << "MetaObjectDocument()" << std::endl;
    }

  this->Clear();
  m_FileName[0] = '\0';
}


MetaObjectDocument::
~MetaObjectDocument()
{
}


void MetaObjectDocument::
PrintInfo() const
{
  MetaDocument::PrintInfo();
  ObjectListType::const_iterator  it = m_objects.begin();
  int index = 1;
  while( it != m_objects.end() )
    {
    std::cout << "object Number: " << index <<std::endl;
    std::cout << "objectName = " << (*it)->GetObjectName() << std::endl;
    std::cout << "NumberOfTransforms = " << (*it)->GetNumberOfTransforms() << std::endl;
    std::cout << "\n";
    ++it;
    }
}


void MetaObjectDocument::
AddObject(ObjectDocumentType::Pointer object )
{
  m_objects.push_back( object );
  m_NObjects++;
}


void MetaObjectDocument::
SetObjectList( ObjectListType& list )
{
  m_objects = list;
  m_NObjects = static_cast<int>(list.size());
}


MetaObjectDocument::ObjectListType *
MetaObjectDocument::
GetObjectList(void)
{
  return &m_objects;
}


void MetaObjectDocument::
Clear(void)
{
  if(META_DEBUG)
    {
    std::cout << "MetaObjectDocument: Clear" << std::endl;
    }

  MetaDocument::Clear();
  m_objects.clear();
}


bool MetaObjectDocument::
Read( const std::string & _fileName )
{
  if(META_DEBUG)
    {
    std::cout << "MetaObjectDocument: Read" << std::endl;
    }

  if(!_fileName.empty())
    {
    m_FileName = _fileName;
    }
  Clear();

  //Setup Document read fields and set Number of Objects read
  M_SetupReadFields();
  M_PrepareNewReadStream();

  m_ReadStream.open(m_FileName.c_str(), std::ios::binary | std::ios::in);
  m_ReadStream.seekg(0,std::ios::beg);
  if(!m_ReadStream.is_open())
    {
    std::cout << "MetaObjectDocument: Read(): Cannot open file: "
      << m_FileName << std::endl;
    return false;
    }
  bool result = M_Read();

  m_ReadStream.close();
  m_ReadStream.clear();
  return result;
}


bool MetaObjectDocument::
Write(const std::string &_fileName)
{
  if(!_fileName.empty())
    {
    m_FileName = _fileName;
    }

  M_SetupWriteFields();
  M_PrepareNewWriteStream();

#ifdef __sgi
  // Create the file. This is required on some older sgi's
  std::ofstream tFile(m_FileName.c_str(),std::ios::out);
  tFile.close();
#endif
  m_WriteStream.open(m_FileName.c_str(), std::ios::binary | std::ios::out);
  if(!m_WriteStream.is_open())
    {
    return false;
    }

  bool result = M_Write();

  m_WriteStream.close();
  m_WriteStream.clear();
  return result;
}


void MetaObjectDocument::
M_SetupReadFields(void)
{
  if(META_DEBUG)
    {
    std::cout << "MetaObjectDocument: M_SetupReadFields" << std::endl;
    }

  MetaDocument::ClearFields();
  MetaDocument::M_SetupReadFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, LABEL_NOBJECTS.c_str(), MET_INT, true);
  mF->required = true;
  mF->terminateRead = true;
  m_Fields.push_back(mF);
}


void MetaObjectDocument::
M_SetupWriteFields(void)
{
  MetaDocument::ClearFields();
  MetaDocument::M_SetupWriteFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, LABEL_NOBJECTS.c_str(), MET_INT, m_NObjects);
  mF->required = true;
  m_Fields.push_back(mF);

  for( unsigned int i = 0; i < m_objects.size(); i++ )
    {
    M_SetupObjectWriteFields(i);
    }
}


void MetaObjectDocument::
M_SetupObjectReadFields(void)
{
  MetaDocument::ClearFields();

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, LABEL_TYPE.c_str(), MET_STRING, true);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, LABEL_NAME.c_str(), MET_STRING, false);
  m_Fields.push_back(mF);

  for( unsigned int i = 0; i < m_MaxNumTransforms; i++ )
    {
    std::stringstream labelTransform;
    labelTransform << LABEL_TRANSFORM << i;
    mF = new MET_FieldRecordType();
    MET_InitReadField(mF, labelTransform.str().c_str(), MET_STRING, false);
    m_Fields.push_back(mF);
    }

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "EndObject" , MET_NONE, true);
  mF->terminateRead = true;
  mF->required = true;
  m_Fields.push_back(mF);
}


void MetaObjectDocument::
M_SetupObjectWriteFields( unsigned int object_idx )
{
  MET_FieldRecordType * mF;

  //Record the type of object
  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, LABEL_TYPE.c_str(), MET_STRING,
    m_objects[object_idx]->GetObjectType().length(),
    m_objects[object_idx]->GetObjectType().c_str() );
  m_Fields.push_back(mF);

  //Record the object Name
  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, LABEL_NAME.c_str(), MET_STRING,
    m_objects[object_idx]->GetObjectName().length(),
    m_objects[object_idx]->GetObjectName().c_str() );
  m_Fields.push_back(mF);

  //Record Names of each Transform
  for( unsigned int i=0; i < m_objects[object_idx]->GetNumberOfTransforms(); i++ )
    {
    std::stringstream label;
    label << LABEL_TRANSFORM << i;
    mF = new MET_FieldRecordType();
    MET_InitWriteField(mF, label.str().c_str(), MET_STRING,
      m_objects[object_idx]->GetTransformNames()[i].length(),
      m_objects[object_idx]->GetTransformNames()[i].c_str() );
    m_Fields.push_back(mF);
  }

  //Place end object tag
  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, "EndObject", MET_NONE );
  m_Fields.push_back(mF);
}


bool MetaObjectDocument::
M_Read(void)
{
  if(META_DEBUG)
    std::cout << "MetaObjectDocument: M_Read: Loading Header" << std::endl;

  if(!MetaDocument::M_Read())
    {
    std::cout << "MetaObjectDocument: M_Read: Error parsing file" << std::endl;
    return false;
    }

  if(META_DEBUG)
    std::cout << "MetaObjectDocument: M_Read: Parsing Header" << std::endl;

  MET_FieldRecordType * mF;

  mF = MET_GetFieldRecord(LABEL_NOBJECTS.c_str(), &m_Fields);
  if(mF != NULL && mF->defined)
    {
    m_NObjects = (int)mF->value[0];
    }


  //Iterate through the objects
  for( int i = 0; i < m_NObjects; i++ )
    {
    //Setup the Object read fields
    M_SetupObjectReadFields();

    MET_Read(m_ReadStream, & m_Fields, '=');

    ObjectDocumentType::Pointer object = ObjectDocumentType::New();

    mF = MET_GetFieldRecord( LABEL_TYPE.c_str(), &m_Fields );
    if(mF != NULL && mF->defined)
      {
      const std::string objectType = (const char *)mF->value;
      if( objectType == ID_LABEL_IMAGETYPE )
          {
          if(META_DEBUG)
            std::cout << "Reading in an image: " << objectType << std::endl;
          object = ImageDocumentType::New();
          }
        else if( objectType == ID_LABEL_BLOBTYPE )
          {
          if(META_DEBUG)
            std::cout << "Reading in a blob: "  << objectType << std::endl;
          object = BlobSpatialObjectDocumentType::New();
          }
        else if( objectType == ID_LABEL_SPATIALOBJTYPE )
          {
          if(META_DEBUG)
            std::cout << "Reading in a spatial object: "  << objectType << std::endl;
          object = SpatialObjectDocumentType::New();
          }
        else
          {
          std::cerr << "Error: Object field type does not match any existing list of types for Object #" << i << std::endl;
          return false;
          }
      }

    //Read Object Name
    mF = MET_GetFieldRecord( LABEL_NAME.c_str(), &m_Fields );
    if(mF != NULL && mF->defined)
      {
      object->SetObjectName( (const char *)mF->value );
      }

    //Read Transform
    for( unsigned int j = 0; j < m_MaxNumTransforms; j++ )
      {
      std::stringstream labelTransform;
      labelTransform << LABEL_TRANSFORM << j;
      mF = MET_GetFieldRecord( labelTransform.str().c_str(), &m_Fields );
      if(mF != NULL && mF->defined)
        {
        // Vector reference count starts 1, not 0
        object->AddTransformNameToBack( (const char *)mF->value );
        if(META_DEBUG) std::cout <<" Transform : " << (const char *)mF->value <<std::endl;
        }
      }
    m_objects.push_back( object );
    }
  return true;
}


bool MetaObjectDocument::
M_Write(void)
{
  if(!MetaDocument::M_Write())
    {
    std::cout << "MetaObjectDocument: M_Read: Error parsing file" << std::endl;
    return false;
    }

  std::ofstream fp;

  fp.open(m_FileName.c_str());
  if(!fp.is_open())
    {
    std::cout << "can't open file " << m_FileName << std::endl;
    return false;
    }

  if(!MET_Write(fp, &m_Fields))
    {
    std::cout << "MetaObject: Write: MET_Write Failed" << std::endl;
    return false;
    }

  fp.close();
  return true;
}

} // End namespace tube
