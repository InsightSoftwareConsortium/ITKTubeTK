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


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>

#include "MetaDocument.h"

namespace tube
{

MetaDocument::
MetaDocument(void)
{
  this->ClearFields();
  MetaDocument::Clear();
  m_ReadStream = NULL;
  m_WriteStream = NULL;
  m_FileName[0] = '\0';
}


MetaDocument::
MetaDocument(const char * _fileName)
{
  this->ClearFields();
  MetaDocument::Clear();
  m_ReadStream = NULL;
  m_WriteStream = NULL;
  this->Read(_fileName);
}


MetaDocument::
~MetaDocument(void)
{
  delete m_ReadStream;
  delete m_WriteStream;

  this->ClearFields();
}


void MetaDocument::
ClearFields()
{
  if(META_DEBUG)
  {
  std::cout << "MetaDocument:ClearFields" << std::endl;
  }
  m_Fields.clear();
}


void MetaDocument::
FileName(const char *_fileName)
{
  if(_fileName != NULL)
    {
    if(_fileName[0] != '\0')
      {
      strcpy(m_FileName, _fileName);
      }
    }
}


const char * MetaDocument::
FileName(void) const
{
  return m_FileName;
}


void MetaDocument::
CopyInfo(const MetaDocument * _object)
{
  DateLastModified(_object->DateLastModified());
  Comment(_object->Comment());
  Name(_object->Name());

}


bool MetaDocument::
Read(const char *_fileName)
{
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: Read" << std::endl;
    }

  if(_fileName != NULL)
    {
    strcpy(m_FileName, _fileName);
    }

  Clear();

  M_SetupReadFields();
  M_PrepareNewReadStream();

  m_ReadStream->open(m_FileName);
  if(!m_ReadStream->is_open())
    {
    std::cout << "MetaDocument: Read: Cannot open file" << std::endl;
    return false;
    }

  bool result = M_Read();

  m_ReadStream->close();
  m_ReadStream->clear();
  return result;
}


bool MetaDocument::
Write(const char *_fileName)
{
  if(_fileName != NULL)
    {
    FileName(_fileName);
    }

  M_SetupWriteFields();

  if(!m_WriteStream)
    {
    m_WriteStream = new std::ofstream;
    }

#ifdef __sgi
  // Create the file. This is required on some older sgi's
  std::ofstream tFile(m_FileName,std::ios::out);
  tFile.close();
#endif
  m_WriteStream->open(m_FileName,std::ios::binary | std::ios::out);
  if(!m_WriteStream->is_open())
    {
    return false;
    }

  bool result = M_Write();

  m_WriteStream->close();
  delete m_WriteStream;
  m_WriteStream = 0;

  return result;
}


void MetaDocument::
PrintInfo(void) const
{
  std::cout << "Date Modified = _" << m_DateLastModified << "_" << std::endl;
  std::cout << "Comment = _" << m_Comment << "_" << std::endl;
  std::cout << "Name = " << m_Name << std::endl;
}


const char* MetaDocument::
DateLastModified(void) const
{
  return m_DateLastModified;
}


void MetaDocument::
DateLastModified( const char* _date )
{
  strcpy(m_DateLastModified, _date );
}


const char * MetaDocument::
Comment(void) const
{
  return m_Comment;
}


void MetaDocument::
Comment(const char * _comment)
{
  strcpy(m_Comment, _comment);
}


void  MetaDocument::
Name(const char *_Name)
{
  if(_Name != NULL)
    {
    strcpy(m_Name, _Name);
    }
}


const char  * MetaDocument::
Name(void) const
{
  return m_Name;
}

void MetaDocument::
Clear(void)
{
  if(META_DEBUG)  std::cout << "MetaDocument: Clear()" << std::endl;
  strcpy(m_Comment, "");
  strcpy(m_DateLastModified, "");
  strcpy(m_Name, "");

  if(META_DEBUG)
    {
    std::cout << "MetaDocument: Clear: m_Name=" << m_Name << std::endl;
    }
  this->ClearFields();
}


void MetaDocument::
M_SetupReadFields(void)
{
  this->ClearFields();
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: M_SetupReadFields" << std::endl;
    }

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Comment", MET_STRING, false);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "DateLastModified", MET_STRING, false);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType;
  MET_InitReadField(mF, "Name", MET_STRING, false);
  m_Fields.push_back(mF);
}


void MetaDocument::
M_SetupWriteFields(void)
{
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: M_SetupWriteFields" << std::endl;
    }

  this->ClearFields();

  if(META_DEBUG)
    {
    std::cout << "MetaDocument: M_SetupWriteFields: Creating Fields"
              << std::endl;
    }

  MET_FieldRecordType * mF;

  if(strlen(m_Comment)>0)
    {
    mF = new MET_FieldRecordType;
    MET_InitWriteField(mF, "Comment", MET_STRING, strlen(m_Comment), m_Comment);
    m_Fields.push_back(mF);
    }

  if(strlen(m_DateLastModified)>0)
    {
    mF = new MET_FieldRecordType;
    MET_InitWriteField(mF, "DateLastModified", MET_STRING, strlen(m_DateLastModified),
                      m_DateLastModified);
    m_Fields.push_back(mF);
    }

  if(strlen(m_Name)>0)
    {
    mF = new MET_FieldRecordType;
    MET_InitWriteField(mF, "Name", MET_STRING, strlen(m_Name),m_Name);
    m_Fields.push_back(mF);
    }
}


bool MetaDocument::
M_Read(void)
{
  if(!MET_Read(*m_ReadStream, & m_Fields, '='))
    {
    std::cout << "MetaDocument: Read: MET_Read Failed" << std::endl;
    return false;
    }

  MET_FieldRecordType * mF;

  mF = MET_GetFieldRecord("DateLastModified", &m_Fields);
  if(mF && mF->defined)
    {
    strcpy(m_DateLastModified, (char *)(mF->value));
    }

  mF = MET_GetFieldRecord("Comment", &m_Fields);
  if(mF && mF->defined)
    {
    strcpy(m_Comment, (char *)(mF->value));
    }


  mF = MET_GetFieldRecord("Name", &m_Fields);
  if(mF && mF->defined)
    {
    strcpy(m_Name, (char *)(mF->value));
    }
  return true;
}


bool MetaDocument::
M_Write(void)
{
  if(!MET_Write(*m_WriteStream, & m_Fields))
    {
    std::cout << "MetaDocument: Write: MET_Write Failed" << std::endl;
    return false;
    }

  return true;
}


void MetaDocument::
M_PrepareNewReadStream()
{
  if(m_ReadStream)
    {
    if(m_ReadStream->is_open())
      {
      m_ReadStream->close();
      }
    m_ReadStream->clear();
    }
  else
    {
    m_ReadStream = new std::ifstream;
    }
}

} // End namespace tube
