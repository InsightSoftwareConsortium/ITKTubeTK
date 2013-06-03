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

#include "tubeMetaDocument.h"

namespace tube
{

MetaDocument::
MetaDocument( void )
{
  this->Clear();
}


MetaDocument::
MetaDocument(const std::string & _fileName)
{
  this->Clear();
  this->Read(_fileName);
}


MetaDocument::
~MetaDocument( void )
{
  this->ClearFields();
}


void MetaDocument::
ClearFields( void )
{
  if(META_DEBUG)
    {
    std::cout << "MetaDocument:ClearFields" << std::endl;
    }

  for(std::vector<MET_FieldRecordType *>::iterator it = m_Fields.begin();
    it != m_Fields.end(); ++it)
    {
    delete *it;
    *it = NULL;
    }

  m_Fields.clear();
}


void MetaDocument::
FileName(const std::string & _fileName)
{
  if(!_fileName.empty())
    {
    m_FileName = _fileName;
    }
}


std::string MetaDocument::
FileName( void ) const
{
  return m_FileName;
}


void MetaDocument::
CopyInfo(const MetaDocument * _object)
{
  m_DateLastModified = _object->DateLastModified();
  m_Comment = _object->Comment();
  m_Name = _object->Name();
}


bool MetaDocument::
Read(const std::string & _fileName)
{
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: Read" << std::endl;
    }

  if(!_fileName.empty())
  {
    m_FileName = _fileName;
  }

  Clear();

  M_SetupReadFields();
  M_PrepareNewReadStream();

  m_ReadStream.open(m_FileName.c_str());
  if(!m_ReadStream.is_open())
    {
    std::cout << "MetaDocument: Read: Cannot open file" << std::endl;
    return false;
    }

  bool result = M_Read();

  m_ReadStream.close();
  m_ReadStream.clear();
  return result;
}


bool MetaDocument::
Write(const std::string & _fileName)
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
  m_WriteStream.open(m_FileName.c_str(),std::ios::binary | std::ios::out);
  if(!m_WriteStream.is_open())
    {
    return false;
    }

  bool result = M_Write();

  m_WriteStream.close();
  m_WriteStream.clear();

  return result;
}


void MetaDocument::
PrintInfo( void ) const
{
  std::cout << "Date Modified = _" << m_DateLastModified << "_" << std::endl;
  std::cout << "Comment = _" << m_Comment << "_" << std::endl;
  std::cout << "Name = " << m_Name << std::endl;
}


std::string MetaDocument::
DateLastModified( void ) const
{
  return m_DateLastModified;
}


void MetaDocument::
DateLastModified( const std::string & _date )
{
  m_DateLastModified = _date;
}


std::string MetaDocument::
Comment( void ) const
{
  return m_Comment;
}


void MetaDocument::
Comment(const std::string & _comment)
{
  m_Comment = _comment;
}


void  MetaDocument::
Name(const std::string & _Name)
{
  if(!_Name.empty())
    {
    m_Name = _Name;
    }
}


std::string MetaDocument::
Name( void ) const
{
  return m_Name;
}

void MetaDocument::
Clear( void )
{
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: Clear()" << std::endl;
    }

  m_Comment.clear();
  m_DateLastModified.clear();
  m_Name.clear();

  if(META_DEBUG)
    {
    std::cout << "MetaDocument: Clear: m_Name=" << m_Name << std::endl;
    }

  this->ClearFields();
}


void MetaDocument::
M_SetupReadFields( void )
{
  this->ClearFields();
  if(META_DEBUG)
    {
    std::cout << "MetaDocument: M_SetupReadFields" << std::endl;
    }

  MET_FieldRecordType * mF;

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "Comment", MET_STRING, false);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "DateLastModified", MET_STRING, false);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "Name", MET_STRING, false);
  m_Fields.push_back(mF);
}


void MetaDocument::
M_SetupWriteFields( void )
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

  if(!m_Comment.empty())
    {
    mF = new MET_FieldRecordType();
    MET_InitWriteField(mF, "Comment", MET_STRING, m_Comment.length(), m_Comment.c_str());
    m_Fields.push_back(mF);
    }

  if(!m_DateLastModified.empty())
    {
    mF = new MET_FieldRecordType();
    MET_InitWriteField(mF, "DateLastModified", MET_STRING, m_DateLastModified.length(),
                      m_DateLastModified.c_str());
    m_Fields.push_back(mF);
    }

  if(!m_Name.empty())
    {
    mF = new MET_FieldRecordType();
    MET_InitWriteField(mF, "Name", MET_STRING, m_Name.length(), m_Name.c_str());
    m_Fields.push_back(mF);
    }
}


bool MetaDocument::
M_Read( void )
{
  if(!MET_Read(m_ReadStream, & m_Fields, '='))
    {
    std::cout << "MetaDocument: Read: MET_Read Failed" << std::endl;
    return false;
    }

  MET_FieldRecordType * mF;

  mF = MET_GetFieldRecord("DateLastModified", &m_Fields);
  if(mF != NULL && mF->defined)
    {
    m_DateLastModified = (const char *)mF->value;
    }

  mF = MET_GetFieldRecord("Comment", &m_Fields);
  if(mF != NULL && mF->defined)
    {
    m_Comment = (const char *)mF->value;
    }


  mF = MET_GetFieldRecord("Name", &m_Fields);
  if(mF != NULL && mF->defined)
    {
    m_Name = (const char *)mF->value;
    }
  return true;
}


bool MetaDocument::
M_Write( void )
{
  if(!MET_Write(m_WriteStream, & m_Fields))
    {
    std::cout << "MetaDocument: Write: MET_Write Failed" << std::endl;
    return false;
    }

  return true;
}


void MetaDocument::
M_PrepareNewReadStream( void )
{
  if(m_ReadStream.is_open())
    {
    m_ReadStream.close();
    }
  m_ReadStream.clear();
}


void MetaDocument::
M_PrepareNewWriteStream( void )
{
  if(m_WriteStream.is_open())
    {
    m_WriteStream.close();
    }
  m_WriteStream.clear();
}

} // End namespace tube
