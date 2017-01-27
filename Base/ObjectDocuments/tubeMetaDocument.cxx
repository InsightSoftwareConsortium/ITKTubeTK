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

#include "tubeMetaDocument.h"

namespace tube
{

MetaDocument::MetaDocument( void )
{
}

MetaDocument::~MetaDocument( void )
{
  this->ClearFields();
}

void MetaDocument::Clear( void )
{
  this->SetComment( "" );
  this->SetDateModified( "" );
  this->SetName( "" );
  this->ClearFields();
}

void MetaDocument::CopyInformation( const Self * data )
{
  this->SetComment( data->GetComment() );
  this->SetDateModified( data->GetDateModified() );
  this->SetName( data->GetName() );
}

bool MetaDocument::Read( const std::string & fileName )
{
  if( !fileName.empty() )
    {
    this->SetFileName( fileName );
    }

  this->Clear();
  this->SetupReadFields();

  if( m_ReadStream.is_open() )
    {
    m_ReadStream.close();
    }

  m_ReadStream.clear();
  m_ReadStream.open( m_FileName.c_str(), std::ios::binary | std::ios::in );

  if( !m_ReadStream.is_open() )
    {
    return false;
    }

  const bool result = this->ReadFields();

  m_ReadStream.close();
  m_ReadStream.clear();

  return result;
}

bool MetaDocument::Write( const std::string & fileName )
{
  if( !fileName.empty() )
    {
    this->SetFileName( fileName );
    }

  this->SetupWriteFields();

  if( m_WriteStream.is_open() )
    {
    m_WriteStream.close();
    }

  m_WriteStream.clear();
  m_WriteStream.open( m_FileName.c_str(), std::ios::binary | std::ios::out );

  if( !m_WriteStream.is_open() )
    {
    return false;
    }

  const bool result = this->WriteFields();

  m_WriteStream.close();
  m_WriteStream.clear();

  return result;
}

void MetaDocument::ClearFields( void )
{
  for( FieldListType::iterator it = m_FieldList.begin();
       it != m_FieldList.end(); ++it )
    {
    delete *it;
    *it = NULL;
    }

  m_FieldList.clear();
}

bool MetaDocument::ReadFields( void )
{
  if( !MET_Read( m_ReadStream, & m_FieldList, '=' ) )
    {
    return false;
    }

  FieldType * const dateModifiedField = MET_GetFieldRecord(
    "DateLastModified", &m_FieldList );

  if( dateModifiedField != NULL && dateModifiedField->defined )
    {
    this->SetDateModified( ( const char * )dateModifiedField->value );
    }

  FieldType * const commentField = MET_GetFieldRecord( "Comment",
    &m_FieldList );

  if( commentField != NULL && commentField->defined )
    {
    this->SetComment( ( const char * )commentField->value );
    }

  FieldType * const nameField = MET_GetFieldRecord( "Name", &m_FieldList );

  if( nameField != NULL && nameField->defined )
    {
    this->SetName( ( const char * )nameField->value );
    }

  return true;
}

void MetaDocument::SetupReadFields( void )
{
  this->ClearFields();

  FieldType * const commentField = new FieldType();
  MET_InitReadField( commentField, "Comment", MET_STRING, false );
  m_FieldList.push_back( commentField );

  FieldType * const dateModifiedField = new FieldType();
  MET_InitReadField( dateModifiedField, "DateLastModified", MET_STRING,
    false );
  m_FieldList.push_back( dateModifiedField );

  FieldType * const nameField = new FieldType();
  MET_InitReadField( nameField, "Name", MET_STRING, false );
  m_FieldList.push_back( nameField );
}

void MetaDocument::SetupWriteFields( void )
{
  this->ClearFields();

  if( std::strlen( this->GetComment() ) > 0 )
    {
    FieldType * const commentField = new FieldType();
    MET_InitWriteField( commentField, "Comment", MET_STRING,
      std::strlen( this->GetComment() ), this->GetComment() );
    m_FieldList.push_back( commentField );
    }

  if( std::strlen( this->GetDateModified() ) > 0 )
    {
    FieldType * const dateModifiedField = new FieldType();
    MET_InitWriteField( dateModifiedField, "DateLastModified", MET_STRING,
      std::strlen( this->GetDateModified() ), this->GetDateModified() );
    m_FieldList.push_back( dateModifiedField );
    }

  if( std::strlen( this->GetName() ) > 0 )
    {
    FieldType * const nameField = new FieldType();
    MET_InitWriteField( nameField, "Name", MET_STRING,
      std::strlen( this->GetName() ), this->GetName() );
    m_FieldList.push_back( nameField );
    }
}

bool MetaDocument::WriteFields( void )
{
  if( !MET_Write( m_WriteStream, & m_FieldList ) )
    {
    return false;
    }

  return true;
}

void MetaDocument::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Comment:      " << m_Comment << std::endl;
  os << indent << "DateModified: " << m_Comment << std::endl;
  os << indent << "Name:         " << m_Name << std::endl;
  os << indent << "FileName:     " << m_FileName << std::endl;
  os << indent << "FieldList:    " << std::endl;

  for( FieldListType::const_iterator it = m_FieldList.begin();
       it != m_FieldList.end(); ++it )
    {
    os << indent << *it << std::endl;
    }
}

} // End namespace tube
