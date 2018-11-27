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

#include "tubeMetaObjectDocument.h"

namespace tube
{

MetaObjectDocument
::MetaObjectDocument( void )
{
  m_MaximumNumberOfTransforms = 20;
  m_NumberOfObjectDocuments = 0;
}


MetaObjectDocument
::~MetaObjectDocument( void )
{
  this->Clear();
}


MetaObjectDocument::ObjectDocumentListType &
MetaObjectDocument
::GetObjectDocumentList( void )
{
  return m_ObjectDocumentList;
}


void
MetaObjectDocument
::SetObjectDocumentList( ObjectDocumentListType & objectDocumentList )
{
  m_ObjectDocumentList = objectDocumentList;
  m_NumberOfObjectDocuments = static_cast< int >( m_ObjectDocumentList.size() );
}


void
MetaObjectDocument
::AddObjectDocument( ObjectDocumentType::Pointer object )
{
  m_ObjectDocumentList.push_back( object );
  ++m_NumberOfObjectDocuments;
}


void
MetaObjectDocument
::Clear( void )
{
  this->Superclass::Clear();

  m_ObjectDocumentList.clear();
}


bool
MetaObjectDocument
::ReadFields( void )
{
  if( !this->Superclass::ReadFields() )
    {
    return false;
    }

  FieldType * mF = MET_GetFieldRecord( "NumberOfObjects", &m_FieldList );

  if( mF != NULL && mF->defined )
    {
    m_NumberOfObjectDocuments = ( int )mF->value[0];
    }

  for( int i = 0; i < m_NumberOfObjectDocuments; ++i )
    {
    this->SetupObjectReadFields();

    MET_Read( m_ReadStream, & m_FieldList, '=' );

    ObjectDocumentType::Pointer object = ObjectDocumentType::New();

    mF = MET_GetFieldRecord( "Type", &m_FieldList );

    if( mF != NULL && mF->defined )
      {
      const std::string objectType = ( const char * )mF->value;

      if( objectType == "Image" )
        {
        object = ImageDocumentType::New();
        }
      else if( objectType == "Blob" )
        {
        object = BlobSpatialObjectDocumentType::New();
        }
      else if( objectType == "SpatialObject" )
        {
        object = SpatialObjectDocumentType::New();
        }
      else
        {
        return false;
        }
      }

    mF = MET_GetFieldRecord( "Name", &m_FieldList );

    if( mF != NULL && mF->defined )
      {
      object->SetObjectName( ( const char * )mF->value );
      }

    // Read transform.
    for( unsigned int j = 0; j < m_MaximumNumberOfTransforms; ++j )
      {
      std::stringstream labelTransform;
      labelTransform << "Transform" << j;

      mF = MET_GetFieldRecord( labelTransform.str().c_str(), &m_FieldList );

      if( mF != NULL && mF->defined )
        {
        // Vector reference count starts 1, not 0.
        object->AddTransformNameToBack( ( const char * )mF->value );
        }
      }
    m_ObjectDocumentList.push_back( object );
    }

  return true;
}


void
MetaObjectDocument
::SetupObjectReadFields( void )
{
  this->Superclass::ClearFields();

  FieldType * const typeField = new FieldType();
  MET_InitReadField( typeField, "Type", MET_STRING, true );
  m_FieldList.push_back( typeField );

  FieldType * const nameField = new FieldType();
  MET_InitReadField( nameField, "Name", MET_STRING, false );
  m_FieldList.push_back( nameField );

  for( unsigned int i = 0; i < m_MaximumNumberOfTransforms; ++i )
    {
    std::stringstream labelTransform;
    labelTransform << "Transform" << i;

    FieldType * const transformNameField = new FieldType();
    MET_InitReadField( transformNameField, labelTransform.str().c_str(),
      MET_STRING, false );
    m_FieldList.push_back( transformNameField );
    }

  FieldType * const endObjectField = new FieldType();
  MET_InitReadField( endObjectField, "EndObject", MET_NONE, true );
  endObjectField->terminateRead = true;
  endObjectField->required = true;
  m_FieldList.push_back( endObjectField );
}


void
MetaObjectDocument
::SetupObjectWriteFields( unsigned int index )
{
  FieldType * const typeField = new FieldType();
  MET_InitWriteField( typeField, "Type", MET_STRING,
    std::strlen( m_ObjectDocumentList[index]->GetObjectType() ),
    m_ObjectDocumentList[index]->GetObjectType() );
  m_FieldList.push_back( typeField );

  FieldType * const nameField = new FieldType();
  MET_InitWriteField( nameField, "Name", MET_STRING,
    std::strlen( m_ObjectDocumentList[index]->GetObjectName().c_str() ),
    m_ObjectDocumentList[index]->GetObjectName().c_str() );
  m_FieldList.push_back( nameField );

  for( unsigned int i = 0;
       i < m_ObjectDocumentList[index]->GetNumberOfTransformNames(); ++i )
    {
    std::stringstream label;
    label << "Transform" << i;
    FieldType * const transformNameField = new FieldType();
    MET_InitWriteField( transformNameField, label.str().c_str(), MET_STRING,
      m_ObjectDocumentList[index]->GetTransformNames()[i].length(),
      m_ObjectDocumentList[index]->GetTransformNames()[i].c_str() );
    m_FieldList.push_back( transformNameField );
  }

  FieldType * const endObjectField = new FieldType();
  MET_InitWriteField( endObjectField, "EndObject", MET_NONE );
  m_FieldList.push_back( endObjectField );
}


void
MetaObjectDocument
::SetupReadFields( void )
{
  this->Superclass::SetupReadFields();

  FieldType * const numberOfObjectDocumentsField = new FieldType();
  MET_InitReadField( numberOfObjectDocumentsField, "NumberOfObjects",
    MET_INT, true );
  numberOfObjectDocumentsField->required = true;
  numberOfObjectDocumentsField->terminateRead = true;
  m_FieldList.push_back( numberOfObjectDocumentsField );
}


void
MetaObjectDocument
::SetupWriteFields( void )
{
  this->Superclass::SetupWriteFields();

  FieldType * const numberOfObjectDocumentsField = new FieldType();
  MET_InitWriteField( numberOfObjectDocumentsField, "NumberOfObjects",
    MET_INT, m_NumberOfObjectDocuments );
  numberOfObjectDocumentsField->required = true;
  m_FieldList.push_back( numberOfObjectDocumentsField );

  for( unsigned int i = 0; i < m_ObjectDocumentList.size(); ++i )
    {
    this->SetupObjectWriteFields( i );
    }
}


bool
MetaObjectDocument
::WriteFields( void )
{
  if( !this->Superclass::WriteFields() )
    {
    return false;
    }

  std::ofstream fp;

  fp.open( this->GetFileName() );

  if( !fp.is_open() )
    {
    return false;
    }

  if( !MET_Write( fp, &m_FieldList ) )
    {
    return false;
    }

  fp.close();

  return true;
}


void
MetaObjectDocument
::PrintSelf( std::ostream & os, Indent indent ) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "MaximumNumberOfTransforms: "
    << m_MaximumNumberOfTransforms << std::endl;
  os << indent << "NumberOfObjectDocuments:   "
    << m_NumberOfObjectDocuments << std::endl;
  os << indent << "ObjectDocumentList:" << std::endl;

  for( ObjectDocumentListType::const_iterator it =
    m_ObjectDocumentList.begin();
    it != m_ObjectDocumentList.end(); ++it )
    {
    os << indent << *it << std::endl;
    }
}

} // End namespace tube
