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

#include <vtkObjectFactory.h>

// MRML includes
#include "vtkMRMLSpatialObjectsStorageNode.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsDisplayNode.h"

// VTK includes
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStringArray.h>

// ITK includes
#include <itkTubeSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro( vtkMRMLSpatialObjectsStorageNode );

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os,indent );
}

//------------------------------------------------------------------------------
bool vtkMRMLSpatialObjectsStorageNode::
CanReadInReferenceNode( vtkMRMLNode *refNode )
{
  return refNode->IsA( "vtkMRMLSpatialObjectsNode" );
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsStorageNode::ReadDataInternal( vtkMRMLNode *refNode )
{
  vtkMRMLSpatialObjectsNode* spatialObjectsNode =
    vtkMRMLSpatialObjectsNode::SafeDownCast( refNode );

  if( Superclass::ReadDataInternal( refNode ) != 0 )
    {
    return 0;
    }

  std::string fullName = this->GetFullNameFromFileName();
  if( fullName == std::string( "" ) )
    {
    vtkErrorMacro( "ReadData: File name not specified" );
    return 0;
    }

  // compute file prefix
  std::string name( fullName );
  std::string::size_type loc = name.find_last_of( "." );
  if( loc == std::string::npos )
    {
    vtkErrorMacro( "ReadData: no file extension specified: " << name.c_str() );
    return 0;
    }
  std::string extension = name.substr( loc );
  vtkDebugMacro( "ReadData: extension = " << extension.c_str() );

  int result = 1;
  try
  {
    if( extension == std::string( ".tre" ) )
      {
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName( fullName );
      reader->Update();

      spatialObjectsNode->SetSpatialObject( reader->GetGroup() );
      }
  }
  catch( ... )
  {
    result = 0;
  }

  if( spatialObjectsNode->GetPolyData() != NULL )
    {
    // is there an active scalar array?
    if( spatialObjectsNode->GetDisplayNode() )
      {
      double *scalarRange = spatialObjectsNode->GetPolyData()->GetScalarRange();
      if( scalarRange )
        {
        vtkDebugMacro( "ReadData: setting scalar range " << scalarRange[0]
                      << ", " << scalarRange[1] );
        spatialObjectsNode->GetDisplayNode()->SetScalarRange( scalarRange );
        }
      }

    spatialObjectsNode->GetPolyData()->Modified();
    }

  return result;
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsStorageNode::WriteDataInternal( vtkMRMLNode *refNode )
{
  vtkMRMLSpatialObjectsNode* spatialObjects =
    vtkMRMLSpatialObjectsNode::SafeDownCast( refNode );

  const std::string fullName = this->GetFullNameFromFileName();
  if( fullName == std::string( "" ) )
    {
    vtkErrorMacro( "vtkMRMLModelNode: File name not specified" );
    return 0;
    }

  const std::string extension =
    itksys::SystemTools::GetFilenameLastExtension( fullName );

  int result = 1;
  if( extension == ".tre" )
    {
    try
      {
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName( fullName.c_str() );
      writer->SetInput( spatialObjects->GetSpatialObject() );
      writer->Update();
      }
    catch( ... )
      {
      result = 0;
      vtkErrorMacro( "Error occured writing Spatial Objects: "
                    << fullName.c_str() );
      }
    }
  else
    {
    result = 0;
    vtkErrorMacro( << "No file extension recognized: " << fullName.c_str() );
    }

  return result;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::InitializeSupportedReadFileTypes( void )
{
  this->SupportedReadFileTypes->InsertNextValue( "SpatialObject ( .tre )" );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::InitializeSupportedWriteFileTypes( void )
{
  this->SupportedWriteFileTypes->InsertNextValue( "SpatialObject ( .tre )" );
}

//------------------------------------------------------------------------------
const char* vtkMRMLSpatialObjectsStorageNode::GetDefaultWriteFileExtension( void )
{
  return "tre";
}
