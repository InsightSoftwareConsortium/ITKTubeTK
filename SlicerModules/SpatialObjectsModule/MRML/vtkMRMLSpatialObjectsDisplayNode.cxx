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

// MRML includes
#include <vtkMRMLDisplayableNode.h>
#include <vtkMRMLScene.h>
#include "vtkMRMLSpatialObjectsDisplayNode.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsDisplayPropertiesNode.h"

// VTK includes
#include <vtkCommand.h>

// STD includes
#include <sstream>

//------------------------------------------------------------------------------
vtkCxxSetReferenceStringMacro( vtkMRMLSpatialObjectsDisplayNode,
                              SpatialObjectsDisplayPropertiesNodeID );

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode::vtkMRMLSpatialObjectsDisplayNode( void )
{
  this->BackfaceCulling = 0;

  // Enumerated
  this->ColorMode = this->colorModeSolid;
  this->SetColor( 0.9, 0.2, 0.1 );

  this->SpatialObjectsDisplayPropertiesNode = NULL;
  this->SpatialObjectsDisplayPropertiesNodeID = NULL;

  this->ScalarRange[0] = 0.;
  this->ScalarRange[1] = 1.;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode::~vtkMRMLSpatialObjectsDisplayNode( void )
{
  this->SetAndObserveSpatialObjectsDisplayPropertiesNodeID( NULL );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::WriteXML( ostream& of, int nIndent )
{
  Superclass::WriteXML( of, nIndent );

  vtkIndent indent( nIndent );
  of << indent << " colorMode =\"" << this->ColorMode << "\"";

  if( this->SpatialObjectsDisplayPropertiesNodeID != NULL )
    {
    of << indent << " SpatialObjectsDisplayPropertiesNodeRef=\""
       << this->SpatialObjectsDisplayPropertiesNodeID << "\"";
    }
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::ReadXMLAttributes( const char** atts )
{
  int disabledModify = this->StartModify();

  Superclass::ReadXMLAttributes( atts );

  const char* attName;
  const char* attValue;
  while( *atts != NULL )
    {
    attName = *( atts++ );
    attValue = *( atts++ );

    if( !std::strcmp( attName, "colorMode" ) )
      {
      std::stringstream ss;
      ss << attValue;
      int colorMode;
      ss >> colorMode;
      this->SetColorMode( colorMode );
      }

    else if( !std::strcmp( attName, "SpatialObjectsDisplayPropertiesNodeRef" ) )
      {
      this->SetSpatialObjectsDisplayPropertiesNodeID( attValue );
      }
    }

  this->EndModify( disabledModify );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::Copy( vtkMRMLNode *anode )
{
  int disabledModify = this->StartModify();

 vtkMRMLSpatialObjectsDisplayNode *node =
   vtkMRMLSpatialObjectsDisplayNode::SafeDownCast( anode );
 this->SetColorMode( node->ColorMode );

  Superclass::Copy( anode );

  this->SetSpatialObjectsDisplayPropertiesNodeID( 
    node->SpatialObjectsDisplayPropertiesNodeID );

  this->EndModify( disabledModify );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::PrintSelf( ostream& os, vtkIndent indent )
{
  Superclass::PrintSelf( os,indent );
  os << indent << "ColorMode: " << this->ColorMode << "\n";
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::UpdateScene( vtkMRMLScene *scene )
{
   Superclass::UpdateScene( scene );

   this->SetAndObserveSpatialObjectsDisplayPropertiesNodeID( 
     this->GetSpatialObjectsDisplayPropertiesNodeID() );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::UpdateReferences( void )
{
  Superclass::UpdateReferences();

  if( this->SpatialObjectsDisplayPropertiesNodeID != NULL &&
      this->Scene->GetNodeByID( 
        this->SpatialObjectsDisplayPropertiesNodeID ) == NULL )
    {
    this->SetAndObserveSpatialObjectsDisplayPropertiesNodeID( NULL );
    }
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::UpdateReferenceID( const char *oldID,
                                                         const char *newID )
{
  if( this->SpatialObjectsDisplayPropertiesNodeID &&
      !std::strcmp( oldID, this->SpatialObjectsDisplayPropertiesNodeID ) )
    {
    this->SetSpatialObjectsDisplayPropertiesNodeID( newID );
    }
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayPropertiesNode* vtkMRMLSpatialObjectsDisplayNode::
GetSpatialObjectsDisplayPropertiesNode( void )
{
  vtkMRMLSpatialObjectsDisplayPropertiesNode* node = NULL;

  // Find the node corresponding to the ID we have saved.
  if  ( this->GetScene() && this->GetSpatialObjectsDisplayPropertiesNodeID() )
    {
    vtkMRMLNode* cnode = this->
      GetScene()->GetNodeByID( this->SpatialObjectsDisplayPropertiesNodeID );
    node = vtkMRMLSpatialObjectsDisplayPropertiesNode::SafeDownCast( cnode );
    }

  return node;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsDisplayNode::
SetAndObserveSpatialObjectsDisplayPropertiesNodeID( const char *id )
{
  if( 
      ( id != this->GetSpatialObjectsDisplayPropertiesNodeID() )
      && id != NULL && this->GetSpatialObjectsDisplayPropertiesNodeID() != NULL
      && ( std::strcmp( id, this->GetSpatialObjectsDisplayPropertiesNodeID() ) == 0 ) )
    {
    return;
    }

  // Stop observing any old node
  vtkSetAndObserveMRMLObjectMacro( this->SpatialObjectsDisplayPropertiesNode,
                                  NULL );

  // Set the ID. This is the "ground truth" reference to the node.
  this->SetSpatialObjectsDisplayPropertiesNodeID( id );

  // Get the node corresponding to the ID.
  // This pointer is only to observe the object.
  vtkMRMLNode *cnode = this->GetSpatialObjectsDisplayPropertiesNode();

  // Observe the node using the pointer.
  vtkSetAndObserveMRMLObjectMacro( this->SpatialObjectsDisplayPropertiesNode ,
                                  cnode );

  //The new SpatialObjectsDisplayPropertiesNode can have a different setting
  // on the properties so we emit the event that the polydata has been modified.
  if( cnode )
    {
    this->InvokeEvent( vtkMRMLModelNode::PolyDataModifiedEvent, this );
    }
}

//------------------------------------------------------------------------------
std::vector<int> vtkMRMLSpatialObjectsDisplayNode::GetSupportedColorModes( void )
{
  std::vector<int> modes;

  modes.clear();
  modes.push_back( vtkMRMLSpatialObjectsDisplayPropertiesNode::ColorOrientation );
  modes.push_back( vtkMRMLSpatialObjectsDisplayPropertiesNode::LinearMeasure );
  modes.push_back( 
    vtkMRMLSpatialObjectsDisplayPropertiesNode::RelativeAnisotropy );

  return modes;
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsDisplayNode::GetNumberOfScalarInvariants( void )
{
  static std::vector<int> modes =
    vtkMRMLSpatialObjectsDisplayNode::GetSupportedColorModes();

  return modes.size();
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsDisplayNode::GetNthScalarInvariant( int i )
{
  static std::vector<int> modes =
    vtkMRMLSpatialObjectsDisplayNode::GetSupportedColorModes();

  return modes[i];
}
