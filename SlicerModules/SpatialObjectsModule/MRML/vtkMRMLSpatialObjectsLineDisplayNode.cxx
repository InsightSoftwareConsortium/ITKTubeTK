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

// VTK includes
#include <vtkAssignAttribute.h>
#include <vtkCallbackCommand.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLNode.h>
#include "vtkMRMLSpatialObjectsLineDisplayNode.h"
#include "vtkMRMLSpatialObjectsDisplayPropertiesNode.h"

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro( vtkMRMLSpatialObjectsLineDisplayNode );


//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsLineDisplayNode::vtkMRMLSpatialObjectsLineDisplayNode( void )
{
  this->ColorMode = vtkMRMLSpatialObjectsDisplayNode::colorModeSolid;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsLineDisplayNode::~vtkMRMLSpatialObjectsLineDisplayNode( void )
{
  this->RemoveObservers ( vtkCommand::ModifiedEvent, this->MRMLCallbackCommand );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::WriteXML( ostream& of, int nIndent )
{
  Superclass::WriteXML( of, nIndent );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::ReadXMLAttributes( const char** atts )
{
  Superclass::ReadXMLAttributes( atts );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::Copy( vtkMRMLNode *anode )
{
  Superclass::Copy( anode );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::PrintSelf( ostream& os,
                                                     vtkIndent indent )
{
  Superclass::PrintSelf( os,indent );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsLineDisplayNode::UpdatePolyDataPipeline( void )
{
  if( !this->GetInputPolyData() || !this->Visibility )
    {
    return;
    }

  this->Superclass::UpdatePolyDataPipeline();

  // Set display properties according to the
  // line display properties node
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    SpatialObjectsDisplayPropertiesNode =
      this->GetSpatialObjectsDisplayPropertiesNode();

  if( SpatialObjectsDisplayPropertiesNode != NULL )
    {
    const int colorMode = this->GetColorMode();
    if( colorMode ==
          vtkMRMLSpatialObjectsDisplayNode::colorModeSolid )
      {
      this->ScalarVisibilityOff();

      vtkMRMLNode* colorNode =
        this->GetScene()->GetNodeByID( "vtkMRMLColorTableNodeFullRainbow" );
      if( colorNode )
        {
        this->SetAndObserveColorNodeID( colorNode->GetID() );
        }

      this->AutoScalarRangeOff();
      this->SetScalarRange( 0, 255 );
      }
    else if( colorMode ==
               vtkMRMLSpatialObjectsDisplayNode::colorModeScalarData )
      {
      this->ScalarVisibilityOn();
      this->AssignAttribute->Update();
      }
    }
  else
    {
    this->ScalarVisibilityOff();
    }
}
