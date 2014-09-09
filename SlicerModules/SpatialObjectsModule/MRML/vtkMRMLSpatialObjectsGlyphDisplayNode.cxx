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

#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>

// VTK includes
#include <vtkConeSource.h>
#include <vtkGlyph3DMapper.h>

// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLNode.h>
#include "vtkMRMLSpatialObjectsGlyphDisplayNode.h"
#include <vtkMRMLDiffusionTensorDisplayPropertiesNode.h>

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSpatialObjectsGlyphDisplayNode);

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsGlyphDisplayNode::vtkMRMLSpatialObjectsGlyphDisplayNode( void )
{
  this->Glyph3DMapper = vtkGlyph3DMapper::New();
  this->ColorMode = vtkMRMLSpatialObjectsDisplayNode::colorModeScalar;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsGlyphDisplayNode::~vtkMRMLSpatialObjectsGlyphDisplayNode( void )
{
  this->RemoveObservers(vtkCommand::ModifiedEvent, this->MRMLCallbackCommand);
  this->Glyph3DMapper->Delete();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsGlyphDisplayNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsGlyphDisplayNode::ReadXMLAttributes(const char** atts)
{
  int disabledModify = this->StartModify();
  Superclass::ReadXMLAttributes(atts);
  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsGlyphDisplayNode::Copy(vtkMRMLNode *anode)
{
  int disabledModify = this->StartModify();
  Superclass::Copy(anode);
  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsGlyphDisplayNode::PrintSelf(ostream& os,
                                                      vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------
vtkAlgorithmOutput* vtkMRMLSpatialObjectsGlyphDisplayNode::GetOutputPort( void )
{
  return this->Glyph3DMapper->GetOutputPort();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsGlyphDisplayNode::UpdatePolyDataPipeline( void )
{
  if(!this->GetInputPolyData()|| !this->Visibility)
    {
    return;
    }

  this->Superclass::UpdatePolyDataPipeline();

  /*if(this->Glyph3DMapper)
    {
    this->Glyph3DMapper->SetInputConnection(
      this->PolyData->GetSource()->GetOutputPort());
    }*/

  // Set display properties according to the
  // tensor-specific display properties node for glyphs
  // TODO Create own prooperties nodes

}
