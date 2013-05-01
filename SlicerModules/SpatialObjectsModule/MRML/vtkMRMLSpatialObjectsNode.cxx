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

// VTK includes
#include <vtkCleanPolyData.h>
#include <vtkCommand.h>
#include <vtkEventBroker.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>

// TractographyMRML includes
#include "vtkMRMLSpatialObjectsGlyphDisplayNode.h"
#include "vtkMRMLSpatialObjectsLineDisplayNode.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsStorageNode.h"
#include "vtkMRMLSpatialObjectsTubeDisplayNode.h"

// MRML includes
#include <vtkMRMLSpatialObjectsDisplayPropertiesNode.h>
#include <vtkMRMLScene.h>

// STD includes
#include <math.h>
#include <vector>
#include <algorithm>

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSpatialObjectsNode);


//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::vtkMRMLSpatialObjectsNode( void )
{
  this->PrepareSubsampling();
  this->SubsamplingRatio = 1;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::~vtkMRMLSpatialObjectsNode( void )
{
  this->CleanSubsampling();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::WriteXML(ostream& of, int nIndent)
{
  Superclass::WriteXML(of, nIndent);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::ReadXMLAttributes(const char** atts)
{
  int disabledModify = this->StartModify();

  Superclass::ReadXMLAttributes(atts);

  const char* attName;
  const char* attValue;
  while (*atts != NULL)
    {
    attName = *(atts++);
    attValue = *(atts++);

    if (!strcmp(attName, "SubsamplingRatio"))
      {
      this->SubsamplingRatio = atof(attValue);
      }
    }

  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::Copy(vtkMRMLNode *anode)
{
  int disabledModify = this->StartModify();

  Superclass::Copy(anode);

  vtkMRMLSpatialObjectsNode *node =
    vtkMRMLSpatialObjectsNode::SafeDownCast(anode);

  if (node)
    {
    this->SetSubsamplingRatio(node->SubsamplingRatio);
    }

  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdateScene(vtkMRMLScene *scene)
{
  Superclass::UpdateScene(scene);

  int disabledModify = this->StartModify();

  // We are forcing the update of the fields as UpdateScene
  // should only be called after loading data
  double ActualSubsamplingRatio = this->SubsamplingRatio;
  this->SubsamplingRatio = 0.;
  this->SetSubsamplingRatio(ActualSubsamplingRatio);

  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdateReferences( void )
{
  for(int ii = 0; ii < this->GetNumberOfDisplayNodes(); ++ii)
    {
    vtkMRMLSpatialObjectsDisplayNode *node = vtkMRMLSpatialObjectsDisplayNode::
      SafeDownCast(this->GetNthDisplayNode(ii));
    if (node)
      {
      node->SetInputPolyData(this->GetFilteredPolyData());
      }
    }

  Superclass::UpdateReferences();
}

//------------------------------------------------------------------------------
vtkPolyData* vtkMRMLSpatialObjectsNode::GetFilteredPolyData( void )
{
  return this->PolyData;
  //return this->CleanPolyDataPostSubsampling->GetOutput();
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetLineDisplayNode( void )
{
  int nnodes = this->GetNumberOfDisplayNodes();
  vtkMRMLSpatialObjectsLineDisplayNode *node = NULL;

  for (int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsLineDisplayNode::SafeDownCast(
             this->GetNthDisplayNode(n));
    if (node)
      {
      break;
      }
    }

  return node;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetTubeDisplayNode( void )
{
  int nnodes = this->GetNumberOfDisplayNodes();
  vtkMRMLSpatialObjectsTubeDisplayNode *node = NULL;

  for (int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsTubeDisplayNode::SafeDownCast(
             this->GetNthDisplayNode(n));
    if (node)
      {
      break;
      }
    }

  return node;
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetGlyphDisplayNode( void )
{
  int nnodes = this->GetNumberOfDisplayNodes();
  vtkMRMLSpatialObjectsGlyphDisplayNode *node = NULL;

  for (int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsGlyphDisplayNode::SafeDownCast(
            this->GetNthDisplayNode(n));
    if (node)
      {
      break;
      }
    }

  return node;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddLineDisplayNode( void )
{
  vtkMRMLSpatialObjectsDisplayNode *node = this->GetLineDisplayNode();
  if (node == NULL)
    {
    node = vtkMRMLSpatialObjectsLineDisplayNode::New();

    if (this->GetScene())
      {
      this->GetScene()->AddNode(node);
      node->Delete();

      vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> glyphSOPN;
      this->GetScene()->AddNode(glyphSOPN.GetPointer());
      node->
        SetAndObserveSpatialObjectsDisplayPropertiesNodeID(glyphSOPN->GetID());
      node->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");

      this->AddAndObserveDisplayNodeID(node->GetID());
      node->SetInputPolyData(this->GetFilteredPolyData());
      }
    }

  return node;
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddTubeDisplayNode( void )
{
  vtkMRMLSpatialObjectsDisplayNode *node = this->GetTubeDisplayNode();
  if (node == NULL)
    {
    node = vtkMRMLSpatialObjectsTubeDisplayNode::New();
    if (this->GetScene())
      {
      this->GetScene()->AddNode(node);
      node->Delete();

      vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> glyphSOPN;
      this->GetScene()->AddNode(glyphSOPN.GetPointer());
      node->
        SetAndObserveSpatialObjectsDisplayPropertiesNodeID(glyphSOPN->GetID());
      node->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");

      this->AddAndObserveDisplayNodeID(node->GetID());
      node->SetInputPolyData(this->GetFilteredPolyData());
      }
    }

  return node;
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddGlyphDisplayNode( void )
{
  vtkMRMLSpatialObjectsDisplayNode *node = this->GetGlyphDisplayNode();
  if (node == NULL)
    {
    node = vtkMRMLSpatialObjectsGlyphDisplayNode::New();
    if (this->GetScene())
      {
      this->GetScene()->AddNode(node);
      node->Delete();

      vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> glyphSOPN;
      this->GetScene()->AddNode(glyphSOPN.GetPointer());
      node->
        SetAndObserveSpatialObjectsDisplayPropertiesNodeID(glyphSOPN->GetID());
      node->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");

      this->AddAndObserveDisplayNodeID(node->GetID());
      node->SetInputPolyData(this->GetFilteredPolyData());
      }
    }

  return node;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::SetAndObservePolyData(vtkPolyData* polyData)
{
  vtkMRMLModelNode::SetAndObservePolyData(polyData);

  if (!polyData)
    {
    return;
    }

  const vtkIdType numberOfPairs = polyData->GetNumberOfLines();

  std::vector<vtkIdType> idVector;
  for(vtkIdType i = 0; i < numberOfPairs; ++i )
    {
    idVector.push_back(i);
    }

  random_shuffle(idVector.begin(), idVector.end());

  this->ShuffledIds->Initialize();
  this->ShuffledIds->SetNumberOfTuples(numberOfPairs);
  for(vtkIdType i = 0;  i < numberOfPairs; ++i)
    {
    this->ShuffledIds->SetValue(i, idVector[i]);
    }

  float subsamplingRatio = 1.f;
  this->SetSubsamplingRatio(subsamplingRatio);
  this->UpdateSubsampling();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::SetSubsamplingRatio(float ratio)
{
  vtkDebugMacro(<< this->GetClassName()
                << " (" << this << "): setting subsamplingRatio to " << ratio);

  const float oldSubsampling = this->SubsamplingRatio;
  // Clamp
  const float newSubsamplingRatio =
    (ratio < 0. ? 0. : (ratio > 1. ? 1.: ratio));
  if (oldSubsampling != newSubsamplingRatio)
    {
    this->SubsamplingRatio = newSubsamplingRatio;
    this->UpdateSubsampling();
    this->Modified();
    }
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::PrepareSubsampling( void )
{
  this->ShuffledIds = vtkIdTypeArray::New();
  this->CleanPolyDataPostSubsampling = vtkCleanPolyData::New();
  this->CleanPolyDataPostSubsampling->ConvertLinesToPointsOff();
  this->CleanPolyDataPostSubsampling->ConvertPolysToLinesOff();
  this->CleanPolyDataPostSubsampling->ConvertStripsToPolysOff();
  this->CleanPolyDataPostSubsampling->PointMergingOff();

  this->CleanPolyDataPostSubsampling->SetInput(this->PolyData);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdateSubsampling( void )
{
  if (!this->GetPolyData())
    {
    return;
    }

  vtkDebugMacro(<< this->GetClassName() << "Updating the subsampling");

  vtkMRMLSpatialObjectsDisplayNode *node = this->GetLineDisplayNode();
  if (node != NULL)
    {
    node->SetInputPolyData(this->GetFilteredPolyData());
    }

  node = this->GetTubeDisplayNode();
  if (node != NULL)
    {
    node->SetInputPolyData(this->GetFilteredPolyData());
    }
  node = this->GetGlyphDisplayNode();
  if (node != NULL)
    {
    node->SetInputPolyData(this->GetFilteredPolyData());
    }

  this->InvokeEvent(vtkMRMLModelNode::PolyDataModifiedEvent, this);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::CleanSubsampling( void )
{
  this->CleanPolyDataPostSubsampling->Delete();
  this->ShuffledIds->Delete();
}

//------------------------------------------------------------------------------
vtkMRMLStorageNode* vtkMRMLSpatialObjectsNode::CreateDefaultStorageNode( void )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::CreateDefaultStorageNode");

  return
    vtkMRMLStorageNode::SafeDownCast(vtkMRMLSpatialObjectsStorageNode::New());
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::CreateDefaultDisplayNodes( void )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::CreateDefaultDisplayNodes");

  vtkMRMLSpatialObjectsDisplayNode *sodn = this->AddLineDisplayNode();
  sodn->SetVisibility(1);
  sodn = this->AddTubeDisplayNode();
  sodn->SetVisibility(0);
  sodn = this->AddGlyphDisplayNode();
  sodn->SetVisibility(0);
}
