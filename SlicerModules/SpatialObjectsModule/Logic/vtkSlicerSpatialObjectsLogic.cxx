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

#include "vtkSlicerSpatialObjectsLogic.h"

// MRML includes
#include <vtkMRMLConfigure.h>
#include <vtkMRMLScene.h>
#include "vtkMRMLSpatialObjectsDisplayPropertiesNode.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsStorageNode.h"
#include "vtkMRMLSpatialObjectsLineDisplayNode.h"
#include "vtkMRMLSpatialObjectsTubeDisplayNode.h"
#include "vtkMRMLSpatialObjectsGlyphDisplayNode.h"

// VTK includes
#include <vtkNew.h>
#include <vtkObjectFactory.h>

// ITK includes
#include <itksys/Directory.hxx>
#include <itksys/SystemTools.hxx>

vtkStandardNewMacro(vtkSlicerSpatialObjectsLogic);

//------------------------------------------------------------------------------
vtkSlicerSpatialObjectsLogic::vtkSlicerSpatialObjectsLogic( void )
{}

//------------------------------------------------------------------------------
vtkSlicerSpatialObjectsLogic::~vtkSlicerSpatialObjectsLogic( void )
{}

//------------------------------------------------------------------------------
int vtkSlicerSpatialObjectsLogic::AddSpatialObjects(const char* dirname,
                                                    const char* suffix )
{
  std::string ssuf = suffix;
  itksys::Directory dir;
  dir.Load(dirname);

  int nfiles = dir.GetNumberOfFiles();
  int res = 1;
  for(int i = 0; i < nfiles; ++i) {
    const char* filename = dir.GetFile(i);
    std::string sname = filename;
    if(!itksys::SystemTools::FileIsDirectory(filename))
      {
      if(sname.find(ssuf) != std::string::npos)
        {
        std::string fullPath = std::string(dir.GetPath()) + "/" + filename;
        if(this->AddSpatialObject(fullPath.c_str()) == NULL)
          {
          res = 0;
          }
        }
      }
  }

  return res;
}

//------------------------------------------------------------------------------
int vtkSlicerSpatialObjectsLogic::
AddSpatialObjects(const char* dirname, std::vector<std::string> suffix)
{
  itksys::Directory dir;
  dir.Load(dirname);

  int nfiles = dir.GetNumberOfFiles();
  int res = 1;

  for(int i = 0; i < nfiles; ++i) {
    const char* filename = dir.GetFile(i);
    std::string name = filename;

    if(!itksys::SystemTools::FileIsDirectory(filename))
      {
      for(unsigned int s = 0; s < suffix.size(); ++s)
        {
        std::string ssuf = suffix[s];
        if(name.find(ssuf) != std::string::npos)
          {
          std::string fullPath = std::string(dir.GetPath()) + "/" + filename;

          if(this->AddSpatialObject(fullPath.c_str()) == NULL)
            {
            res = 0;
            }
          }
        }
      }
  }

  return res;
}

//------------------------------------------------------------------------------
void vtkSlicerSpatialObjectsLogic
::SetSpatialObject(vtkMRMLSpatialObjectsNode* spatialObjectNode,
                   const char* filename)
{
  if (!spatialObjectNode)
    {
    return;
    }

  vtkMRMLSpatialObjectsStorageNode* storageNode =
    vtkMRMLSpatialObjectsStorageNode::SafeDownCast(
      spatialObjectNode->GetStorageNode());
  bool addStorageNodeToScene = false;
  if (!storageNode)
    {
    addStorageNodeToScene = true;
    storageNode = vtkMRMLSpatialObjectsStorageNode::New();
    }
  if (filename)
    {
    storageNode->SetFileName(filename);

    if(storageNode->ReadData(spatialObjectNode) != 0)
      {
      const itksys_stl::string fname(filename);
      itksys_stl::string name =
        itksys::SystemTools::GetFilenameWithoutExtension(fname);
      std::string uname(
        this->GetMRMLScene()->GetUniqueNameByString(name.c_str()));
      spatialObjectNode->SetName(uname.c_str());
      }
    }
  else
    {
    spatialObjectNode->Reset();
    }

  this->GetMRMLScene()->SaveStateForUndo();
  if (storageNode && addStorageNodeToScene)
    {
    this->GetMRMLScene()->AddNode(storageNode);
    storageNode->Delete();
    }
  this->Modified();
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode*
vtkSlicerSpatialObjectsLogic::AddSpatialObject(const char* filename)
{
  vtkDebugMacro("Adding spatial objects from filename ");

  vtkNew<vtkMRMLSpatialObjectsNode> spatialObjectsNode;
  this->SetSpatialObject(spatialObjectsNode.GetPointer(), filename);

  this->GetMRMLScene()->SaveStateForUndo();
  this->GetMRMLScene()->AddNode(spatialObjectsNode.GetPointer());
  this->AddDisplayNodes(spatialObjectsNode.GetPointer());
  this->Modified();

  return spatialObjectsNode.GetPointer();
}


//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode*
vtkSlicerSpatialObjectsLogic::MergeSpatialObjectFromFilename(
  vtkMRMLSpatialObjectsNode* recipient, const char* filename)
{
  if (!filename)
    {
    return recipient;
    }

  vtkNew<vtkMRMLSpatialObjectsNode> donor;
  vtkNew<vtkMRMLSpatialObjectsStorageNode> donorReader;
  donorReader->SetFileName(filename);
  donorReader->ReadData(donor.GetPointer());
  this->SetSpatialObject(donor.GetPointer(), filename);

  this->MergeSpatialObject(recipient, donor.GetPointer());
  return recipient;
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode*
vtkSlicerSpatialObjectsLogic::MergeSpatialObject(
  vtkMRMLSpatialObjectsNode* recipient, vtkMRMLSpatialObjectsNode* donor)
{
  if (!recipient || !donor)
    {
    return recipient;
    }

  recipient->GetSpatialObject()->AddSpatialObject(donor->GetSpatialObject());
  recipient->UpdatePolyDataFromSpatialObject();
  return recipient;
}

//------------------------------------------------------------------------------
void vtkSlicerSpatialObjectsLogic::
AddDisplayNodes(vtkMRMLSpatialObjectsNode* spatialObjectsNode)
{
  if (!spatialObjectsNode)
    {
    return;
    }

  vtkNew<vtkMRMLSpatialObjectsLineDisplayNode> displayLineNode;
  vtkNew<vtkMRMLSpatialObjectsTubeDisplayNode> displayTubeNode;
  vtkNew<vtkMRMLSpatialObjectsGlyphDisplayNode> displayGlyphNode;

  vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> lineProperties;
  vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> tubeProperties;
  vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> glyphProperties;
  glyphProperties->SetGlyphGeometry(
    vtkMRMLSpatialObjectsDisplayPropertiesNode::Lines);

  displayLineNode->SetVisibility(1);
  displayTubeNode->SetVisibility(0);
  displayGlyphNode->SetVisibility(0);

  this->GetMRMLScene()->SaveStateForUndo();
  this->GetMRMLScene()->AddNode(lineProperties.GetPointer());
  this->GetMRMLScene()->AddNode(tubeProperties.GetPointer());
  this->GetMRMLScene()->AddNode(glyphProperties.GetPointer());

  displayLineNode->
    SetAndObserveSpatialObjectsDisplayPropertiesNodeID(
      lineProperties->GetID());
  displayTubeNode->
    SetAndObserveSpatialObjectsDisplayPropertiesNodeID(
      tubeProperties->GetID());
  displayGlyphNode->
    SetAndObserveSpatialObjectsDisplayPropertiesNodeID(
      glyphProperties->GetID());

  this->GetMRMLScene()->AddNode(displayLineNode.GetPointer());
  this->GetMRMLScene()->AddNode(displayTubeNode.GetPointer());
  this->GetMRMLScene()->AddNode(displayGlyphNode.GetPointer());

  displayLineNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");
  displayTubeNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");
  displayGlyphNode->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");

  spatialObjectsNode->SetAndObserveDisplayNodeID(displayLineNode->GetID());
  spatialObjectsNode->AddAndObserveDisplayNodeID(displayTubeNode->GetID());
  spatialObjectsNode->AddAndObserveDisplayNodeID(displayGlyphNode->GetID());
}

//------------------------------------------------------------------------------
int vtkSlicerSpatialObjectsLogic::
SaveSpatialObject(const char* filename,
                  vtkMRMLSpatialObjectsNode *spatialObjectsNode)
{
   if(spatialObjectsNode == NULL || filename == NULL)
    {
    return 0;
    }

  vtkMRMLSpatialObjectsStorageNode* storageNode = NULL;
  vtkMRMLStorageNode* snode = spatialObjectsNode->GetStorageNode();
  if(snode != NULL)
    {
    storageNode = vtkMRMLSpatialObjectsStorageNode::SafeDownCast(snode);
    }

  if(storageNode == NULL)
    {
    storageNode = vtkMRMLSpatialObjectsStorageNode::New();
    storageNode->SetScene(this->GetMRMLScene());
    this->GetMRMLScene()->AddNode(storageNode);
    spatialObjectsNode->SetAndObserveStorageNodeID(storageNode->GetID());
    storageNode->Delete();
    }

  storageNode->SetFileName(filename);

  return storageNode->WriteData(spatialObjectsNode);
}

//------------------------------------------------------------------------------
void vtkSlicerSpatialObjectsLogic
::CopySpatialObject(vtkMRMLSpatialObjectsNode* from,
                    vtkMRMLSpatialObjectsNode* to)
{
  if (from && to)
    {
    to->SetSpatialObject(from->GetSpatialObject());
    }
}

//------------------------------------------------------------------------------
void vtkSlicerSpatialObjectsLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->vtkObject::PrintSelf(os, indent);

  os << indent << "vtkSlicerSpatialObjectsLogic: "
     << this->GetClassName() << "\n";
}

//------------------------------------------------------------------------------
void vtkSlicerSpatialObjectsLogic::RegisterNodes( void )
{
  if(!this->GetMRMLScene())
    {
    return;
    }

  this->GetMRMLScene()->RegisterNodeClass(
    vtkNew<vtkMRMLSpatialObjectsNode>().GetPointer());
  this->GetMRMLScene()->RegisterNodeClass(
    vtkNew<vtkMRMLSpatialObjectsLineDisplayNode>().GetPointer());
  this->GetMRMLScene()->RegisterNodeClass(
    vtkNew<vtkMRMLSpatialObjectsTubeDisplayNode>().GetPointer());
  this->GetMRMLScene()->RegisterNodeClass(
    vtkNew<vtkMRMLSpatialObjectsGlyphDisplayNode>().GetPointer());
  this->GetMRMLScene()->RegisterNodeClass(
    vtkNew<vtkMRMLSpatialObjectsStorageNode>().GetPointer());
}
