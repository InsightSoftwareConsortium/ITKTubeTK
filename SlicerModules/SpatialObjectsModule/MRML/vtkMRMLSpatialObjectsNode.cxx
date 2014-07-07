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
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkCommand.h>
#include <vtkDoubleArray.h>
#include <vtkEventBroker.h>
#include <vtkExtractPolyDataGeometry.h>
#include <vtkExtractSelectedPolyDataIds.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPlanes.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
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
#include <cmath>
#include <vector>
#include <algorithm>

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSpatialObjectsNode);


//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::vtkMRMLSpatialObjectsNode( void )
{
  this->PrepareSubsampling();
  this->Reset();
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
  while(*atts != NULL)
    {
    attName = *(atts++);
    attValue = *(atts++);

    if(!std::strcmp(attName, "SubsamplingRatio"))
      {
      this->SubsamplingRatio = std::atof(attValue);
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

  if(node)
    {
    this->SetSubsamplingRatio(node->SubsamplingRatio);
    }

  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::Reset()
{
  this->SubsamplingRatio = 1;

  this->SpatialObject = TubeNetType::New();
  this->SpatialObject->Initialize();
  this->UpdatePolyDataFromSpatialObject();
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::TubeNetPointerType
vtkMRMLSpatialObjectsNode::GetSpatialObject()
{
  return this->SpatialObject;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::SetSpatialObject(TubeNetPointerType object)
{
  if (this->SpatialObject == object)
    {
    return;
    }

  this->SpatialObject = object;
  this->UpdatePolyDataFromSpatialObject();
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
    if(node)
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

  for(int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsLineDisplayNode::SafeDownCast(
             this->GetNthDisplayNode(n));
    if(node)
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

  for(int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsTubeDisplayNode::SafeDownCast(
             this->GetNthDisplayNode(n));
    if(node)
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

  for(int n = 0; n < nnodes; ++n)
    {
    node = vtkMRMLSpatialObjectsGlyphDisplayNode::SafeDownCast(
            this->GetNthDisplayNode(n));
    if(node)
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
  if(node == NULL)
    {
    node = vtkMRMLSpatialObjectsLineDisplayNode::New();

    if(this->GetScene())
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
  if(node == NULL)
    {
    node = vtkMRMLSpatialObjectsTubeDisplayNode::New();
    if(this->GetScene())
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
  if(node == NULL)
    {
    node = vtkMRMLSpatialObjectsGlyphDisplayNode::New();
    if(this->GetScene())
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

  if(!polyData)
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
  if(oldSubsampling != newSubsamplingRatio)
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
  if(!this->GetPolyData())
    {
    return;
    }

  vtkDebugMacro(<< this->GetClassName() << "Updating the subsampling");

  vtkMRMLSpatialObjectsDisplayNode *node = this->GetLineDisplayNode();
  if(node != NULL)
    {
    node->SetInputPolyData(this->GetFilteredPolyData());
    }

  node = this->GetTubeDisplayNode();
  if(node != NULL)
    {
    node->SetInputPolyData(this->GetFilteredPolyData());
    }
  node = this->GetGlyphDisplayNode();
  if(node != NULL)
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
void vtkMRMLSpatialObjectsNode::UpdatePolyDataFromSpatialObject( void )
{
  typedef itk::Point<double, 3>                PointType;
  typedef itk::VesselTubeSpatialObject<3>      VesselTubeType;
  typedef VesselTubeType::TubePointType        VesselTubePointType;

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    this->SpatialObject->GetChildren(
      this->SpatialObject->GetMaximumDepth(), childName);

  // -----------------------------------------------------------------------
  // Copy skeleton points from vessels into polydata structure
  // -----------------------------------------------------------------------

  // Initialize the SpatialObject
  // Count number of points && remove dupplicate
  int totalNumberOfPoints = 0;
  for(TubeNetType::ChildrenListType::iterator tubeIT = tubeList->begin();
        tubeIT != tubeList->end();
        ++tubeIT )
    {
    VesselTubeType* currTube =
      dynamic_cast<VesselTubeType*>((*tubeIT).GetPointer());
    if (!currTube)
      {
      continue;
      }

    currTube->RemoveDuplicatePoints();
    const itk::SizeValueType numberOfPoints = currTube->GetNumberOfPoints();
    if( numberOfPoints < 2 )
      {
      continue;
      }

    totalNumberOfPoints += numberOfPoints;
    }

  // Create the points
  vtkNew<vtkPoints> vesselsPoints;
  vesselsPoints->SetNumberOfPoints(totalNumberOfPoints);

  // Create the Lines
  vtkNew<vtkCellArray> vesselLinesCA;

  // Create scalar array that indicates the radius at each
  // centerline point.
  vtkNew<vtkDoubleArray> tubeRadius;
  tubeRadius->SetName("TubeRadius");
  tubeRadius->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates TubeID.
  vtkNew<vtkDoubleArray> tubeIDs;
  tubeIDs->SetName("TubeIDs");
  tubeIDs->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates both tangents at each
  // centerline point.
  vtkNew<vtkDoubleArray> tan1;
  tan1->SetName("Tan1");
  tan1->SetNumberOfTuples(3 * totalNumberOfPoints);
  tan1->SetNumberOfComponents(3);

  vtkNew<vtkDoubleArray> tan2;
  tan2->SetName("Tan2");
  tan2->SetNumberOfTuples(3 * totalNumberOfPoints);
  tan2->SetNumberOfComponents(3);

  // Create scalar array that indicates Ridgness and medialness at each
  // centerline point.
  bool containsMidialnessInfo = false;
  vtkNew<vtkDoubleArray> medialness;
  medialness->SetName("Medialness");
  medialness->SetNumberOfTuples(totalNumberOfPoints);

  bool containsRidgnessInfo = false;
  vtkNew<vtkDoubleArray> ridgeness;
  ridgeness->SetName("Ridgeness");
  ridgeness->SetNumberOfTuples(totalNumberOfPoints);

  int pointID = 0;
  for(TubeNetType::ChildrenListType::iterator tubeIT = tubeList->begin();
        tubeIT != tubeList->end(); ++tubeIT )
    {
    VesselTubeType* currTube =
      dynamic_cast<VesselTubeType*>((*tubeIT).GetPointer());
    if (!currTube)
      {
      continue;
      }

    const itk::SizeValueType tubeSize = currTube->GetNumberOfPoints();
    if( tubeSize < 2 )
      {
      continue;
      }

    currTube->ComputeTangentAndNormals();

    // Create a pointID list [linear for a polyline]
    vtkIdType* pointIDs = new vtkIdType[tubeSize];
    vtkNew<vtkPolyLine> vesselLine;

    // Get the tube element spacing information.
    //const double* axesRatio = currTube->GetSpacing();

    currTube->ComputeObjectToWorldTransform();

    size_t numberOfPoints = currTube->GetPoints().size();
    for(size_t index = 0; index < numberOfPoints; ++pointID, ++index)
      {
      VesselTubePointType* tubePoint =
        dynamic_cast<VesselTubePointType*>(currTube->GetPoint(index));
      assert(tubePoint);

      PointType inputPoint = tubePoint->GetPosition();

      inputPoint =
        currTube->GetIndexToWorldTransform()->TransformPoint( inputPoint );
      inputPoint[0] = -inputPoint[0];
      inputPoint[1] = -inputPoint[1];

      pointIDs[index] = pointID;

      // Insert points using the element spacing information.
      vesselsPoints->SetPoint( pointID, inputPoint[0], inputPoint[1],
        inputPoint[2] );
        //inputPoint[1] * axesRatio[1] / axesRatio[0],
        //inputPoint[2] * axesRatio[2] / axesRatio[0]);

      // TubeID
      tubeIDs->SetTuple1(pointID, currTube->GetId());

      // Radius
      tubeRadius->SetTuple1(pointID, tubePoint->GetRadius());

      // Tangeantes
      tan1->SetTuple3(pointID,
                      tubePoint->GetNormal1()[0],
                      tubePoint->GetNormal1()[1],
                      tubePoint->GetNormal1()[2]);

      tan2->SetTuple3(pointID,
                      tubePoint->GetNormal2()[0],
                      tubePoint->GetNormal2()[1],
                      tubePoint->GetNormal2()[2]);

      // Medialness & Ridgness
      if(tubePoint->GetMedialness() != 0)
        {
        containsMidialnessInfo = true;
        medialness->SetTuple1(pointID, tubePoint->GetMedialness());
        }

      if(tubePoint->GetRidgeness() != 0)
        {
        containsRidgnessInfo = true;
        ridgeness->SetTuple1(pointID, tubePoint->GetRidgeness());
        }
      }

    vesselLine->Initialize(tubeSize,
                            pointIDs,
                            vesselsPoints.GetPointer());
    vesselLinesCA->InsertNextCell(vesselLine.GetPointer());
    delete [] pointIDs;
    }

  // Convert spatial objects to a PolyData
  vtkNew<vtkPolyData> vesselsPD;
  vesselsPD->SetLines(vesselLinesCA.GetPointer());
  vesselsPD->SetPoints(vesselsPoints.GetPointer());

  // Add the Radius information
  vesselsPD->GetPointData()->AddArray(tubeRadius.GetPointer());
  vesselsPD->GetPointData()->SetActiveScalars("TubeRadius");

  // Add the TudeID information
  vesselsPD->GetPointData()->AddArray(tubeIDs.GetPointer());

  // Add Tangeantes information
  vesselsPD->GetPointData()->AddArray(tan1.GetPointer());
  vesselsPD->GetPointData()->AddArray(tan2.GetPointer());

  // Add Medialness & Ridgness if contains information
  if(containsMidialnessInfo == true)
    {
    vesselsPD->GetPointData()->AddArray(medialness.GetPointer());
    }

  if(containsRidgnessInfo == true)
    {
    vesselsPD->GetPointData()->AddArray(ridgeness.GetPointer());
    }

  // Remove any duplicate points from polydata.
  // The tubes generation will fails if any duplicates points are present.
  // Cleaned before, could create degeneration problems with the cells
  //vtkNew<vtkCleanPolyData> cleanedVesselPD;
  //cleanedVesselPD->SetInput(vesselsPD.GetPointer());

  vtkDebugMacro("Points: " << totalNumberOfPoints);

  this->SetAndObservePolyData(vesselsPD.GetPointer());
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
