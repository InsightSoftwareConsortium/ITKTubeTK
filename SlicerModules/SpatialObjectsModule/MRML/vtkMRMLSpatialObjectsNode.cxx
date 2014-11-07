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
  this->PrepareCleaning();
  this->Reset();
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::~vtkMRMLSpatialObjectsNode( void )
{
  this->RemoveCleaning();
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

  //const char* attName;
  //const char* attValue;
  //while(*atts != NULL)
  //  {
  //  attName = *(atts++);
  //  attValue = *(atts++);
  //  }

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
    this->SetSpatialObject(node->GetSpatialObject());
    }

  this->EndModify(disabledModify);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::Reset()
{
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
void vtkMRMLSpatialObjectsNode::UpdateReferences( void )
{
  for(int ii = 0; ii < this->GetNumberOfDisplayNodes(); ++ii)
    {
    vtkMRMLSpatialObjectsDisplayNode *node = vtkMRMLSpatialObjectsDisplayNode::
      SafeDownCast(this->GetNthDisplayNode(ii));
    if(node)
      {
#if VTK_MAJOR_VERSION <= 5
      node->SetInputPolyData(this->GetFilteredPolyData());
#else
      node->SetInputPolyDataConnection(this->GetFilteredPolyDataConnection());
#endif
      }
    }

  Superclass::UpdateReferences();
}

#if VTK_MAJOR_VERSION <= 5
//------------------------------------------------------------------------------
vtkPolyData* vtkMRMLSpatialObjectsNode::GetFilteredPolyData( void )
{
  this->CleanPolyData->Update();
  return this->CleanPolyData->GetOutput();
}
#else
//------------------------------------------------------------------------------
vtkAlgorithmOutput* vtkMRMLSpatialObjectsNode::
GetFilteredPolyDataConnection( void )
{
  this->CleanPolyData->Update();
  return this->CleanPolyData->GetOutputPort();
}
#endif

namespace
{
//------------------------------------------------------------------------------
// Helper method for factorizing the GetDisplayXNode methods
// (X = Line, Tube or Glyph)
template<typename T>
vtkMRMLSpatialObjectsDisplayNode*
  TemplatedGetDisplayNode(vtkMRMLSpatialObjectsNode* self)
{
  int nnodes = self->GetNumberOfDisplayNodes();
  T *node = NULL;

  for(int n = 0; n < nnodes; ++n)
    {
    node = T::SafeDownCast(self->GetNthDisplayNode(n));
    if(node)
      {
      break;
      }
    }

  return node;
}

//------------------------------------------------------------------------------
// Helper method for factorizing the AddDisplayXNode methods
// (X = Line, Tube or Glyph)
template<typename T>
vtkMRMLSpatialObjectsDisplayNode*
  TemplatedAddDisplayNode(vtkMRMLSpatialObjectsNode* self)
{
  vtkMRMLSpatialObjectsDisplayNode* node = TemplatedGetDisplayNode<T>(self);
  if(node == NULL && self->GetScene())
    {
    node = T::New(); // No smart pointer, it's returned at the end
    vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> properties;
    self->GetScene()->SaveStateForUndo();

    self->GetScene()->AddNode(properties.GetPointer());
    node->SetAndObserveSpatialObjectsDisplayPropertiesNodeID(
      properties->GetID());

    self->GetScene()->AddNode(node);
    node->Delete(); // See vtkMRMLScene::AddNode
    node->SetAndObserveColorNodeID("vtkMRMLColorTableNodeRainbow");
    self->AddAndObserveDisplayNodeID(node->GetID());
    }

  return node;
}
} // end namespace

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetLineDisplayNode( void )
{
  return TemplatedGetDisplayNode<vtkMRMLSpatialObjectsLineDisplayNode>(this);
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetTubeDisplayNode( void )
{
  return TemplatedGetDisplayNode<vtkMRMLSpatialObjectsTubeDisplayNode>(this);
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
GetGlyphDisplayNode( void )
{
  return TemplatedGetDisplayNode<vtkMRMLSpatialObjectsGlyphDisplayNode>(this);
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddLineDisplayNode( void )
{
  return TemplatedAddDisplayNode<vtkMRMLSpatialObjectsLineDisplayNode>(this);
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddTubeDisplayNode( void )
{
  return TemplatedAddDisplayNode<vtkMRMLSpatialObjectsTubeDisplayNode>(this);
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddGlyphDisplayNode( void )
{
  return TemplatedAddDisplayNode<vtkMRMLSpatialObjectsGlyphDisplayNode>(this);
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

  this->UpdateCleaning();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::PrepareCleaning( void )
{
  this->ShuffledIds = vtkIdTypeArray::New();
  this->CleanPolyData = vtkCleanPolyData::New();
  this->CleanPolyData->ConvertLinesToPointsOff();
  this->CleanPolyData->ConvertPolysToLinesOff();
  this->CleanPolyData->ConvertStripsToPolysOff();
  this->CleanPolyData->PointMergingOff();

#if VTK_MAJOR_VERSION <= 5
  this->CleanPolyData->SetInput(this->PolyData);
#else
  this->CleanPolyData->SetInputConnection(this->GetPolyDataConnection());
#endif
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdateCleaning( void )
{
  if(!this->GetPolyData())
    {
    return;
    }

  vtkDebugMacro(<< this->GetClassName() << "Updating the subsampling");

  vtkMRMLSpatialObjectsDisplayNode *node = this->GetLineDisplayNode();
  if(node != NULL)
    {
#if VTK_MAJOR_VERSION <= 5
    node->SetInputPolyData(this->GetFilteredPolyData());
#else
    node->SetInputPolyDataConnection(this->GetFilteredPolyDataConnection());
#endif
    }

  node = this->GetTubeDisplayNode();
  if(node != NULL)
    {
#if VTK_MAJOR_VERSION <= 5
    node->SetInputPolyData(this->GetFilteredPolyData());
#else
    node->SetInputPolyDataConnection(this->GetFilteredPolyDataConnection());
#endif
    }
  node = this->GetGlyphDisplayNode();
  if(node != NULL)
    {
#if VTK_MAJOR_VERSION <= 5
    node->SetInputPolyData(this->GetFilteredPolyData());
#else
    node->SetInputPolyDataConnection(this->GetFilteredPolyDataConnection());
#endif
    }

  this->InvokeEvent(vtkMRMLModelNode::PolyDataModifiedEvent, this);
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::RemoveCleaning( void )
{
  this->CleanPolyData->Delete();
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

    double spacingX = currTube->GetSpacing()[0];
    if( spacingX == 0 )
      {
      spacingX = 1;
      }

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
      tubeRadius->SetTuple1(pointID, tubePoint->GetRadius() * spacingX );

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
  sodn->SetVisibility(1);
  sodn = this->AddGlyphDisplayNode();
  sodn->GetSpatialObjectsDisplayPropertiesNode()->SetGlyphGeometry(
    vtkMRMLSpatialObjectsDisplayPropertiesNode::Lines);
  sodn->SetVisibility(0);
}
