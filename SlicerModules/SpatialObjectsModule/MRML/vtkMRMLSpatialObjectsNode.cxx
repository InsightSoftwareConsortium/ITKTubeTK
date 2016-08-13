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
#include <vtkUnsignedCharArray.h>

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
  this->m_SpatialObject = TubeNetType::New();
  this->m_SpatialObject->Initialize();
  this->UpdatePolyDataFromSpatialObject();
  this->BuildDefaultColorMap();
}

//------------------------------------------------------------------------------
vtkMRMLSpatialObjectsNode::TubeNetPointerType
vtkMRMLSpatialObjectsNode::GetSpatialObject()
{
  return this->m_SpatialObject;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::SetSpatialObject(TubeNetPointerType object)
{
  if (this->m_SpatialObject == object)
    {
    return;
    }

  this->m_SpatialObject = object;
  this->UpdatePolyDataFromSpatialObject();
  this->BuildDefaultColorMap();
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
  this->m_CleanPolyData->Update();
  return this->m_CleanPolyData->GetOutput();
}
#else
//------------------------------------------------------------------------------
vtkAlgorithmOutput* vtkMRMLSpatialObjectsNode::
GetFilteredPolyDataConnection( void )
{
  return this->m_CleanPolyData->GetOutputPort();
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
  vtkMRMLSpatialObjectsDisplayNode* displayNode =
    TemplatedAddDisplayNode<vtkMRMLSpatialObjectsLineDisplayNode>(this);
  this->UpdateCleaning();
  return displayNode;
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddTubeDisplayNode( void )
{
  vtkMRMLSpatialObjectsDisplayNode* displayNode =
    TemplatedAddDisplayNode<vtkMRMLSpatialObjectsTubeDisplayNode>(this);
  this->UpdateCleaning();
  return displayNode;
}

//----------------------------------------------------------------------------
vtkMRMLSpatialObjectsDisplayNode* vtkMRMLSpatialObjectsNode::
AddGlyphDisplayNode( void )
{
  vtkMRMLSpatialObjectsDisplayNode* displayNode =
    TemplatedAddDisplayNode<vtkMRMLSpatialObjectsGlyphDisplayNode>(this);
  this->UpdateCleaning();
  return displayNode;
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

  this->m_ShuffledIds->Initialize();
  this->m_ShuffledIds->SetNumberOfTuples(numberOfPairs);
  for(vtkIdType i = 0;  i < numberOfPairs; ++i)
    {
    this->m_ShuffledIds->SetValue(i, idVector[i]);
    }

  this->UpdateCleaning();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::PrepareCleaning( void )
{
  this->m_ShuffledIds = vtkIdTypeArray::New();
  this->m_CleanPolyData = vtkCleanPolyData::New();
  this->m_CleanPolyData->ConvertLinesToPointsOff();
  this->m_CleanPolyData->ConvertPolysToLinesOff();
  this->m_CleanPolyData->ConvertStripsToPolysOff();
  this->m_CleanPolyData->PointMergingOff();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdateCleaning( void )
{
  if(!this->GetPolyData())
    {
    return;
    }

#if VTK_MAJOR_VERSION <= 5
  this->m_CleanPolyData->SetInput(this->GetPolyData());
#else
  this->m_CleanPolyData->SetInputConnection(this->GetPolyDataConnection());
#endif

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
  this->m_CleanPolyData->Delete();
  this->m_ShuffledIds->Delete();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::UpdatePolyDataFromSpatialObject( void )
{
  typedef itk::Point<double, 3>                PointType;
  typedef itk::VesselTubeSpatialObject<3>      VesselTubeType;
  typedef VesselTubeType::TubePointType        VesselTubePointType;
  typedef itk::IndexValueType                  TubeIdType;

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    this->m_SpatialObject->GetChildren(
      this->m_SpatialObject->GetMaximumDepth(), childName);

  // -----------------------------------------------------------------------
  // Copy skeleton points from vessels into polydata structure
  // -----------------------------------------------------------------------

  // Initialize the SpatialObject
  // Count number of points && remove dupplicate
  int totalNumberOfPoints = 0;
  TubeIdType maxTubeId = 0;
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
    if( maxTubeId < currTube->GetId() )
      {
      maxTubeId = currTube->GetId();
      }
    }

  //Making sure tubeId is unique to a tube
  std::set< TubeIdType > tubeIds;
  std::set< TubeIdType >::iterator it;
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

    TubeIdType currTubeId = currTube->GetId();
    it = tubeIds.find( currTubeId );
    if( it ==  tubeIds.end() )
      {
      tubeIds.insert( currTubeId );
      }
    else
      {
      currTube->SetId( maxTubeId + 1 );
      maxTubeId++;
      }
    }
  // Create the points
  vtkNew<vtkPoints> vesselsPoints;
  vesselsPoints->SetDataTypeToDouble();
  vesselsPoints->SetNumberOfPoints(totalNumberOfPoints);

  // Create the Lines
  vtkNew<vtkCellArray> vesselLinesCA;

  // Create scalar array that indicates TubeID.
  vtkNew<vtkDoubleArray> tubeIDs;
  tubeIDs->SetName("TubeIDs");
  tubeIDs->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates Parent TubeID.
  vtkNew<vtkDoubleArray> parentTubeIDs;
  parentTubeIDs->SetName("ParentTubeIDs");
  parentTubeIDs->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates root TubeID.
  vtkNew<vtkDoubleArray> rootTubeIDs;
  rootTubeIDs->SetName("RootTubeIDs");
  rootTubeIDs->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates whether or not a tube is root
  vtkNew<vtkUnsignedCharArray> rootIndicator;
  rootIndicator->SetName("IsRoot");
  rootIndicator->SetNumberOfTuples(totalNumberOfPoints);

  // Create scalar array that indicates native point colors given in tre file
  vtkNew< vtkUnsignedCharArray > tubeColors;
  tubeColors->SetName("TubeColor");
  tubeColors->SetNumberOfTuples(4 * totalNumberOfPoints);
  tubeColors->SetNumberOfComponents( 4 );

  // Create scalar array that indicates native point colors given in tre file
  vtkNew< vtkUnsignedCharArray > tubePointColors;
  tubePointColors->SetName("TubePointColor");
  tubePointColors->SetNumberOfTuples(4 * totalNumberOfPoints);
  tubePointColors->SetNumberOfComponents( 4 );

  // Create scalar array that indicates the radius at each
  // centerline point.
  vtkNew<vtkDoubleArray> tubeRadius;
  tubeRadius->SetName("TubeRadius");
  tubeRadius->SetNumberOfTuples(totalNumberOfPoints);

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

    // compute root tube ID
    VesselTubeType* curRootTube = currTube;

    while( !curRootTube->GetRoot() )
      {
      VesselTubeType* curParentTube =
        dynamic_cast<VesselTubeType*>( curRootTube->GetParent() );
      if( !curParentTube )
        {
        break;
        }
      curRootTube = curParentTube;
      }

    VesselTubeType::PropertyType::PixelType curTubeColor =
      currTube->GetProperty()->GetColor();

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

      // Parent TubeID
      parentTubeIDs->SetTuple1(pointID, currTube->GetParentId());

      // Root TubeID
      rootTubeIDs->SetTuple1(pointID, curRootTube->GetId());

      // Is the current tube a root
      rootIndicator->SetTuple1(pointID, currTube->GetRoot());

      // Native tube color from tre file
      tubeColors->SetTuple4(pointID,
        curTubeColor.GetRed(),
        curTubeColor.GetGreen(),
        curTubeColor.GetBlue(),
        curTubeColor.GetAlpha() );

      // Native tube point color from tre file
      tubePointColors->SetTuple4(pointID,
        tubePoint->GetRed(),
        tubePoint->GetGreen(),
        tubePoint->GetBlue(),
        tubePoint->GetAlpha() );

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

  // Add the Parent TudeID information
  vesselsPD->GetPointData()->AddArray(parentTubeIDs.GetPointer());

  // Add the Root TudeID information
  vesselsPD->GetPointData()->AddArray(rootTubeIDs.GetPointer());

  // Add the info about whether or not a tube is a root
  vesselsPD->GetPointData()->AddArray(rootIndicator.GetPointer());

  // Add tube point color information
  vesselsPD->GetPointData()->AddArray(tubeColors.GetPointer());

  // Add tube point color information
  vesselsPD->GetPointData()->AddArray(tubePointColors.GetPointer());

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

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::BuildDefaultColorMap( void )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::BuildDefaultColorMap");

  typedef itk::VesselTubeSpatialObject<3>      VesselTubeType;

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList =
    this->m_SpatialObject->GetChildren(
      this->m_SpatialObject->GetMaximumDepth(), childName);

  for (TubeNetType::ChildrenListType::iterator tubeIt = tubeList->begin();
      tubeIt != tubeList->end(); ++tubeIt)
    {
    VesselTubeType* currTube =
      dynamic_cast<VesselTubeType*>((*tubeIt).GetPointer());
    if (!currTube || currTube->GetNumberOfPoints() < 1)
      {
      continue;
      }
    std::map< int, std::vector<double> >::iterator it;
    it = this->m_DefaultColorMap.find(currTube->GetId());
    if (it == this->m_DefaultColorMap.end())
      {
      std::vector<double> color;
      color.push_back(currTube->GetProperty()->GetColor().GetRed());
      color.push_back(currTube->GetProperty()->GetColor().GetGreen());
      color.push_back(currTube->GetProperty()->GetColor().GetBlue());
      this->m_DefaultColorMap[currTube->GetId()] = color;
      }
    }
}

//------------------------------------------------------------------------------
bool vtkMRMLSpatialObjectsNode::GetColorFromDefaultColorMap
  ( int TubeId, std::vector<double> &color )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::GetColorFromDefaultColorMap");

  std::map< int, std::vector<double> >::iterator itDefaultColorMap;
  itDefaultColorMap = this->m_DefaultColorMap.find( TubeId );
  if ( itDefaultColorMap != this->m_DefaultColorMap.end() )
    {
    color = this->m_DefaultColorMap.find( TubeId )->second;
    return true;
    }
  return false;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::InsertSelectedTube( int TubeId )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::InsertSelectedTube");

  this->m_SelectedTubeIds.insert( TubeId );
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::ClearSelectedTubes()
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::ClearSelectedTubes");

  this->m_SelectedTubeIds.clear();
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsNode::EraseSelectedTube( int TubeId )
{
  vtkDebugMacro("vtkMRMLSpatialObjectsNode::EraseSelectedTube");
  std::set<int>::iterator it = this->m_SelectedTubeIds.find( TubeId );
  if ( it != this->m_SelectedTubeIds.end() )
    {
    this->m_SelectedTubeIds.erase( it );
    }
}
