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

// MRML includes
#include "vtkMRMLSpatialObjectsStorageNode.h"
#include "vtkMRMLSpatialObjectsNode.h"
#include "vtkMRMLSpatialObjectsDisplayNode.h"

// VTK includes
#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkCleanPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkStringArray.h>

// ITK includes
#include <itkTubeSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>

//------------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSpatialObjectsStorageNode);

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------
bool vtkMRMLSpatialObjectsStorageNode::
CanReadInReferenceNode(vtkMRMLNode *refNode)
{
  return refNode->IsA("vtkMRMLSpatialObjectsNode");
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsStorageNode::ReadDataInternal(vtkMRMLNode *refNode)
{
  vtkMRMLSpatialObjectsNode* spatialObjectsNode =
    vtkMRMLSpatialObjectsNode::SafeDownCast(refNode);

  if(Superclass::ReadDataInternal(refNode) != 0)
    {
    return 0;
    }

  std::string fullName = this->GetFullNameFromFileName();
  if(fullName == std::string(""))
    {
    vtkErrorMacro("ReadData: File name not specified");
    return 0;
    }

  // compute file prefix
  std::string name(fullName);
  std::string::size_type loc = name.find_last_of(".");
  if( loc == std::string::npos )
    {
    vtkErrorMacro("ReadData: no file extension specified: " << name.c_str());
    return 0;
    }
  std::string extension = name.substr(loc);
  vtkDebugMacro("ReadData: extension = " << extension.c_str());

  int result = 1;
  try
  {
    if(extension == std::string(".tre"))
      {
      ReaderType::Pointer reader = ReaderType::New();
      reader->SetFileName(fullName);
      reader->Update();
      // WARNING : Should we check if the tube contains less than 2 points...

      char childName[] = "Tube";
      TubeNetType::ChildrenListType* tubeList =
        reader->GetGroup()->GetChildren(reader->GetGroup()->GetMaximumDepth(), childName);

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
        TubeType* currTube =
          static_cast<TubeType*>((*tubeIT).GetPointer());

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
        TubeType* currTube =
          static_cast<TubeType*>((*tubeIT).GetPointer());

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
        const double* axesRatio = currTube->GetSpacing();

        int index = 0;
        std::vector<TubePointType>::iterator  tubePointIterator;
        for(tubePointIterator = currTube->GetPoints().begin();
             tubePointIterator != currTube->GetPoints().end();
             ++tubePointIterator, ++pointID, ++index)
          {
          PointType inputPoint = tubePointIterator->GetPosition();
          pointIDs[index] = pointID;

          // Insert points using the element spacing information.
          vesselsPoints->SetPoint(pointID,
                                  inputPoint[0],
                                  inputPoint[1] * axesRatio[1] / axesRatio[0],
                                  inputPoint[2] * axesRatio[2] / axesRatio[0]);

          // TubeID
          tubeIDs->SetTuple1(pointID, currTube->GetId());

          // Radius
          tubeRadius->SetTuple1(pointID, tubePointIterator->GetRadius());

          // Tangeantes
          tan1->SetTuple3(pointID,
                          (*tubePointIterator).GetNormal1()[0],
                          (*tubePointIterator).GetNormal1()[1],
                          (*tubePointIterator).GetNormal1()[2]);

          tan2->SetTuple3(pointID,
                          (*tubePointIterator).GetNormal2()[0],
                          (*tubePointIterator).GetNormal2()[1],
                          (*tubePointIterator).GetNormal2()[2]);

          // Medialness & Ridgness
          if(tubePointIterator->GetMedialness() != 0)
            {
            containsMidialnessInfo = true;
            }
          medialness->SetTuple1(pointID, tubePointIterator->GetMedialness());

          if(tubePointIterator->GetRidgeness() != 0)
            {
            containsRidgnessInfo = true;
            }
          ridgeness->SetTuple1(pointID, tubePointIterator->GetRidgeness());
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

      spatialObjectsNode->SetAndObservePolyData(vesselsPD.GetPointer());
      spatialObjectsNode->SetSpatialObject(reader->GetGroup());
      delete tubeList;
    }
  }
  catch(...)
  {
    result = 0;
  }

  if(spatialObjectsNode->GetPolyData() != NULL)
    {
    // is there an active scalar array?
    if(spatialObjectsNode->GetDisplayNode())
      {
      double *scalarRange = spatialObjectsNode->GetPolyData()->GetScalarRange();
      if(scalarRange)
        {
        vtkDebugMacro("ReadData: setting scalar range " << scalarRange[0]
                      << ", " << scalarRange[1]);
        spatialObjectsNode->GetDisplayNode()->SetScalarRange(scalarRange);
        }
      }

    spatialObjectsNode->GetPolyData()->Modified();
    }

  return result;
}

//------------------------------------------------------------------------------
int vtkMRMLSpatialObjectsStorageNode::WriteDataInternal(vtkMRMLNode *refNode)
{
  vtkMRMLSpatialObjectsNode* spatialObjects =
    vtkMRMLSpatialObjectsNode::SafeDownCast(refNode);

  const std::string fullName = this->GetFullNameFromFileName();
  if(fullName == std::string(""))
    {
    vtkErrorMacro("vtkMRMLModelNode: File name not specified");
    return 0;
    }

  const std::string extension =
    itksys::SystemTools::GetFilenameLastExtension(fullName);

  int result = 1;
  if(extension == ".tre")
    {
    try
      {
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(fullName.c_str());
      writer->SetInput(spatialObjects->GetSpatialObject());
      writer->Update();
      }
    catch(...)
      {
      result = 0;
      vtkErrorMacro("Error occured writing Spatial Objects: "
                    << fullName.c_str());
      }
    }
  else
    {
    result = 0;
    vtkErrorMacro( << "No file extension recognized: " << fullName.c_str());
    }

  return result;
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::InitializeSupportedReadFileTypes( void )
{
  this->SupportedReadFileTypes->InsertNextValue("SpatialObject (.tre)");
}

//------------------------------------------------------------------------------
void vtkMRMLSpatialObjectsStorageNode::InitializeSupportedWriteFileTypes( void )
{
  this->SupportedWriteFileTypes->InsertNextValue("SpatialObject (.tre)");
}

//------------------------------------------------------------------------------
const char* vtkMRMLSpatialObjectsStorageNode::GetDefaultWriteFileExtension( void )
{
  return "tre";
}
