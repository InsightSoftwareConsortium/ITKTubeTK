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

// .NAME vtkSlicerSpatialObjectsLogic -
// slicer logic class for spatial objects manipulation.
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the spatial objects.

#ifndef __vtkSlicerSpatialObjectsLogic_h
#define __vtkSlicerSpatialObjectsLogic_h

#include <vtkSlicerModuleLogic.h>
#include <vtkSlicerSpatialObjectsModuleLogicExport.h>

// STD includes
#include <cstdlib>

class vtkMRMLSpatialObjectsNode;


class VTK_SLICER_SPATIALOBJECTS_MODULE_LOGIC_EXPORT vtkSlicerSpatialObjectsLogic
 : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerSpatialObjectsLogic *New( void );
  vtkTypeMacro(vtkSlicerSpatialObjectsLogic,vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Read and update the given spatial object from the new spatial object
  // stored in the filename. If not storage node exists for the given spatial
  // object, one will be created.
  void SetSpatialObject(
    vtkMRMLSpatialObjectsNode* spatialObject, const char* filename);

  // Description:
  // Create SpatialObjectsNode and
  // read their polydata from a specified directory if a filename is specified.
  // Also create the logic object for its display.
  vtkMRMLSpatialObjectsNode* AddSpatialObject(const char* filename = 0);

  // Description:
  // Create SpatialObjectsNode and
  // read their polydata from a specified directory.
  // Files matching suffix are read
  // Internally calls AddSpatialObjects for each file.
  int AddSpatialObjects(const char* dirname, const char* suffix);

  // Description:
  // Create SpatialObjectsNode and
  // read their polydata from a specified directory.
  // Files matching all suffixes are read
  // Internally calls AddSpatialObjects for each file.
  int AddSpatialObjects(const char* dirname, std::vector<std::string> suffix);

  // Description:
  // Write SpatialObjectsNode's polydata  to a specified file.
  int SaveSpatialObject(const char* filename,
                        vtkMRMLSpatialObjectsNode *spatialObjectsNode);

  // Description:
  // Utility function to copy a spatial object from a node to another
  // (implemented since python doesn't wrap itk objects)
  void CopySpatialObject(
    vtkMRMLSpatialObjectsNode* from, vtkMRMLSpatialObjectsNode* to);

  // Description:
  // Given a filename pointing to a spatial object file and a spatial object
  // node, merge the spatial object read from file in the recipient node.
  // The updated recipient is returned.
  // \sa MergeSpatialObject
  vtkMRMLSpatialObjectsNode* MergeSpatialObjectFromFilename(
    vtkMRMLSpatialObjectsNode* recipient, const char* filename);

  // Description:
  // Given two spatial object nodes, add the donor spatial object to the
  // recipient. The donor stays intact and the updated recipient is returned.
  // \sa MergeSpatialObjectFromFilename
  vtkMRMLSpatialObjectsNode* MergeSpatialObject(
    vtkMRMLSpatialObjectsNode* recipient, vtkMRMLSpatialObjectsNode* donor);

  // Description:
  // Utility function to add the display nodes to a given spatial object node.
  void AddDisplayNodes(vtkMRMLSpatialObjectsNode* spatialObjectsNode);

  // Description:
  // Register MRML Node classes to Scene.
  // Called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes( void );

protected:
  vtkSlicerSpatialObjectsLogic( void );
  ~vtkSlicerSpatialObjectsLogic( void );
  vtkSlicerSpatialObjectsLogic(const vtkSlicerSpatialObjectsLogic&);
  void operator=(const vtkSlicerSpatialObjectsLogic&);

  // Description:
  // Collection of pointers to display logic objects
  // for spatial objects nodes in the scene.
  vtkCollection *DisplayLogicCollection;

}; // End class vtkSlicerSpatialObjectsLogic

#endif // End !defined(__vtkSlicerSpatialObjectsLogic_h)
