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

#include "vtkSlicerModuleLogic.h"
#include "vtkSlicerSpatialObjectsModuleLogicExport.h"

// STD includes
#include <cstdlib>

class vtkMRMLSpatialObjectsNode;


class VTK_SLICER_SPATIALOBJECTS_MODULE_LOGIC_EXPORT vtkSlicerSpatialObjectsLogic
 : public vtkSlicerModuleLogic
{
public:
  static vtkSlicerSpatialObjectsLogic *New();
  vtkTypeRevisionMacro(vtkSlicerSpatialObjectsLogic,vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Create SpatialObjectsNode and
  // read their polydata from a specified directory.
  // Also create the logic object for its display.
  vtkMRMLSpatialObjectsNode* AddSpatialObject(const char* filename);

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
  // Register MRML Node classes to Scene.
  // Called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();

protected:
  vtkSlicerSpatialObjectsLogic();
  ~vtkSlicerSpatialObjectsLogic();
  vtkSlicerSpatialObjectsLogic(const vtkSlicerSpatialObjectsLogic&);
  void operator=(const vtkSlicerSpatialObjectsLogic&);

  // Description:
  // Collection of pointers to display logic objects
  // for spatial objects nodes in the scene.
  vtkCollection *DisplayLogicCollection;
};

#endif
