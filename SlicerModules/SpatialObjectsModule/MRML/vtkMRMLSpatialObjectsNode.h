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

/// vtkMRMLSpatialObjectsNode -
/// MRML node to represent SpatialObject structures.
///
/// vtkMRMLSpatialObjects nodes contain trajectories ("tubes")
/// from vessels, internally represented as vtkPolyData.
/// A SpatialObjects node contains many tubes and forms the smallest
/// logical unit of the vessel network.
/// that MRML will manage/read/write
/// each vessel has accompanying spatial object data.
/// Visualization parameters for these nodes are controlled by the
/// vtkMRMLSpatialObjectsDisplayNode class.

#ifndef __vtkMRMLSpatialObjectsNode_h
#define __vtkMRMLSpatialObjectsNode_h

#include <vtkMRMLModelNode.h>

// Spatial Objects includes
#include <vtkSlicerSpatialObjectsModuleMRMLExport.h>
#include <itkGroupSpatialObject.h>

class vtkMRMLSpatialObjectsDisplayNode;
class vtkExtractSelectedPolyDataIds;
class vtkMRMLAnnotationNode;
class vtkIdTypeArray;
class vtkExtractPolyDataGeometry;
class vtkPlanes;
class vtkCleanPolyData;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT vtkMRMLSpatialObjectsNode
  : public vtkMRMLModelNode
{
public:
  typedef itk::GroupSpatialObject<3> TubeNetType;
  typedef TubeNetType::Pointer TubeNetPointerType;

  static vtkMRMLSpatialObjectsNode* New( void );
  vtkTypeMacro(vtkMRMLSpatialObjectsNode, vtkMRMLModelNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //----------------------------------------------------------------------------
  /// MRMLNode methods
  //----------------------------------------------------------------------------
  virtual vtkMRMLNode* CreateNodeInstance( void );

  ///
  /// Read node attributes from XML (MRML) file
  virtual void ReadXMLAttributes(const char** atts);

  ///
  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  ///
  /// Copy the node's attributes to this object
  /// Does NOT copy: ID, FilePrefix, Name, ID
  virtual void Copy(vtkMRMLNode *node);

  ///
  /// Reset the node to an empty spatial object
  virtual void Reset();

  ///
  /// Updates this node if it depends on other nodes
  /// when the node is deleted in the scene
  virtual void UpdateReferences( void );

  ///
  /// Finds the storage node and read the data
  virtual void UpdateScene(vtkMRMLScene *scene);

  ///
  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName( void )
  {return "SpatialObjects";}

  /// Get the subsampling ratio for the polydata
  vtkGetMacro(SubsamplingRatio, float);

  /// Set the subsampling ratio for the polydata
  //
  virtual void SetSubsamplingRatio(float);
  virtual float GetSubsamplingRatioMinValue( void ) {return 0.;}
  virtual float GetSubsamplingRatioMaxValue( void ) {return 1.;}

  ///
  /// Get the subsampled PolyData converted from the real data in the node.
  virtual vtkPolyData* GetFilteredPolyData( void );

  ///
  /// Get associated line display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetLineDisplayNode( void );

  ///
  /// Get associated tube display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetTubeDisplayNode( void );

  ///
  /// Get associated glyph display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetGlyphDisplayNode( void );

  ///
  /// Add line display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddLineDisplayNode( void );

  ///
  /// Add tube display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddTubeDisplayNode( void );

  ///
  /// Add glyph display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddGlyphDisplayNode( void );

  ///
  /// Create and return default storage node or NULL if does not have one.
  virtual vtkMRMLStorageNode* CreateDefaultStorageNode( void );

  ///
  /// Create default display nodes
  virtual void CreateDefaultDisplayNodes( void );

  // Description:
  // Get/Set the SpatialObject when a new node is set
  TubeNetPointerType GetSpatialObject( void );
  void SetSpatialObject(TubeNetPointerType object);

  /// Set and observe poly data for this model
  virtual void SetAndObservePolyData(vtkPolyData* polyData);

protected:
  vtkMRMLSpatialObjectsNode( void );
  ~vtkMRMLSpatialObjectsNode( void );
  vtkMRMLSpatialObjectsNode(const vtkMRMLSpatialObjectsNode&);
  void operator=(const vtkMRMLSpatialObjectsNode&);

  // Description
  // Contains the SpatialObject structure used to generate the differents
  // PolyData for visualization and allow keeping further informations
  // for object processing and editions.
  TubeNetPointerType SpatialObject;

  vtkIdTypeArray* ShuffledIds;

  virtual void PrepareSubsampling( void );
  virtual void UpdateSubsampling( void );
  virtual void CleanSubsampling( void );
  virtual void UpdatePolyDataFromSpatialObject( void );

  vtkCleanPolyData* CleanPolyDataPostSubsampling;
  float SubsamplingRatio;

}; // End class vtkMRMLSpatialObjectsNode

#endif // End !defined(__vtkMRMLSpatialObjectsNode_h)
