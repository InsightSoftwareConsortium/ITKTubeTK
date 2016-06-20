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
#include <vtkVersion.h>
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
  typedef itk::GroupSpatialObject<3>  TubeNetType;
  typedef TubeNetType::Pointer        TubeNetPointerType;

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
  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName( void )
    {return "SpatialObjects";}

  ///
  /// Get the subsampled PolyData converted from the real data in the node.
#if VTK_MAJOR_VERSION <= 5
   virtual vtkPolyData* GetFilteredPolyData( void );
#else
  virtual vtkAlgorithmOutput* GetFilteredPolyDataConnection( void );
#endif

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

  // Description:
  // Reset and recompute the polydata from the spatial object.
  virtual void UpdatePolyDataFromSpatialObject( void );

  //Description
  //Build the colormap for the default color of the tubes of the spatial object
  void BuildDefaultColorMap( void );
  //Get color
  bool GetColorFromDefaultColorMap( int TubeId, std::vector<double> &color );

  //Get SelectedTubeIds
  std::set<int> const & GetSelectedTubeIds() const
    { return m_SelectedTubeIds; }

  //Set SelectedTubeIds
  void InsertSelectedTube( int TubeId );
  //Clear Selected tubeIds
  void ClearSelectedTubes();
  //Erase a tube from selectetubeIds
  void EraseSelectedTube( int TubeId );

protected:
  vtkMRMLSpatialObjectsNode( void );
  ~vtkMRMLSpatialObjectsNode( void );
  vtkMRMLSpatialObjectsNode(const vtkMRMLSpatialObjectsNode&);
  void operator=(const vtkMRMLSpatialObjectsNode&);

  // Description
  // Contains the SpatialObject structure used to generate the differents
  // PolyData for visualization and allow keeping further informations
  // for object processing and editions.
  TubeNetPointerType m_SpatialObject;

  vtkIdTypeArray* m_ShuffledIds;

  virtual void PrepareCleaning( void );
  virtual void UpdateCleaning( void );
  virtual void RemoveCleaning( void );

  vtkCleanPolyData* m_CleanPolyData;

  std::set<int>                           m_SelectedTubeIds;
  std::map< int, std::vector<double> >    m_DefaultColorMap;

private:
  void Reset( vtkMRMLNode* ) {};
}; // End class vtkMRMLSpatialObjectsNode

#endif // End !defined(__vtkMRMLSpatialObjectsNode_h)
