/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Michael Jeulin-L, Kitware Inc.

==============================================================================*/

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

#include "vtkMRMLModelNode.h"

// Spatial Objects includes
#include "vtkSlicerSpatialObjectsModuleMRMLExport.h"
#include <itkGroupSpatialObject.h>

class vtkMRMLSpatialObjectsDisplayNode;
class vtkExtractSelectedPolyDataIds;
class vtkMRMLAnnotationNode;
class vtkIdTypeArray;
class vtkExtractPolyDataGeometry;
class vtkPlanes;
class vtkCleanPolyData;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT vtkMRMLSpatialObjectsNode :
  public vtkMRMLModelNode
{
public:
  typedef itk::GroupSpatialObject<3> TubeNetType;

  static vtkMRMLSpatialObjectsNode* New();
  vtkTypeMacro(vtkMRMLSpatialObjectsNode, vtkMRMLModelNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //----------------------------------------------------------------------------
  /// MRMLNode methods
  //----------------------------------------------------------------------------
  virtual vtkMRMLNode* CreateNodeInstance();

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
  /// Updates this node if it depends on other nodes
  /// when the node is deleted in the scene
  virtual void UpdateReferences();

  ///
  /// Finds the storage node and read the data
  virtual void UpdateScene(vtkMRMLScene *scene);

  ///
  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName()
  {return "SpatialObjects";}

  /// Get the subsampling ratio for the polydata
  vtkGetMacro(SubsamplingRatio, float);

  /// Set the subsampling ratio for the polydata
  //
  virtual void SetSubsamplingRatio(float);
  virtual float GetSubsamplingRatioMinValue(){return 0.;}
  virtual float GetSubsamplingRatioMaxValue(){return 1.;}

  ///
  /// Get the subsampled PolyData converted from the real data in the node.
  virtual vtkPolyData* GetFilteredPolyData();

  ///
  /// Get associated line display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetLineDisplayNode();

  ///
  /// Get associated tube display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetTubeDisplayNode();

  ///
  /// Get associated glyph display node or NULL if not set.
  vtkMRMLSpatialObjectsDisplayNode* GetGlyphDisplayNode();

  ///
  /// Add line display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddLineDisplayNode();

  ///
  /// Add tube display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddTubeDisplayNode();

  ///
  /// Add glyph display node if not already present and return it.
  vtkMRMLSpatialObjectsDisplayNode* AddGlyphDisplayNode();

  ///
  /// Create and return default storage node or NULL if does not have one.
  virtual vtkMRMLStorageNode* CreateDefaultStorageNode();

  ///
  /// Create default display nodes
  virtual void CreateDefaultDisplayNodes();

  // Description:
  // Get/Set the SpatialObject when a new node is set
  vtkGetMacro(SpatialObject, TubeNetType*);
  vtkSetMacro(SpatialObject, TubeNetType*);

  /// Set and observe poly data for this model
  virtual void SetAndObservePolyData(vtkPolyData* polyData);

protected:
  vtkMRMLSpatialObjectsNode();
  ~vtkMRMLSpatialObjectsNode();
  vtkMRMLSpatialObjectsNode(const vtkMRMLSpatialObjectsNode&);
  void operator=(const vtkMRMLSpatialObjectsNode&);

  // Description
  // Contains the SpatialObject structure used to generate the differents
  // PolyData for visualization and allow keeping further informations
  // for object processing and editions.
  TubeNetType* SpatialObject;

  vtkIdTypeArray* ShuffledIds;

  virtual void PrepareSubsampling();
  virtual void UpdateSubsampling();
  virtual void CleanSubsampling();

  vtkCleanPolyData* CleanPolyDataPostSubsampling;
  float SubsamplingRatio;
};

#endif
