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
/// vtkMRMLSpatialObjectsGlyphDisplayNode -
/// MRML node to represent display properties for vessels.
///
/// vtkMRMLSpatialObjectsTubeDisplayNode nodes store
/// display properties of vessels including color type, radius of the tube,
/// number of sides, display on/off for glyphs and display of
/// trajectory as a line or tube.

#ifndef __vtkMRMLSpatialObjectsGlyphDisplayNode_h
#define __vtkMRMLSpatialObjectsGlyphDisplayNode_h

#include "vtkMRMLSpatialObjectsDisplayNode.h"

class vtkGlyph3DMapper;
class vtkPolyData;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT
vtkMRMLSpatialObjectsGlyphDisplayNode : public vtkMRMLSpatialObjectsDisplayNode
{
public:
  static vtkMRMLSpatialObjectsGlyphDisplayNode* New();
  vtkTypeMacro(vtkMRMLSpatialObjectsGlyphDisplayNode,
               vtkMRMLSpatialObjectsDisplayNode);
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
  /// Get node XML tag name (like Volume, UnstructuredGrid)
  virtual const char* GetNodeTagName()
  {return "SpatialObjectsGlyphDisplayNode";}

  ///
  /// Update the pipeline based on this node attributes
  virtual void UpdatePolyDataPipeline();

 protected:
  vtkMRMLSpatialObjectsGlyphDisplayNode();
  ~vtkMRMLSpatialObjectsGlyphDisplayNode();
  vtkMRMLSpatialObjectsGlyphDisplayNode(
    const vtkMRMLSpatialObjectsGlyphDisplayNode&);
  void operator=(const vtkMRMLSpatialObjectsGlyphDisplayNode&);

  /// Return the polydata that is processed by the display node.
  /// This is the polydata that needs to be connected with the mappers.
  virtual vtkAlgorithmOutput* GetOutputPort();

  /// Pipeline
  vtkGlyph3DMapper* Glyph3DMapper;
};

#endif
