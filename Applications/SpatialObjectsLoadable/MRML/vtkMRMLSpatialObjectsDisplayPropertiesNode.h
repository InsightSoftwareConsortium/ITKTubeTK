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

/// vtkMRMLSpatialObjectsDisplayPropertiesNode -
/// MRML node for display of a spatial objects.
/// 
/// This node describes display properties at the single-skeleton point level.
/// A spatial object can be displayed using various
/// scalar invariants and glyphs.
/// This class is used by classes (e.g. vtkMRMLSpatialObjectsDisplayNode)
/// that handle higher-level display concepts for many spatial objects,
/// such as choosing between scalars/glyphs/etc.
/// For specific display needs.
/// This class inherits from the vtkMRMLColorNode->vtkMRMLColorTableNode
/// superclasses, used for vtkMRMLModelNodes and vtkMRMLVolumeNodes, in order
/// toprovide specific lookup tables for the scalar invariant display.

#ifndef __vtkMRMLSpatialObjectsDisplayPropertiesNode_h
#define __vtkMRMLSpatialObjectsDisplayPropertiesNode_h

#include "vtkMRMLSpatialObjectsDisplayNode.h"
#include "vtkMRMLColorTableNode.h"

//
// Set built-in type. Creates member Set"name"() (e.g., SetVisibility());
//
#define SpatialObjectsPropertySetMacro(name,type) \
virtual void Set##name (type _arg) \
  { \
  vtkDebugMacro(<< this->GetClassName() << " (" << this << "): setting " #name " to " << _arg); \
  if (this->name != _arg) \
    { \
    this->name = _arg; \
    if (this->GlyphGeometry == this->Lines || \
        this->GlyphGeometry == this->Tubes) \
      { \
      this->UpdateGlyphSource(); \
      } \
    this->Modified(); \
    } \
  }

class vtkPolyData;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT
vtkMRMLSpatialObjectsDisplayPropertiesNode : public vtkMRMLColorTableNode
{
public:
  static vtkMRMLSpatialObjectsDisplayPropertiesNode* New();
  vtkTypeMacro(vtkMRMLSpatialObjectsDisplayPropertiesNode,
               vtkMRMLColorTableNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  //----------------------------------------------------------------------------
  /// MRMLNode methods
  //----------------------------------------------------------------------------
  virtual vtkMRMLNode* CreateNodeInstance();

  ///
  /// Read node attributes from a MRML file in XML format.
  virtual void ReadXMLAttributes(const char** atts);

  ///
  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML(ostream& of, int indent);

  ///
  /// Copy the node's attributes to this object.
  /// Does NOT copy: ID, FilePrefix, Name, ID
  virtual void Copy(vtkMRMLNode *node);

  ///
  /// Get node XML tag name (like Volume, Model)
  virtual const char* GetNodeTagName()
  {return "SpatialObjectsDisplayProperties";}

  //----------------------------------------------------------------------------
  /// Display Information:
  /// Types of scalars that may be generated from spatial objects.
  //----------------------------------------------------------------------------
  enum
  {
    LinearMeasure = 0,
    ColorOrientation = 1,
    ColorMode = 2,
    RelativeAnisotropy = 3
  };

  static bool ScalarInvariantHasKnownScalarRange(int ScalarInvariant);
  static void ScalarInvariantKnownScalarRange(int ScalarInvariant,
                                              double range[2]);

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose scalar invariant
  //----------------------------------------------------------------------------

  ///
  /// Get type of scalar invariant (tensor-derived scalar, invariant to tensor 
  /// rotation) selected for display.
  vtkGetMacro(ScalarInvariant, int);

  ///
  /// Get type of scalar invariant (tensor-derived scalar, invariant to tensor 
  /// rotation) selected for display.
  vtkSetMacro(ScalarInvariant, int);
 
  // TODO add all scalar

  ///
  /// Set scalar invariant to relative anisotropy
  void SetScalarInvariantToRelativeAnisotropy()
  {this->SetScalarInvariant(this->RelativeAnisotropy);}

  ///
  /// Return a text string describing the ScalarInvariant variable
  virtual const char * GetScalarInvariantAsString();

  //----------------------------------------------------------------------------
  /// Display Information: Types of glyph geometry that can be displayed
  //----------------------------------------------------------------------------
  enum
  {
    Lines = 0,
    Tubes = 1,
    Cones = 2,
    Disks = 3
  };

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose the type of glyph geometry
  //----------------------------------------------------------------------------
  ///
  /// Get the type of glyph geometry (line, tubes, cones etc.)
  vtkGetMacro(GlyphGeometry, int);

  ///
  /// Set the type of glyph geometry (line, ellipsoid, etc.)
  /// Update the glyph polydata source
  void SetGlyphGeometry(int geometry);

  void SetGlyphGeometryToLines()
  {this->SetGlyphGeometry(this->Lines);}

  void SetGlyphGeometryToTubes()
  {this->SetGlyphGeometry(this->Tubes);}

  void SetGlyphGeometryToCones()
  {this->SetGlyphGeometry(this->Cones);}

  void SetGlyphGeometryToDisks()
  {this->SetGlyphGeometry(this->Disks);}

  ///
  /// Return the lowest and highest integers, for use in looping
  int GetFirstGlyphGeometry() {return this->Lines;}
  int GetLastGlyphGeometry() {return this->Disks;}

  ///
  /// Return a text string describing the GlyphGeometry variable
  virtual const char * GetGlyphGeometryAsString();
  virtual const char * GetGlyphGeometryAsString(int);
  
  //----------------------------------------------------------------------------
  /// Display Information: Parameters of the different geometries
  //----------------------------------------------------------------------------

  ///
  /// Get/Set the scale factor applied to the glyphs.
  vtkGetMacro(GlyphScaleFactor, double);
  vtkSetMacro(GlyphScaleFactor, double);

  ///
  /// Get/Set the resolution of lines displayed
  vtkGetMacro(LineGlyphResolution, int);
  SpatialObjectsPropertySetMacro(LineGlyphResolution, int);

  ///
  /// Get/Set the radius of the tube
  vtkGetMacro(TubeGlyphRadius, double);
  SpatialObjectsPropertySetMacro(TubeGlyphRadius, double);

  ///
  /// Get/Set Number of sides of tube glyph (3 gives a triangular tube, etc.)
  vtkGetMacro(TubeGlyphNumberOfSides, int);
  SpatialObjectsPropertySetMacro(TubeGlyphNumberOfSides, int);

  // TODO other representation properties

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose the type of glyph coloring
  //----------------------------------------------------------------------------

  ///
  /// Get/Set type of scalar invariant selected for display.
  vtkGetMacro(ColorGlyphBy, int);
  vtkSetMacro(ColorGlyphBy, int);

  ///
  /// Return the lowest and highest integers, for use in looping
  static int GetFirstColorGlyphBy();
  static int GetLastColorGlyphBy();
  
  ///
  /// Return a text string describing the ColorGlyphBy
  virtual const char* GetColorGlyphByAsString();
 
  ///
  /// Set scalar invariant to LinearMeasure.
  void ColorGlyphByLinearMeasure()
  {this->SetColorGlyphBy(this->LinearMeasure);}

  ///
  /// Set scalar invariant to ColorOrientation.
  void ColorGlyphByColorOrientation()
  {this->SetColorGlyphBy(this->ColorOrientation);}

  ///
  /// Set scalar invariant to ColorMode.
  void ColorGlyphByColorMode()
  {this->SetColorGlyphBy(this->ColorMode);}

  ///
  /// Set scalar invariant to RelativeAnisotropy.
  void ColorGlyphByRelativeAnisotropy()
  {this->SetColorGlyphBy(this->RelativeAnisotropy);}

  //--------------------------------------------------------------------------
  /// Convenient functions to get an appropriate glyph source
  //--------------------------------------------------------------------------

  ///
  /// Get a polydata object according to current glyph display settings
  /// (so a line, sphere, or tube) to use as a source for a glyphing filter.
  vtkGetObjectMacro(GlyphSource, vtkPolyData);

  ///
  /// Return a text string describing the GlyphScalar variable
  static const char* GetScalarEnumAsString(int val);

  /// Return the lowest and highest integers, for use in looping
  static int GetFirstScalarInvariant();
  static int GetLastScalarInvariant();

protected:
  vtkMRMLSpatialObjectsDisplayPropertiesNode();
  ~vtkMRMLSpatialObjectsDisplayPropertiesNode();
  vtkMRMLSpatialObjectsDisplayPropertiesNode(
    const vtkMRMLSpatialObjectsDisplayPropertiesNode&);
  void operator=(const vtkMRMLSpatialObjectsDisplayPropertiesNode&);

  virtual void SetGlyphSource(vtkPolyData* glyphSource);
  virtual void UpdateGlyphSource();

  /// ---- Parameters that should be written to MRML --- //
  /// Scalar display parameters
  int ScalarInvariant;

  /// Glyph general parameters
  int GlyphGeometry;
  int ColorGlyphBy;
  double GlyphScaleFactor;

  /// Line Glyph parameters
  int LineGlyphResolution;

  /// Tube Glyph parameters
  double TubeGlyphRadius;
  int TubeGlyphNumberOfSides;
  /// ---- End of parameters that should be written to MRML --- //
 
  /// Pipeline
  vtkPolyData* GlyphSource;
};

#endif
