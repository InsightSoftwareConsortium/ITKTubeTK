/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

/// vtkMRMLSpatialObjectsDisplayPropertiesNode -
/// MRML node for display of a spatial objects.
///
/// This node describes display properties at the single-skeleton point level.
/// A spatial object can be displayed using various
/// scalar invariants and glyphs.
/// This class is used by classes ( e.g. vtkMRMLSpatialObjectsDisplayNode )
/// that handle higher-level display concepts for many spatial objects,
/// such as choosing between scalars/glyphs/etc.
/// For specific display needs.
/// This class inherits from the vtkMRMLColorNode->vtkMRMLColorTableNode
/// superclasses, used for vtkMRMLModelNodes and vtkMRMLVolumeNodes, in order
/// toprovide specific lookup tables for the scalar invariant display.

#ifndef __vtkMRMLSpatialObjectsDisplayPropertiesNode_h
#define __vtkMRMLSpatialObjectsDisplayPropertiesNode_h

#include "vtkMRMLSpatialObjectsDisplayNode.h"
#include <vtkMRMLColorTableNode.h>

//
// Set built-in type. Creates member Set"name"() ( e.g., SetVisibility() );
//
#define SpatialObjectsPropertySetMacro( name,type ) \
virtual void Set##name ( type _arg ) \
  { \
  vtkDebugMacro( << this->GetClassName() << " ( " << this << " ): setting " #name " to " << _arg ); \
  if( this->name != _arg ) \
    { \
    this->name = _arg; \
    if( this->GlyphGeometry == this->Lines || \
        this->GlyphGeometry == this->Tubes ) \
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
  static vtkMRMLSpatialObjectsDisplayPropertiesNode* New( void );
  vtkTypeMacro( vtkMRMLSpatialObjectsDisplayPropertiesNode,
               vtkMRMLColorTableNode );
  void PrintSelf( ostream& os, vtkIndent indent );

  //----------------------------------------------------------------------------
  /// MRMLNode methods
  //----------------------------------------------------------------------------
  virtual vtkMRMLNode* CreateNodeInstance( void );

  ///
  /// Read node attributes from a MRML file in XML format.
  virtual void ReadXMLAttributes( const char** atts );

  ///
  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML( ostream& of, int indent );

  ///
  /// Copy the node's attributes to this object.
  /// Does NOT copy: ID, FilePrefix, Name, ID
  virtual void Copy( vtkMRMLNode *node );

  ///
  /// Get node XML tag name ( like Volume, Model )
  virtual const char* GetNodeTagName( void )
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

  static bool ScalarInvariantHasKnownScalarRange( int ScalarInvariant );
  static void ScalarInvariantKnownScalarRange( int ScalarInvariant,
                                              double range[2] );

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose scalar invariant
  //----------------------------------------------------------------------------

  ///
  /// Get type of scalar invariant ( tensor-derived scalar, invariant to tensor
  /// rotation ) selected for display.
  vtkGetMacro( ScalarInvariant, int );

  ///
  /// Get type of scalar invariant ( tensor-derived scalar, invariant to tensor
  /// rotation ) selected for display.
  vtkSetMacro( ScalarInvariant, int );

  // TODO add all scalar

  ///
  /// Set scalar invariant to relative anisotropy
  void SetScalarInvariantToRelativeAnisotropy( void )
  {this->SetScalarInvariant( this->RelativeAnisotropy );}

  ///
  /// Return a text string describing the ScalarInvariant variable
  virtual const char * GetScalarInvariantAsString( void );

  //----------------------------------------------------------------------------
  /// Display Information: Types of glyph geometry that can be displayed
  //----------------------------------------------------------------------------
  enum { Lines = 0, Tubes = 1, Cones = 2, Disks = 3 };

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose the type of glyph geometry
  //----------------------------------------------------------------------------
  ///
  /// Get the type of glyph geometry ( line, tubes, cones etc. )
  vtkGetMacro( GlyphGeometry, int );

  ///
  /// Set the type of glyph geometry ( line, ellipsoid, etc. )
  /// Update the glyph polydata source
  void SetGlyphGeometry( int geometry );

  void SetGlyphGeometryToLines( void )
  {this->SetGlyphGeometry( this->Lines );}

  void SetGlyphGeometryToTubes( void )
  {this->SetGlyphGeometry( this->Tubes );}

  void SetGlyphGeometryToCones( void )
  {this->SetGlyphGeometry( this->Cones );}

  void SetGlyphGeometryToDisks( void )
  {this->SetGlyphGeometry( this->Disks );}

  ///
  /// Return the lowest and highest integers, for use in looping
  int GetFirstGlyphGeometry( void ) {return this->Lines;}
  int GetLastGlyphGeometry( void ) {return this->Disks;}

  ///
  /// Return a text string describing the GlyphGeometry variable
  virtual const char * GetGlyphGeometryAsString( void );
  virtual const char * GetGlyphGeometryAsString( int );

  //----------------------------------------------------------------------------
  /// Display Information: Parameters of the different geometries
  //----------------------------------------------------------------------------

  ///
  /// Get/Set the scale factor applied to the glyphs.
  vtkGetMacro( GlyphScaleFactor, double );
  vtkSetMacro( GlyphScaleFactor, double );

  ///
  /// Get/Set the resolution of lines displayed
  vtkGetMacro( LineGlyphResolution, int );
  SpatialObjectsPropertySetMacro( LineGlyphResolution, int );

  ///
  /// Get/Set the radius of the tube
  vtkGetMacro( TubeGlyphRadius, double );
  SpatialObjectsPropertySetMacro( TubeGlyphRadius, double );

  ///
  /// Get/Set Number of sides of tube glyph ( 3 gives a triangular tube, etc. )
  vtkGetMacro( TubeGlyphNumberOfSides, int );
  SpatialObjectsPropertySetMacro( TubeGlyphNumberOfSides, int );

  // TODO other representation properties

  //----------------------------------------------------------------------------
  /// Display Information: Functions to choose the type of glyph coloring
  //----------------------------------------------------------------------------

  ///
  /// Get/Set type of scalar invariant selected for display.
  vtkGetMacro( ColorGlyphBy, int );
  vtkSetMacro( ColorGlyphBy, int );

  ///
  /// Return the lowest and highest integers, for use in looping
  static int GetFirstColorGlyphBy( void );
  static int GetLastColorGlyphBy( void );

  ///
  /// Return a text string describing the ColorGlyphBy
  virtual const char* GetColorGlyphByAsString( void );

  ///
  /// Set scalar invariant to LinearMeasure.
  void ColorGlyphByLinearMeasure( void )
  {this->SetColorGlyphBy( this->LinearMeasure );}

  ///
  /// Set scalar invariant to ColorOrientation.
  void ColorGlyphByColorOrientation( void )
  {this->SetColorGlyphBy( this->ColorOrientation );}

  ///
  /// Set scalar invariant to ColorMode.
  void ColorGlyphByColorMode( void )
  {this->SetColorGlyphBy( this->ColorMode );}

  ///
  /// Set scalar invariant to RelativeAnisotropy.
  void ColorGlyphByRelativeAnisotropy( void )
  {this->SetColorGlyphBy( this->RelativeAnisotropy );}

  //--------------------------------------------------------------------------
  /// Convenient functions to get an appropriate glyph source
  //--------------------------------------------------------------------------

  ///
  /// Get a polydata object according to current glyph display settings
  /// ( so a line, sphere, or tube ) to use as a source for a glyphing filter.
  vtkGetObjectMacro( GlyphSource, vtkPolyData );

  ///
  /// Return a text string describing the GlyphScalar variable
  static const char* GetScalarEnumAsString( int val );

  /// Return the lowest and highest integers, for use in looping
  static int GetFirstScalarInvariant( void );
  static int GetLastScalarInvariant( void );

protected:
  vtkMRMLSpatialObjectsDisplayPropertiesNode( void );
  ~vtkMRMLSpatialObjectsDisplayPropertiesNode( void );
  vtkMRMLSpatialObjectsDisplayPropertiesNode( 
    const vtkMRMLSpatialObjectsDisplayPropertiesNode& );
  void operator=( const vtkMRMLSpatialObjectsDisplayPropertiesNode& );

  virtual void SetGlyphSource( vtkPolyData* glyphSource );
  virtual void UpdateGlyphSource( void );

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

}; // End class vtkMRMLSpatialObjectsDisplayPropertiesNode

#endif // End !defined( __vtkMRMLSpatialObjectsDisplayPropertiesNode_h )
