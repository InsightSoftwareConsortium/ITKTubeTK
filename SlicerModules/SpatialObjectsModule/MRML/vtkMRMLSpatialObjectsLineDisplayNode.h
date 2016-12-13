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

/// vtkMRMLSpatialObjectsLineDisplayNode -
/// MRML node to represent display properties for spatial objects.
///
/// vtkMRMLSpatialObjectsTubeDisplayNode nodes store
/// display properties of vessels including color type, radius of the tube,
/// number of sides, display on/off for glyphs and display of
/// trajectory as a line or tube.

#ifndef __vtkMRMLSpatialObjectsLineDisplayNode_h
#define __vtkMRMLSpatialObjectsLineDisplayNode_h

#include "vtkMRMLSpatialObjectsDisplayNode.h"

class vtkPolyData;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT
vtkMRMLSpatialObjectsLineDisplayNode : public vtkMRMLSpatialObjectsDisplayNode
{
public:
  static vtkMRMLSpatialObjectsLineDisplayNode* New( void );
  vtkTypeMacro( vtkMRMLSpatialObjectsLineDisplayNode,
               vtkMRMLSpatialObjectsDisplayNode );
  void PrintSelf( ostream& os, vtkIndent indent );

  //----------------------------------------------------------------------------
  /// MRMLNode methods
  //----------------------------------------------------------------------------
  virtual vtkMRMLNode* CreateNodeInstance( void );

  ///
  /// Read node attributes from XML ( MRML ) file
  virtual void ReadXMLAttributes( const char** atts );

  ///
  /// Write this node's information to a MRML file in XML format.
  virtual void WriteXML( ostream& of, int indent );

  ///
  /// Copy the node's attributes to this object
  /// Does NOT copy: ID, FilePrefix, Name, ID
  virtual void Copy( vtkMRMLNode *node );

  ///
  /// Get node XML tag name ( like Volume, UnstructuredGrid )
  virtual const char* GetNodeTagName( void )
  {return "SpatialObjectsLineDisplayNode";}

  ///
  /// Update the pipeline based on this node attributes
  virtual void UpdatePolyDataPipeline( void );

protected:
  vtkMRMLSpatialObjectsLineDisplayNode( void );
  ~vtkMRMLSpatialObjectsLineDisplayNode( void );
  vtkMRMLSpatialObjectsLineDisplayNode( 
    const vtkMRMLSpatialObjectsLineDisplayNode& );
  void operator= ( const vtkMRMLSpatialObjectsLineDisplayNode& );

}; // End class vtkMRMLSpatialObjectsLineDisplayNode

#endif // End !defined( __vtkMRMLSpatialObjectsLineDisplayNode_h )
