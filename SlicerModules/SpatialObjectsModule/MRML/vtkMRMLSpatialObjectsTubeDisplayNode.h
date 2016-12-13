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

/// vtkMRMLSpatialObjectsTubeDisplayNode -
/// MRML node to represent display properties for spatial objects.
///
/// vtkMRMLSpatialObjectsTubeDisplayNode nodes store
/// display properties of vessels including color type, radius of the tube,
/// number of sides, display on/off for glyphs and display of
/// trajectory as a line or tube.

#ifndef __vtkMRMLSpatialObjectsTubeDisplayNode_h
#define __vtkMRMLSpatialObjectsTubeDisplayNode_h

#include "vtkMRMLSpatialObjectsDisplayNode.h"

class vtkAssignAttribute;
class vtkPolyData;
class vtkPolyDataTensorToColor;
class vtkTubeFilter;
class vtkPolyDataColorLinesByOrientation;

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT
vtkMRMLSpatialObjectsTubeDisplayNode : public vtkMRMLSpatialObjectsDisplayNode
{
 public:
  static vtkMRMLSpatialObjectsTubeDisplayNode* New( void );
  vtkTypeMacro( vtkMRMLSpatialObjectsTubeDisplayNode,
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
  {return "SpatialObjectsTubeDisplayNode";}

  ///
  /// Return the polydata that was set by SetInputPolyData
  virtual vtkPolyData* GetInputPolyData( void );

  ///
  /// Update the pipeline based on this node attributes
  virtual void UpdatePolyDataPipeline( void );

  //----------------------------------------------------------------------------
  /// Display Information: Geometry to display ( not mutually exclusive )
  //----------------------------------------------------------------------------

  ///
  /// The minimum tube radius.
  vtkSetMacro( TubeRadius, double );
  vtkGetMacro( TubeRadius, double );

  ///
  /// Number of tube sides.
  vtkSetMacro( TubeNumberOfSides, int );
  vtkGetMacro( TubeNumberOfSides, int );

protected:
  vtkMRMLSpatialObjectsTubeDisplayNode( void );
  ~vtkMRMLSpatialObjectsTubeDisplayNode( void );
  vtkMRMLSpatialObjectsTubeDisplayNode(
    const vtkMRMLSpatialObjectsTubeDisplayNode& );
  void operator=( const vtkMRMLSpatialObjectsTubeDisplayNode& );

  /// To be reimplemented in subclasses if the input of the pipeline changes
#if ( VTK_MAJOR_VERSION <= 5 )
  virtual void SetInputToPolyDataPipeline( vtkPolyData* polyData );
#else
  virtual void SetInputToPolyDataPipeline( vtkAlgorithmOutput* polyDataConnection );
#endif

  /// Return the polydata that is processed by the display node.
  /// This is the polydata that needs to be connected with the mappers.
  virtual vtkAlgorithmOutput* GetOutputPort( void );

  /// Properties
  int    TubeNumberOfSides;
  double TubeRadius;

  /// Pipeline
  vtkAssignAttribute* amontAssignAttribute;
  vtkAssignAttribute* avalAssignAttribute;
  /// Creates with a tube for every point based on the radius.
  vtkTubeFilter*      TubeFilter;

}; // End class vtkMRMLSpatialObjectsTubeDisplayNode

#endif // End !defined( __vtkMRMLSpatialObjectsTubeDisplayNode_h )
