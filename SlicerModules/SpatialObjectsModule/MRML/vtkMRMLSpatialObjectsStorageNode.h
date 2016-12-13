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

///  vtkMRMLSpatialObjectsStorageNode -
/// MRML node for SpatialObjects storage on disk.
///
/// The storage node has methods to read/write itkSpatialObjects from disk and
/// generates the PolyData.

#ifndef __vtkMRMLSpatialObjectsStorageNode_h
#define __vtkMRMLSpatialObjectsStorageNode_h

// MRML includes
#include <vtkMRMLModelStorageNode.h>

// SpatialObjects includes
#include <vtkSlicerSpatialObjectsModuleMRMLExport.h>

#include <itkVesselTubeSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkSpatialObjectWriter.h>
#include <itkGroupSpatialObject.h>
#include <itkPoint.h>

class VTK_SLICER_SPATIALOBJECTS_MODULE_MRML_EXPORT
vtkMRMLSpatialObjectsStorageNode : public vtkMRMLModelStorageNode
{
public:
  typedef itk::SpatialObjectReader<3> ReaderType;
  typedef itk::SpatialObjectWriter<3> WriterType;
  typedef itk::GroupSpatialObject<3>  TubeNetType;

  static vtkMRMLSpatialObjectsStorageNode *New( void );
  vtkTypeMacro( vtkMRMLSpatialObjectsStorageNode, vtkMRMLModelStorageNode );
  void PrintSelf( ostream& os, vtkIndent indent );

  virtual vtkMRMLNode* CreateNodeInstance( void );

  ///
  /// Get node XML tag name ( like Storage, Model )
  virtual const char* GetNodeTagName( void ) {return "SpatialObjectsStorage";}

  ///
  /// Return a default file extension for writting
  virtual const char* GetDefaultWriteFileExtension( void );

  ///
  /// Return true if the node can be read in
  virtual bool CanReadInReferenceNode( vtkMRMLNode *refNode );

protected:
  vtkMRMLSpatialObjectsStorageNode( void ) {}
  ~vtkMRMLSpatialObjectsStorageNode( void ) {}
  vtkMRMLSpatialObjectsStorageNode( const vtkMRMLSpatialObjectsStorageNode& );
  void operator=( const vtkMRMLSpatialObjectsStorageNode& );

  ///
  /// Initialize all the supported read file types
   virtual void InitializeSupportedReadFileTypes( void );

  ///
  /// Initialize all the supported write file types
  virtual void InitializeSupportedWriteFileTypes( void );

  /// Read data and set it in the referenced node
  virtual int ReadDataInternal( vtkMRMLNode *refNode );

  /// Write data from a  referenced node
  virtual int WriteDataInternal( vtkMRMLNode *refNode );

}; // End class vtkMRMLSpatialObjectsStorageNode

#endif // End !defined( __vtkMRMLSpatialObjectsStorageNode_h )
