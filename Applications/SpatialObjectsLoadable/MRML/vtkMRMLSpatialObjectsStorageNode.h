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

///  vtkMRMLSpatialObjectsStorageNode -
/// MRML node for SpatialObjects storage on disk.
///
/// The storage node has methods to read/write itkSpatialObjects from disk and
/// generates the PolyData.

#ifndef __vtkMRMLSpatialObjectsStorageNode_h
#define __vtkMRMLSpatialObjectsStorageNode_h

// MRML includes
#include "vtkMRMLModelStorageNode.h"

// SpatialObjects includes
#include "vtkSlicerSpatialObjectsModuleMRMLExport.h"

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

  typedef itk::VesselTubeSpatialObjectPoint<3> TubePointType;
  typedef itk::Point<double, 3>                PointType;
  typedef itk::VesselTubeSpatialObject<3>      TubeType;

  static vtkMRMLSpatialObjectsStorageNode *New();
  vtkTypeMacro(vtkMRMLSpatialObjectsStorageNode, vtkMRMLModelStorageNode);
  void PrintSelf(ostream& os, vtkIndent indent);

  virtual vtkMRMLNode* CreateNodeInstance();

  ///
  /// Get node XML tag name (like Storage, Model)
  virtual const char* GetNodeTagName() {return "SpatialObjectsStorage";};

  ///
  /// Return a default file extension for writting
  virtual const char* GetDefaultWriteFileExtension();

  ///
  /// Return true if the node can be read in
  virtual bool CanReadInReferenceNode(vtkMRMLNode *refNode);

protected:
  vtkMRMLSpatialObjectsStorageNode(){};
  ~vtkMRMLSpatialObjectsStorageNode(){};
  vtkMRMLSpatialObjectsStorageNode(const vtkMRMLSpatialObjectsStorageNode&);
  void operator=(const vtkMRMLSpatialObjectsStorageNode&);

  ///
  /// Initialize all the supported read file types
   virtual void InitializeSupportedReadFileTypes();

  ///
  /// Initialize all the supported write file types
  virtual void InitializeSupportedWriteFileTypes();

  /// Read data and set it in the referenced node
  virtual int ReadDataInternal(vtkMRMLNode *refNode);

  /// Write data from a  referenced node
  virtual int WriteDataInternal(vtkMRMLNode *refNode);
};

#endif
