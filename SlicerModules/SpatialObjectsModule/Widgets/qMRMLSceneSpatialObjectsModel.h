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

#ifndef __qMRMLSceneSpatialObjectsModel_h
#define __qMRMLSceneSpatialObjectsModel_h

#include <qSlicerSpatialObjectsModuleWidgetsExport.h>

// MRMLWidgets includes
#include <qMRMLSceneDisplayableModel.h>
class qMRMLSceneSpatialObjectsModelPrivate;

/// The Visibility icon is in the same column than the name by default.
class Q_SLICER_MODULE_SPATIALOBJECTS_WIDGETS_EXPORT
qMRMLSceneSpatialObjectsModel : public qMRMLSceneDisplayableModel
{
Q_OBJECT

  /// Control in which column vtkMRMLModelDisplayNode::Color are displayed
  /// ( Qt::DecorationRole ). Even if a vtkMRMLModelNode doesn't have a color
  /// proper, the color of its display node is used. If the model node has
  /// more than one display node and their colors are different, it uses
  /// an invalid color.
  /// A value of -1 ( default ) hides the column
  Q_PROPERTY ( int colorColumn READ colorColumn WRITE setColorColumn )

  // This property holds the column ID where the node tube visisbility is shown.
  /// A value of -1 ( default ) hides the column.
  Q_PROPERTY ( int lineVisibilityColumn READ lineVisibilityColumn
              WRITE setLineVisibilityColumn )

  // This property holds the column ID where the node tube visisbility is shown.
  /// A value of -1 ( default ) hides the column.
  Q_PROPERTY ( int tubeVisibilityColumn READ tubeVisibilityColumn
              WRITE setTubeVisibilityColumn )

  // Property holds the column ID where the node glyph visisbility is shown.
  /// A value of -1 ( default ) hides the column.
  Q_PROPERTY ( int glyphVisibilityColumn READ glyphVisibilityColumn
              WRITE setGlyphVisibilityColumn )

public:
  typedef qMRMLSceneDisplayableModel Superclass;
  qMRMLSceneSpatialObjectsModel( QObject *parent=0 );
  virtual ~qMRMLSceneSpatialObjectsModel();

  int colorColumn()const;
  void setColorColumn( int column );

  int  lineVisibilityColumn()const;
  void setLineVisibilityColumn( int column );

  int tubeVisibilityColumn()const;
  void setTubeVisibilityColumn( int column );

  int glyphVisibilityColumn()const;
  void setGlyphVisibilityColumn( int column );

protected:
  qMRMLSceneSpatialObjectsModel( qMRMLSceneSpatialObjectsModelPrivate* pimpl,
                                QObject *parent=0 );

  /// Reimplemented to listen to the displayable DisplayModifiedEvent event for
  /// visibility check state changes.
  virtual QFlags<Qt::ItemFlag> nodeFlags( vtkMRMLNode* node, int column ) const;
  virtual void updateItemDataFromNode( QStandardItem* item,
                                      vtkMRMLNode* node,
                                      int column );
  virtual void updateNodeFromItemData( vtkMRMLNode* node, QStandardItem* item );
  virtual int maxColumnId()const;

  void updateVilibilityFromNode( QStandardItem* item,
                                vtkMRMLNode* node,
                                bool slice = false );
  void updateVilibilityFromItem( QStandardItem* item,
                                vtkMRMLNode* node,
                                bool slice = false );

private:
  Q_DECLARE_PRIVATE( qMRMLSceneSpatialObjectsModel );
  Q_DISABLE_COPY( qMRMLSceneSpatialObjectsModel );
};

#endif
