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

#ifndef __qSlicerSpatialObjectsWidget_h
#define __qSlicerSpatialObjectsWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

// qMRML includes
#include <qSlicerSpatialObjectsModuleWidgetsExport.h>

class qSlicerSpatialObjectsWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLSpatialObjectsNode;
class vtkMRMLSpatialObjectsDisplayNode;
class vtkMRMLSpatialObjectsDisplayPropertiesNode;

class Q_SLICER_MODULE_SPATIALOBJECTS_WIDGETS_EXPORT qSlicerSpatialObjectsWidget
 : public qSlicerWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:
  typedef qSlicerWidget Superclass;
  qSlicerSpatialObjectsWidget( QWidget *parent=0 );
  virtual ~qSlicerSpatialObjectsWidget();

  QColor color() const;
  double opacity() const;
  double ambient() const;
  double diffuse() const;
  double specular() const;
  double specularPower() const;
  bool backfaceCulling() const;

  vtkMRMLSpatialObjectsNode* SpatialObjectsNode() const;
  vtkMRMLSpatialObjectsDisplayNode* SpatialObjectsDisplayNode() const;
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    SpatialObjectsDisplayPropertiesNode() const;

public slots:
  void setSpatialObjectsNode( vtkMRMLNode *node, int DisplayNodeIndex = 0 );
  void setSpatialObjectsNode( vtkMRMLSpatialObjectsNode* SpatialObjectsNode, int DisplayNodeIndex = 0 );

  void setVisibility( bool );
  void setColorByScalar();
  void onColorByScalarChanged( int );
  void onColorByScalarValuesChanged( double,double );
  void setColorByCellScalars();
  void setColorBySolid();
  void setSliceIntersection( int );
  void setSliceIntersectionThickness( double );
  void clickColorBySolid( bool );
  void onColorBySolidChanged( const QColor& );
  void setColorByCellScalarsColorTable( vtkMRMLNode* );
  void setOpacity( double );

  /// Set the values on the display node
  void setColor( const QColor& );
  void setAmbient( double );
  void setDiffuse( double );
  void setSpecular( double );
  void setSpecularPower( double );
  void setBackfaceCulling( bool );

protected slots:
  void updateWidgetFromMRML();

protected:
  void setSpatialObjectsDisplayNode( vtkMRMLNode *node );
  void setSpatialObjectsDisplayNode( vtkMRMLSpatialObjectsDisplayNode *node );

  QScopedPointer<qSlicerSpatialObjectsWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE( qSlicerSpatialObjectsWidget );
  Q_DISABLE_COPY( qSlicerSpatialObjectsWidget );

  int inProcess;
};

#endif
