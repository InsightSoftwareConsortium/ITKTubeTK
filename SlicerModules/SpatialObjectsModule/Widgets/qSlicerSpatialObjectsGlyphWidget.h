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

#ifndef __qSlicerSpatialObjectsGlyphWidget_h
#define __qSlicerSpatialObjectsGlyphWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

// qMRML includes
#include <qSlicerSpatialObjectsModuleWidgetsExport.h>

class qSlicerSpatialObjectsGlyphWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLSpatialObjectsDisplayNode;
class vtkMRMLSpatialObjectsDisplayPropertiesNode;

class Q_SLICER_MODULE_SPATIALOBJECTS_WIDGETS_EXPORT
qSlicerSpatialObjectsGlyphWidget : public qSlicerWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:
  typedef qSlicerWidget Superclass;
  qSlicerSpatialObjectsGlyphWidget( QWidget *parent=0 );
  virtual ~qSlicerSpatialObjectsGlyphWidget();

  vtkMRMLSpatialObjectsDisplayNode* spatialObjectsDisplayNode() const;
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    spatialObjectsDisplayPropertiesNode() const;

public slots:
  void setSpatialObjectsDisplayNode( vtkMRMLNode* );
  void setSpatialObjectsDisplayNode( vtkMRMLSpatialObjectsDisplayNode* );
  void setSpatialObjectsDisplayPropertiesNode( vtkMRMLNode* );
  void setSpatialObjectsDisplayPropertiesNode
    ( vtkMRMLSpatialObjectsDisplayPropertiesNode* );

  void setGlyphScaleFactor( double );
  void setGlyphSpacing( double );
  void setGlyphType( int );
  void setTubeGlyphNumberOfSides( double );
  void setTubeGlyphRadius( double );

protected slots:
  void updateWidgetFromMRMLDisplayNode();
  void updateWidgetFromMRMLDisplayPropertiesNode();

protected:
  QScopedPointer<qSlicerSpatialObjectsGlyphWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE( qSlicerSpatialObjectsGlyphWidget );
  Q_DISABLE_COPY( qSlicerSpatialObjectsGlyphWidget );

  int updating;
};

#endif
