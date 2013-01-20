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

#ifndef __qSlicerSpatialObjectsGlyphWidget_h
#define __qSlicerSpatialObjectsGlyphWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

// qMRML includes
#include "qSlicerSpatialObjectsModuleWidgetsExport.h"

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
  qSlicerSpatialObjectsGlyphWidget(QWidget *parent=0);
  virtual ~qSlicerSpatialObjectsGlyphWidget();

  vtkMRMLSpatialObjectsDisplayNode* spatialObjectsDisplayNode() const;
  vtkMRMLSpatialObjectsDisplayPropertiesNode*
    spatialObjectsDisplayPropertiesNode() const;

public slots:
  void setSpatialObjectsDisplayNode(vtkMRMLNode*);
  void setSpatialObjectsDisplayNode(vtkMRMLSpatialObjectsDisplayNode*);
  void setSpatialObjectsDisplayPropertiesNode(vtkMRMLNode*);
  void setSpatialObjectsDisplayPropertiesNode
    (vtkMRMLSpatialObjectsDisplayPropertiesNode*);

  void setGlyphScaleFactor(double);
  void setGlyphSpacing(double);
  void setGlyphType(int);
  void setTubeGlyphNumberOfSides(double);
  void setTubeGlyphRadius(double);  
  
protected slots:
  void updateWidgetFromMRMLDisplayNode();
  void updateWidgetFromMRMLDisplayPropertiesNode();

protected:
  QScopedPointer<qSlicerSpatialObjectsGlyphWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSpatialObjectsGlyphWidget);
  Q_DISABLE_COPY(qSlicerSpatialObjectsGlyphWidget);

  int updating;
};

#endif
