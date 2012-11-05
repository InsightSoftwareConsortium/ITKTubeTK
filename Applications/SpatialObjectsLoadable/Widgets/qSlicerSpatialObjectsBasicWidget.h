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

#ifndef __qSlicerSpatialObjectsBasicWidget_h
#define __qSlicerSpatialObjectsBasicWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

// qMRML includes
#include "qSlicerSpatialObjectsModuleWidgetsExport.h"

class qSlicerSpatialObjectsBasicWidgetPrivate;
class vtkMRMLNode;
class vtkMRMLSpatialObjectsNode;

class Q_SLICER_MODULE_SPATIALOBJECTS_WIDGETS_EXPORT
qSlicerSpatialObjectsBasicWidget : public qSlicerWidget
{
  Q_OBJECT
  QVTK_OBJECT

public:
  typedef qSlicerWidget Superclass;
  qSlicerSpatialObjectsBasicWidget(QWidget *parent=0);
  virtual ~qSlicerSpatialObjectsBasicWidget();

  vtkMRMLSpatialObjectsNode* spatialObjectsNode() const;

public slots:
  void setSpatialObjectsNode(vtkMRMLNode*);
  void setSpatialObjectsNode(vtkMRMLSpatialObjectsNode*);
  void setLineVisibility(int);
  void setTubeVisibility(int);
  void setGlyphVisibility(int);

protected slots:
  void updateWidgetFromMRML();

protected:
  QScopedPointer<qSlicerSpatialObjectsBasicWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSpatialObjectsBasicWidget);
  Q_DISABLE_COPY(qSlicerSpatialObjectsBasicWidget);

  int updating;
};

#endif
