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

#ifndef __qSlicerSpatialObjectsWidget_h
#define __qSlicerSpatialObjectsWidget_h

// Qt includes
#include <QWidget>

// CTK includes
#include <ctkVTKObject.h>

// SlicerQt includes
#include <qSlicerWidget.h>

// qMRML includes
#include "qSlicerSpatialObjectsModuleWidgetsExport.h"

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
  qSlicerSpatialObjectsWidget(QWidget *parent=0);
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
  void setSpatialObjectsNode(vtkMRMLNode *node);
  void setSpatialObjectsNode(vtkMRMLSpatialObjectsNode *node);

  void setSpatialObjectsDisplayNode(vtkMRMLNode *node);
  void setSpatialObjectsDisplayNode(vtkMRMLSpatialObjectsDisplayNode *node);

  void setVisibility(bool);
  void setColorByScalar();
  void onColorByScalarChanged(int);
  void onColorByScalarRangeChanged(double,double);
  void setColorByCellScalars();
  void setColorBySolid();
  void clickColorBySolid(bool);
  void onColorBySolidChanged(const QColor&);
  void setColorByCellScalarsColorTable(vtkMRMLNode*);
  void setOpacity(double);
  
  /// Set the values on the display node
  void setColor(const QColor&);
  void setAmbient(double);
  void setDiffuse(double);
  void setSpecular(double);
  void setSpecularPower(double);
  void setBackfaceCulling(bool);
  
protected slots:
  void updateWidgetFromMRML();

protected:
  QScopedPointer<qSlicerSpatialObjectsWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSpatialObjectsWidget);
  Q_DISABLE_COPY(qSlicerSpatialObjectsWidget);

  int inProcess;
};

#endif
