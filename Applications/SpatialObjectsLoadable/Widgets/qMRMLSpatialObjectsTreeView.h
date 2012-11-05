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

#ifndef __qMRMLSpatialObjectsTreeView_h
#define __qMRMLSpatialObjectsTreeView_h

// Qt includes
#include <QTreeView>
#include "qSlicerSpatialObjectsModuleWidget.h"

// CTK includes
#include <ctkPimpl.h>

// qMRML includes
#include "qMRMLTreeView.h"

#include "qSlicerSpatialObjectsModuleWidgetsExport.h"

class qMRMLSpatialObjectsTreeViewPrivate;
class vtkMRMLNode;
class vtkMRMLScene;
class vtkSlicerSpatialObjectsLogic;

/// \ingroup Slicer_QtModules_SpatialObjects
class Q_SLICER_MODULE_SPATIALOBJECTS_WIDGETS_EXPORT
qMRMLSpatialObjectsTreeView : public qMRMLTreeView
{
Q_OBJECT

public:
  typedef qMRMLTreeView Superclass;
  qMRMLSpatialObjectsTreeView(QWidget *parent=0);
  virtual ~qMRMLSpatialObjectsTreeView();

  virtual void setMRMLScene(vtkMRMLScene* scene);

  // Register the logic
  void setLogic(vtkSlicerSpatialObjectsLogic* logic);
  virtual bool clickDecoration(const QModelIndex& index);

signals:
  void onPropertyEditButtonClicked(QString id);
  // type 0-line 1-tube 2-glyph
  void visibilityChanged(int type);

protected:
  QScopedPointer<qMRMLSpatialObjectsTreeViewPrivate> d_ptr;
  #ifndef QT_NO_CURSOR
    bool viewportEvent(QEvent* e);
  #endif

private:
  Q_DECLARE_PRIVATE(qMRMLSpatialObjectsTreeView);
  Q_DISABLE_COPY(qMRMLSpatialObjectsTreeView);

  vtkSlicerSpatialObjectsLogic* Logic;
  
  // Toggle the visibility of an SpatialObjects.
  void onVisibilityColumnClicked(vtkMRMLNode* node);
};

#endif
