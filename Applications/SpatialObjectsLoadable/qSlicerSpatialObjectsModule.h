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

#ifndef __qSlicerSpatialObjectsModule_h
#define __qSlicerSpatialObjectsModule_h

// SlicerQt includes
#include "qSlicerLoadableModule.h"

#include "qSlicerSpatialObjectsModuleExport.h"

class qSlicerAbstractModuleWidget;
class qSlicerSpatialObjectsModulePrivate;

/// \ingroup Slicer_QtModules_SpatialObjects
class Q_SLICER_QTMODULES_SPATIALOBJECTS_EXPORT qSlicerSpatialObjectsModule :
  public qSlicerLoadableModule
{
  Q_OBJECT
  Q_INTERFACES(qSlicerLoadableModule);

public:
  typedef qSlicerLoadableModule Superclass;
  qSlicerSpatialObjectsModule(QObject *parent=0);
  virtual ~qSlicerSpatialObjectsModule();

  virtual QString helpText() const;
  virtual QString acknowledgementText() const;
  virtual QStringList contributors() const;
  virtual QIcon icon() const;
  virtual QStringList categories() const;
  virtual QStringList dependencies() const;
  qSlicerGetTitleMacro(QTMODULE_TITLE);

protected:
  /// Initialize the module.
  /// Register the spatial objects reader/writer.
  virtual void setup();

  /// Create and return the widget representation associated to this module
  virtual qSlicerAbstractModuleRepresentation* createWidgetRepresentation();

  /// Create and return the logic associated to this module
  virtual vtkMRMLAbstractLogic* createLogic();

protected:
  QScopedPointer<qSlicerSpatialObjectsModulePrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerSpatialObjectsModule);
  Q_DISABLE_COPY(qSlicerSpatialObjectsModule);
};

#endif
