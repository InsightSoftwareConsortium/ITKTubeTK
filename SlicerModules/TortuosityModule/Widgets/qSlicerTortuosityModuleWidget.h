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

#ifndef __qSlicerTortuosityModuleWidget_h
#define __qSlicerTortuosityModuleWidget_h

// CTK includes
#include <ctkPimpl.h>

// SlicerQt includes
#include <qSlicerAbstractModuleWidget.h>
#include <qSlicerTortuosityModuleWidgetsExport.h>

class qSlicerTortuosityModuleWidgetPrivate;

class vtkMRMLNode;
class vtkMRMLSpatialObjectsNode;
class vtkSlicerTortuosityLogic;

/// \ingroup Slicer_QtModules_Tortuosity
class Q_SLICER_MODULE_TORTUOSITY_WIDGETS_EXPORT
qSlicerTortuosityModuleWidget : public qSlicerAbstractModuleWidget
{
  Q_OBJECT

public:
  typedef qSlicerAbstractModuleWidget Superclass;
  qSlicerTortuosityModuleWidget( QWidget *parent=0 );
  virtual ~qSlicerTortuosityModuleWidget();

public slots:
  void runMetrics( int flag );
  void setCurrentSpatialObjectsNode( vtkMRMLNode* node );
  void setCurrentSpatialObjectsNode( vtkMRMLSpatialObjectsNode* node );

protected slots:
  void loadColorsFromCSV();
  void runSelectedMetrics( bool run );
  void saveCurrentSpatialObjectAsCSV( bool save );
  void smoothingMethodChanged( int index );

protected:
  virtual void setup();

  QScopedPointer< qSlicerTortuosityModuleWidgetPrivate > d_ptr;

private:
  Q_DECLARE_PRIVATE( qSlicerTortuosityModuleWidget );
  Q_DISABLE_COPY( qSlicerTortuosityModuleWidget );
};

#endif
