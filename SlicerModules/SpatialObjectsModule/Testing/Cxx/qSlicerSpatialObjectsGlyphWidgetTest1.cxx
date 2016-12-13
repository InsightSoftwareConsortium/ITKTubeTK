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

// Qt includes
#include <QApplication>
#include <QTimer>

// SpatialObjects includes
#include "qSlicerSpatialObjectsGlyphWidget.h"

// MRML includes
#include <vtkMRMLSpatialObjectsNode.h>
#include <vtkMRMLSpatialObjectsGlyphDisplayNode.h>
#include <vtkMRMLSpatialObjectsDisplayPropertiesNode.h>
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkNew.h>

//-----------------------------------------------------------------------------
int qSlicerSpatialObjectsGlyphWidgetTest1( int argc, char * argv[] )
{
  QApplication app( argc, argv );

  vtkNew<vtkMRMLScene> scene;

  vtkNew<vtkMRMLSpatialObjectsNode> so;
  scene->AddNode( so.GetPointer() );
  vtkNew<vtkMRMLSpatialObjectsGlyphDisplayNode> soDisplay;
  scene->AddNode( soDisplay.GetPointer() );
  vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> soDisplayProperties;
  scene->AddNode( soDisplayProperties.GetPointer() );

  qSlicerSpatialObjectsGlyphWidget widget;
  widget.setSpatialObjectsDisplayNode( so->GetDisplayNode() );

  so->SetAndObserveDisplayNodeID( soDisplay->GetID() );

  vtkMRMLNode* nullNode = NULL;
  widget.setSpatialObjectsDisplayNode( nullNode );
  widget.setSpatialObjectsDisplayNode( so->GetDisplayNode() );

  soDisplay->SetAndObserveSpatialObjectsDisplayPropertiesNodeID( 
      soDisplayProperties->GetID() );

  widget.setSpatialObjectsDisplayNode( nullNode );
  widget.setSpatialObjectsDisplayNode( so->GetDisplayNode() );
  widget.show();

  if( argc < 2 || QString( argv[1] ) != "-I" )
    {
    QTimer::singleShot( 200, &app, SLOT( quit() ) );
    }

  return app.exec();
}
