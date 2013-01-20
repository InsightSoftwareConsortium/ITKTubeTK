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

// Qt includes
#include <QApplication>
#include <QTimer>

// SpatialObjects includes
#include "qSlicerSpatialObjectsGlyphWidget.h"

// MRML includes
#include <vtkMRMLSpatialObjectsGlyphDisplayNode.h>
#include <vtkMRMLSpatialObjectsDisplayPropertiesNode.h>
#include <vtkMRMLScene.h>

// VTK includes
#include <vtkNew.h>

//-----------------------------------------------------------------------------
int qSlicerSpatialObjectsGlyphWidgetTest1( int argc, char * argv[] )
{
  QApplication app(argc, argv);

  vtkNew<vtkMRMLScene> scene;

  vtkNew<vtkMRMLSpatialObjectsGlyphDisplayNode> soDisplay;
  scene->AddNode(soDisplay.GetPointer());
  vtkNew<vtkMRMLSpatialObjectsDisplayPropertiesNode> soDisplayProperties;
  scene->AddNode(soDisplayProperties.GetPointer());

  qSlicerSpatialObjectsGlyphWidget widget;
  widget.setSpatialObjectsDisplayNode(static_cast<vtkMRMLNode*>(0));
  widget.setSpatialObjectsDisplayNode(soDisplay.GetPointer());

  soDisplay->SetAndObserveSpatialObjectsDisplayPropertiesNodeID(
      soDisplayProperties->GetID());

  widget.setSpatialObjectsDisplayNode(soDisplay.GetPointer());
  widget.show();

  if (argc < 2 || QString(argv[1]) != "-I")
    {
    QTimer::singleShot(200, &app, SLOT(quit()));
    }

  return app.exec();
}
